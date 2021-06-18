#! /usr/bin/env python
from __future__ import division

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 20:11:19 2018

@author: Kevin
"""

# excecute the following commands in the anaconda prompt
# conda install -c anaconda xarray
# conda install -c anaconda shapely
# pip install PyMieScatt
# pip install opencv-python

# required Excel file with CIE color information
# https://www.rit.edu/science/munsell-color-science-lab-educational-resources?ref=rit-search#useful-color-data
# Full set of 1nm data, including all of the following:
#    Illuminant A
#    Illuminant D65
#    VM(λ) 1988 Spectral Luminous Efficiency Function for photopic vision
#    V'(λ) Spectral Luminous Efficiency Function for scotopic vision
#    1931 2° CIE Standard Colorimetric Observer Data
#    1964 10 °CIE Standard Colorimetric Observer Data
#    Excel with all of the above
# save as 'all_1nm_data.xls'

# required Excel file with refractive index information
# https://refractiveindex.info/?shelf=main&book=Au&page=Johnson
# Data  [CSV - comma separated]
# save as 'Johnson.csv'

# P. B. Johnson and R. W. Christy. Optical constants of the noble metals,
#     Phys. Rev. B 6, 4370-4379 (1972)
# M. N. Polyanskiy, "Refractive index database," https://refractiveindex.info.
#     Accessed on 2019-09-30.

# http://www.color.org/srgb04.xalter
# https://pymiescatt.readthedocs.io/en/latest/index.html
# http://cvrl.ioo.ucl.ac.uk/index.htm

import math
import sys
import os
import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from scipy.interpolate import interp1d

import nanocolor

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

figureDPI = 72
savefigureFlag = True

if sys.version_info >= (3, 0):
    import xarray as xr

    storeXarray = True
else:
    storeXarray = False

# http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
XYZtolRGB = np.array(
    [
        [3.2406255, -1.537208, -0.4986286],
        [-0.9689307, 1.8757561, 0.0415175],
        [0.0557101, -0.2040211, 1.0569959],
    ]
)


def single_gauss_func(x, a, x0, sigma):
    gauss = np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    normal_gauss = a / np.sum(gauss) * gauss
    return normal_gauss


def absorbanceToTristim(waves, absorbance, Yr, gammaFlag=True):
    pixel = np.zeros((1, 1, 3), dtype=np.float32)
    pixel[0, 0, 0] = np.trapz(CIEX * illum * 10 ** -absorbance, waves) / Yr
    pixel[0, 0, 1] = np.trapz(CIEY * illum * 10 ** -absorbance, waves) / Yr
    pixel[0, 0, 2] = np.trapz(CIEZ * illum * 10 ** -absorbance, waves) / Yr
    XYZ = pixel[0, 0, :]
    RGB = cv2.cvtColor(pixel, cv2.COLOR_XYZ2RGB)
    RGBg = np.zeros((RGB.shape), dtype=np.float32)
    for cc in range(RGB.shape[2]):
        if RGB[0, 0, cc] <= 0.0031308:
            RGBg[0, 0, cc] = 12.92 * RGB[0, 0, cc]
        else:
            RGBg[0, 0, cc] = 1.055 * RGB[0, 0, cc] ** (1 / 2.4) - 0.055
        if RGBg[0, 0, cc] > 1:
            RGBg[0, 0, cc] = 1
        elif RGBg[0, 0, cc] < 0:
            RGBg[0, 0, cc] = 0
    if gammaFlag:
        HSV = cv2.cvtColor(RGBg, cv2.COLOR_RGB2HSV)[0, 0, :]
        LAB = cv2.cvtColor(RGBg, cv2.COLOR_RGB2LAB)[0, 0, :]
        RGB = RGB[0, 0, :]
        RGBg = RGBg[0, 0, :]
        rgb = np.zeros((RGBg.shape))
        rgb[0] = RGBg[0] / np.sum(RGBg)
        rgb[1] = RGBg[1] / np.sum(RGBg)
        rgb[2] = RGBg[2] / np.sum(RGBg)
        rat = np.zeros((RGBg.shape))
        rat[0] = RGBg[0] / RGBg[1]
        rat[1] = RGBg[0] / RGBg[2]
        rat[2] = RGBg[1] / RGBg[2]
    else:
        HSV = cv2.cvtColor(RGB, cv2.COLOR_RGB2HSV)[0, 0, :]
        LAB = cv2.cvtColor(RGB, cv2.COLOR_RGB2LAB)[0, 0, :]
        RGB = RGB[0, 0, :]
        RGBg = RGBg[0, 0, :]
        rgb = np.zeros((RGB.shape))
        rgb[0] = RGB[0] / np.sum(RGB)
        rgb[1] = RGB[1] / np.sum(RGB)
        rgb[2] = RGB[2] / np.sum(RGB)
        rat = np.zeros((RGB.shape))
        rat[0] = RGB[0] / RGB[1]
        rat[1] = RGB[0] / RGB[2]
        rat[2] = RGB[1] / RGB[2]
    HSV[0] = ShiftHOriginToValue(HSV[0], 360, 360.0 / 3, direction="ccw")
    HSV[0] = HSV[0] / 360.0
    LAB[0] = LAB[0] / 100.0
    LAB[1] = (LAB[1] + 128) / 255.0
    LAB[2] = (LAB[2] + 128) / 255.0
    RGBg = np.rint(RGBg * 255) / 255.0
    return RGB, HSV, LAB, XYZ, rgb, rat, RGBg


def absorbanceToRGB(waves, absorbance, Yr, gammaFlag=True):
    XYZ = np.zeros((3), dtype=np.float32)
    XYZ[0] = np.trapz(CIEX * illum * 10 ** -absorbance, waves) / Yr
    XYZ[1] = np.trapz(CIEY * illum * 10 ** -absorbance, waves) / Yr
    XYZ[2] = np.trapz(CIEZ * illum * 10 ** -absorbance, waves) / Yr
    RGB = np.matmul(XYZtolRGB, XYZ)
    RGBg = np.zeros((RGB.shape), dtype=np.float32)
    for cc in range(3):
        if RGB[cc] <= 0.0031308:
            RGBg[cc] = 12.92 * RGB[cc]
        else:
            RGBg[cc] = 1.055 * RGB[cc] ** (1 / 2.4) - 0.055
        if RGBg[cc] > 1:
            RGBg[cc] = 1
        elif RGBg[cc] < 0:
            RGBg[cc] = 0
    rgb = np.zeros((RGB.shape))
    if gammaFlag:
        rgb[0] = RGBg[0] / np.sum(RGBg)
        rgb[1] = RGBg[1] / np.sum(RGBg)
        rgb[2] = RGBg[2] / np.sum(RGBg)
    else:
        rgb[0] = RGB[0] / np.sum(RGB)
        rgb[1] = RGB[1] / np.sum(RGB)
        rgb[2] = RGB[2] / np.sum(RGB)
    RGBg = np.rint(RGBg * 255) / 255.0
    return RGB, XYZ, rgb, RGBg


def ShiftHOriginToValue(hue, maxHue, newOrigin, direction="cw"):
    shifthsv = np.copy(hue).astype("float")
    shiftAmount = maxHue - newOrigin
    shifthsv[hue < newOrigin] = shifthsv[hue < newOrigin] + shiftAmount
    shifthsv[hue >= newOrigin] = shifthsv[hue >= newOrigin] - newOrigin
    hue = shifthsv
    if direction == "ccw":
        hue = maxHue - hue
    return hue


waves = np.arange(360, 830, 0.1)

colorDataFrame = pd.read_excel("data/all_1nm_data.xls", skiprows=63)
fCIEX = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 5], kind="cubic")
fCIEY = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 6], kind="cubic")
fCIEZ = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 7], kind="cubic")
fD65 = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 2], kind="cubic")
CIEX = fCIEX(waves)
CIEY = fCIEY(waves)
CIEZ = fCIEZ(waves)
D65 = fD65(waves)
illum = D65
Yr = np.trapz(CIEY * illum, waves)

nDataFrame = pd.read_csv("data/Johnson.csv", skiprows=0, nrows=49)
kDataFrame = pd.read_csv("data/Johnson.csv", skiprows=50, nrows=49)
fSplineN = interp1d(
    nDataFrame.values[:, 0] * 1000, nDataFrame.values[:, 1], kind="cubic"
)
fSplineK = interp1d(
    kDataFrame.values[:, 0] * 1000, kDataFrame.values[:, 1], kind="cubic"
)
m = fSplineN(waves) + 1j * fSplineK(waves)

root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
#file_path = askopenfilename(initialdir=filePathImage,filetypes=[('image files', '*.jpg | *.jpeg | *.png'),('video files', '*.mp4 | *.mkv | *.avi | *.MOV'),('all files', '.*')])
file_path = askopenfilename(filetypes=[('lumerical files', '*.txt'),('all files', '.*')])
file_pathSplit = os.path.split(file_path)
file_dir=file_pathSplit[0]
file_file=file_pathSplit[1]
dfLumerical=pd.read_csv(file_path, nrows=100)



# ris=np.arange(1,2,0.005)
# diameters=np.arange(30,120,1)
ris = np.array([1.33, 1.38, 1.43])
diameters = np.array([115, 100, 80, 60, 35])
diameter_stdevs = diameters * 0.1

# gen='G3'
# diameter=34.1
# diameter_stdev=6.7

# gen='G5'
# diameter=59.8
# diameter_stdev=7.7

# gen='G7'
# diameter=81.5
# diameter_stdev=7.6

# gen='G9'
# diameter=115
# diameter_stdev=17

scaleFactorMono = 16
scaleFactorPoly = 0.4
colorDataMono = []
colorDataPoly = []
bextArray = np.zeros((waves.shape[0]))
sizeDistributionDiameterBins = np.arange(20.0, 200.0, 1)
if storeXarray:
    mieEfficenciesDataArray = xr.DataArray(
        dims=("ri", "diameter", "wavelength", "mie"),
        coords={
            "ri": ris,
            "diameter": diameters,
            "wavelength": waves,
            "mie": ["qext", "qsca", "qabs", "g", "qpr", "qback", "qratio"],
        },
    )
else:
    mieEfficencies = []
spectCount = 0
for diameter, _ in zip(diameters, diameter_stdevs):
    for ri in ris:
        print("diameter = " + str(diameter) + " ri = " + str(ri))
        qbk = np.empty(len(waves), dtype=np.float64)
        qpr = np.empty(len(waves), dtype=np.float64)
        albedo = np.empty(len(waves), dtype=np.float64)
        g = np.empty(len(waves), dtype=np.float64)
        qext, qsca, qabs = nanocolor.mie_efficiencies_by_wavelength(
            diameter,
            np.asarray([1.0]),
            waves,
            m.reshape((1, -1)),
            medium=ri,
        )
        lmax = waves[(waves > 450) & (waves < 800)][
            np.argmax(qext[(waves > 450) & (waves < 800)])
        ]
        RGB, HSV, LAB, XYZ, rgb, rat, RGBg = absorbanceToTristim(
            waves, qext / scaleFactorMono, Yr, gammaFlag=True
        )
        # RGB, XYZ, rgb, RGBg = absorbanceToRGB(
        #     waves, qext / scaleFactorMono, Yr, gammaFlag=True
        # )
        colorDataMono.append(
            {
                "lmax": lmax,
                "ri": ri,
                "diameter": diameter,
                "R": RGB[0],
                "G": RGB[1],
                "B": RGB[2],
                "H": HSV[0],
                "S": HSV[1],
                "V": HSV[2],
                "L*": LAB[0],
                "a*": LAB[1],
                "b*": LAB[2],
                "r": rgb[0],
                "g": rgb[1],
                "b": rgb[2],
                "R/G": rat[0],
                "R/B": rat[1],
                "G/B": rat[2],
                "Rg": RGBg[0],
                "Gg": RGBg[1],
                "Bg": RGBg[2],
            }
        )
        # colorDataMono.append(
        #     {
        #         "ri": ri,
        #         "diameter": diameter,
        #         "R": RGB[0],
        #         "G": RGB[1],
        #         "B": RGB[2],
        #         "r": rgb[0],
        #         "g": rgb[1],
        #         "b": rgb[2],
        #         "Rg": RGBg[0],
        #         "Gg": RGBg[1],
        #         "Bg": RGBg[2],
        #     }
        # )
        if storeXarray:
            mieEfficenciesDataArray.loc[
                dict(ri=ri, diameter=diameter, mie="qext")
            ] = qext
            mieEfficenciesDataArray.loc[
                dict(ri=ri, diameter=diameter, mie="qsca")
            ] = qsca
            mieEfficenciesDataArray.loc[
                dict(ri=ri, diameter=diameter, mie="qabs")
            ] = qabs
            mieEfficenciesDataArray.loc[dict(ri=ri, diameter=diameter, mie="g")] = g
            mieEfficenciesDataArray.loc[dict(ri=ri, diameter=diameter, mie="qpr")] = qpr
            # mieEfficenciesDataArray.loc[
            #     dict(ri=ri, diameter=diameter, mie="qback")
            # ] = qback
            # mieEfficenciesDataArray.loc[
            #     dict(ri=ri, diameter=diameter, mie="qratio")
            # ] = qratio
        else:
            mieEfficencies.append(
                {
                    "ri": ri,
                    "diameter": diameter,
                    "qext": qext,
                    "qsca": qsca,
                    "qabs": qabs,
                }
            )

        # sizeDistribution = single_gauss_func(
        #     sizeDistributionDiameterBins, 10, diameter, diameter_stdev
        # )
        # for index in range(waves.shape[0]):
        #     bext, bsca, babs, G, bpr, bback, bratio = ps.Mie_SD(
        #         m[index],
        #         waves[index],
        #         sizeDistributionDiameterBins,
        #         sizeDistribution,
        #         nMedium=ri,
        #         asDict=False,
        #     )
        #     bextArray[index] = bext
        # RGB, HSV, LAB, XYZ, rgb, rat, RGBg = absorbanceToTristim(
        #     waves, bextArray / scaleFactorPoly, Yr, gammaFlag=True
        # )
        # lmax = waves[(waves > 450) & (waves < 800)][
        #     np.argmax(qext[(waves > 450) & (waves < 800)])
        # ]
        # colorDataPoly.append(
        #     {
        #         "lmax": lmax,
        #         "ri": ri,
        #         "diameter": diameter,
        #         "R": RGB[0],
        #         "G": RGB[1],
        #         "B": RGB[2],
        #         "H": HSV[0],
        #         "S": HSV[1],
        #         "V": HSV[2],
        #         "L*": LAB[0],
        #         "a*": LAB[1],
        #         "b*": LAB[2],
        #         "r": rgb[0],
        #         "g": rgb[1],
        #         "b": rgb[2],
        #         "R/G": rat[0],
        #         "R/B": rat[1],
        #         "G/B": rat[2],
        #         "Rg": RGBg[0],
        #         "Gg": RGBg[1],
        #         "Bg": RGBg[2],
        #     }
        # )

dfColorMono = pd.DataFrame(colorDataMono)
dfColorMono = dfColorMono[
    [
        "ri",
        "diameter",
        "R",
        "G",
        "B",
        "H",
        "S",
        "V",
        "L*",
        "a*",
        "b*",
        "r",
        "g",
        "b",
        "R/G",
        "R/B",
        "G/B",
        "Rg",
        "Gg",
        "Bg",
        "lmax",
    ]
]
# dfColorMono=dfColorMono[['ri','diameter','R','G','B','r','g','b','Rg','Gg','Bg']]
if not (storeXarray):
    dfmieEfficencies = pd.DataFrame(mieEfficencies)

# now = datetime.datetime.now()
# writer = pd.ExcelWriter("ModelingSpectra" + now.strftime("%m-%d-%Yat%H-%M") + ".xlsx")
# workbook = writer.book
# dfColorMono.to_excel(writer, "ModelingResultsMono", index=False)
# writer.save()
# if len(colorDataPoly) > 0:
#     dfColorPoly = pd.DataFrame(colorDataPoly)
#     dfColorPoly = dfColorPoly[
#         [
#             "ri",
#             "diameter",
#             "R",
#             "G",
#             "B",
#             "H",
#             "S",
#             "V",
#             "L*",
#             "a*",
#             "b*",
#             "r",
#             "g",
#             "b",
#             "R/G",
#             "R/B",
#             "G/B",
#             "Rg",
#             "Gg",
#             "Bg",
#             "lmax",
#         ]
#     ]

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial"]

fontSizeLarge = 14
fontSizeLargest = 18
fontSizeSmall = 12
fontSizeLegend = 12

figColor, cie = plt.subplots(1, 1, figsize=(6, 6))
cie.plot(waves, CIEX, "r", label="CIE X")
cie.plot(waves, CIEY, "g", label="CIE Y")
cie.plot(waves, CIEZ, "b", label="CIE Z")
cie.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
cie.set_ylabel("Response", fontsize=fontSizeLarge, fontweight="bold")
legend = cie.legend(loc="upper right", fontsize=fontSizeLegend)
cie.set_xlim([360, 830])
cie.set_ylim([0, 2])
cie.set_yticks(np.array([0, 0.5, 1, 1.5, 2]))
cie.set_yticklabels(np.array([0, 0.5, 1, 1.5, 2]), fontsize=fontSizeLarge)
cie.set_xticks(np.array([400, 500, 600, 700, 800]))
cie.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
if savefigureFlag:
    figColor.savefig("CIE_XYZ.png", dpi=figureDPI)

figIll, ill = plt.subplots(1, 1, figsize=(6, 6))
ill.plot(waves, illum, "k", label="D65 illuminant")
ill.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
ill.set_ylabel("Intensity", fontsize=fontSizeLarge, fontweight="bold")
legend = ill.legend(loc="upper right", fontsize=fontSizeLegend)
ill.set_xlim([360, 830])
ill.set_xticks(np.array([400, 500, 600, 700, 800]))
ill.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
ill.set_ylim([0, 120])
ill.set_yticks(np.array([0, 25, 50, 75, 100]))
ill.set_yticklabels(np.array([0, 25, 50, 75, 100]), fontsize=fontSizeLarge)
if savefigureFlag:
    figIll.savefig("D65.png", dpi=figureDPI)

figRI, opt = plt.subplots(1, 1, figsize=(6, 6))
opt.plot(waves, fSplineN(waves), "-k", label="refractive index (n)")
opt.plot(waves, fSplineK(waves), ":k", label="extinction coefficient (k)")
opt.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
opt.set_ylabel("Complex Refractive Index", fontsize=fontSizeLarge, fontweight="bold")
legend = opt.legend(loc="upper left", fontsize=fontSizeLegend)
opt.set_xlim([360, 830])
opt.set_xticks(np.array([400, 500, 600, 700, 800]))
opt.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
opt.set_ylim([0, 5.2])
opt.set_yticks(np.array([0, 1, 2, 3, 4, 5]))
opt.set_yticklabels(np.array([0, 1, 2, 3, 4, 5]), fontsize=fontSizeLarge)
if savefigureFlag:
    figRI.savefig("RI.png", dpi=figureDPI)

figExt, ext = plt.subplots(1, 1, figsize=(6, 6))
ris = np.array([1.33])
diameters = np.array([115, 100, 80, 60, 35])
mie = "qext"
for diameter in diameters:
    for ri in ris:
        if storeXarray:
            signal = mieEfficenciesDataArray.loc[
                dict(ri=ri, diameter=diameter, mie=mie)
            ].values
        else:
            signal = dfmieEfficencies[
                (dfmieEfficencies["ri"] == ri)
                & (dfmieEfficencies["diameter"] == diameter)
            ][mie].values[0]
        dfBool = (dfColorMono["diameter"] == diameter) & (dfColorMono["ri"] == ri)
        color = dfColorMono[dfBool][["Rg", "Gg", "Bg"]].values
        ext.plot(waves, signal, label=str(diameter) + " nm", color=color[0, :])
ext.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
ext.set_ylabel(r"Q$\mathbf{_{ext}}$", fontsize=fontSizeLarge, fontweight="bold")
legend = ext.legend(loc="upper left", fontsize=fontSizeLegend)
ext.set_xlim([360, 830])
ext.set_xticks(np.array([400, 500, 600, 700, 800]))
ext.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
# ext.set_ylim([0, 0.5])
# ext.set_yticks(np.array([0,.1,.2,.3,.4,.5]))
# ext.set_yticklabels(np.array([0,.1,.2,.3,.4,.5]),fontsize = fontSizeLarge)
if savefigureFlag:
    figExt.savefig("Au_Ext.png", dpi=figureDPI)
