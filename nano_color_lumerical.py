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

#import math
#import sys
import os
import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from matplotlib import rcParams
from scipy.interpolate import interp1d
import re
#import nanocolor
import xlsxwriter
import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askdirectory
#from tkinter.filedialog import asksaveasfilename
import xarray as xr
from scipy.optimize import curve_fit

figureDPI = 72
savefigureFlag = True

# http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
XYZtolRGB = np.array(
    [
        [3.2406255, -1.537208, -0.4986286],
        [-0.9689307, 1.8757561, 0.0415175],
        [0.0557101, -0.2040211, 1.0569959],
    ]
)

# defines a single Gaussian function
def gaussian(x, mu, sigma,
             amp):  # this defines a function with the paramaters x=wavelength, mu=peak center, sigma=peak width, amp=peak height
    # Note that x data MUST be the first argument when defining a function that will be used with CurveFit
    y = (amp / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp((-(x - mu) ** 2) / (2 * sigma ** 2))
    return y

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

    else:
        HSV = cv2.cvtColor(RGB, cv2.COLOR_RGB2HSV)[0, 0, :]
        LAB = cv2.cvtColor(RGB, cv2.COLOR_RGB2LAB)[0, 0, :]
        RGB = RGB[0, 0, :]
        RGBg = RGBg[0, 0, :]
        rgb = np.zeros((RGB.shape))
        rgb[0] = RGB[0] / np.sum(RGB)
        rgb[1] = RGB[1] / np.sum(RGB)
        rgb[2] = RGB[2] / np.sum(RGB)
    HSV[0] = ShiftHOriginToValue(HSV[0], 360, 360.0 / 3, direction="ccw")
    HSV[0] = HSV[0] / 360.0
    LAB[0] = LAB[0] / 100.0
    LAB[1] = (LAB[1] + 128) / 255.0
    LAB[2] = (LAB[2] + 128) / 255.0
    RGBg = np.rint(RGBg * 255) / 255.0
    return RGB, HSV, LAB, XYZ, rgb, RGBg

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

def GetFilesToProcess():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    
    # the following line lets the user select a file that can be processed later
    #file_path = askopenfilename(filetypes=[('lumerical files', '*.txt'),('all files', '.*')])
    startpath = askdirectory()
    
    # the following block of code lets the user choose whether to process ONLY the file selected during the askopenfilename (Option 1)
    # or whether to also process all other files in that directory (option 2)
    # or whether to also process all other files in that directory AND all subdirectories (option 3)
    mode = input('root is : ' + startpath + '\n' + '\n' + 'Do you want to process one polarization pair in this folder (1), all files in this directory (2), or all files in this directory and all subdirectories (3)?')  # asks user which files to process
    filestoprocess = []  # initializes "files to process" as an empty list. We can then add files to process to this list based on option selected by user
    folders=set([])
    if mode == "1":
        file_path = askopenfilename(filetypes=[('lumerical files', '*.txt'),('all files', '.*')],initialdir=startpath)
        filestoprocess.append(file_path)  # only file processed will be the file selected during askopenfilename
        folders.add(os.path.basename(startpath))
    if mode == "2":
        folders.add(os.path.basename(startpath))
        for filename in os.listdir(startpath):  # gets all the filenames in the directory
            if filename[-3:] == "txt":  # if last 3 characters in filename are csv, then do the next line
                filestoprocess.append(os.path.join(startpath,filename))  # joins the directory to the filename to create a filepath, then adds that filepath to the "files to process" list
    if mode == "3":
        for (dirpath, dirnames, filenames) in os.walk(startpath):  # os.walk walks through the directory AND all subdirectories where a given file is stored; dirpath dirnames filenames lists all files in the directory and subdirectories
            for dirname in dirnames:
                folders.add(dirname)
            for filename in filenames:
                if filename[-3:] == "txt":  # if last 3 characters in filename are csv, then do the next line

                    filestoprocess.append(os.path.join(os.path.normpath(dirpath), filename))
    return(startpath,list(folders),filestoprocess)

# creating workbook
startpath,folders,filestoprocess=GetFilesToProcess()
workbook = xlsxwriter.Workbook(os.path.join(startpath, 'Lumerical_Data_Summary.xlsx'))

# adding worksheets with names
rawData = workbook.add_worksheet('Raw Data')
condensedData = workbook.add_worksheet('Condensed Data')
summaryData = workbook.add_worksheet('Summary Data')
sensitivityData = workbook.add_worksheet('Sensitivity')
dfLumerical=pd.DataFrame()

ris=set([])
polarizations=set([])
waves=set([])
for file in filestoprocess:
    fileLines = open(file, 'r')
    lines = fileLines.readlines()
    s1=0
    s2=0
    e1=0
    e2=0
    count = 0
    stage=0
    flgContent=False
    for line in lines:
        if (line=='\n') and (flgContent==True):
            stage=stage+1
            flgContent=False
            if stage==2:
                e1 = count-1
            elif stage==4:
                e2 = count-1
        elif (line!='\n') and (flgContent!=True):
            stage=stage+1
            flgContent=True
            if stage==1:
                s1 = count
            elif stage==3:
                s2 = count
        count=count+1
    fileLines.close()

    parentFolder=os.path.basename(os.path.dirname(file))
    #parentFolder=parentFolder.replace(" ", "_")
    ri=re.findall(r"[-+]?\d*\.\d+", file)[0]
    polarization=re.findall('\(\d*?\)', file)[0][1:-1]
    ris.add(ri)
    polarizations.add(polarization)
#file_open=file_dir+r"/"+file[:file.find("(")]+"("+polarization+").txt"
    dfFile=pd.read_csv(file, skiprows=s1, nrows=e1-s1)
    baseLabel=parentFolder+":"+"ri_"+ri+"_p_"+polarization+"_"
    dfLumerical[baseLabel+"l"]=dfFile["lambda"]*1e9
    waves.add(tuple(dfFile["lambda"]*1e9))
    dfLumerical[baseLabel+"A"]=dfFile[" Y"]
    dfFile=pd.read_csv(file, skiprows=s2, nrows=e2-s2)
    dfLumerical[baseLabel+"S"]=dfFile[" Y"]
    dfLumerical[baseLabel+"E"]=dfLumerical[baseLabel+"A"]+dfLumerical[baseLabel+"S"]
ris=list(ris)
polarizations=list(polarizations)
polarizations.append('m')
waves=list(waves)
# waveIncrement=np.ceil(np.max(np.abs(np.diff(waves))))
# for wave,index in zip(waves,range(len(waves))): 
#     waves[index]=np.array(wave)
#     smallestWaveIncrement=np.floor(np.min(np.abs(np.diff(waves[index]))))
#     if smallestWaveIncrement<waveIncrement:
#         waveIncrement=smallestWaveIncrement
waveIncrement=1.0
waveMin=360
waveMax=830
wavelengths = np.arange(waveMin, waveMax+waveIncrement, waveIncrement)
xrDataArray = xr.DataArray(
    dims=("folder", "ri", "polarization", "signalType", "signal"),
    coords={
        "folder" : folders,
        "ri": ris,
        "polarization": polarizations,
#        "wavelength": wavelengths,
        "signalType": ["A", "S", "E"],
        "signal": np.zeros(len(wavelengths)),
    },
)

for (columnName, columnData) in dfLumerical.iteritems():
    #print('Colunm Name : ', columnName)
    if columnName[-1:]!="l":
        folder=columnName[:columnName.find(":")]
        ri=columnName[columnName.find("ri_")+3:columnName.find("_p")]
        polarization=columnName[columnName.find("p_")+2:-2]
        signalType=columnName[-1:]
        baseLabel=folder+":"+"ri_"+ri+"_p_"+polarization+"_"   
        wavesIn=np.array(dfLumerical[baseLabel+"l"])
        dataIn=np.array(columnData.values)
        if np.min(wavesIn)>360:
            wavesIn=np.append(wavesIn,360)
            dataIn=np.append(dataIn,0)
        if np.max(wavesIn)<830:
            wavesIn=np.append(wavesIn,830)
            dataIn=np.append(dataIn,0)
        functionData=interp1d(wavesIn,dataIn, kind="cubic")
#        functionData=interp1d(wavesIn,dataIn, kind="linear")
        dataOut=functionData(wavelengths)
        xrDataArray.loc[dict(folder=folder,ri=ri,polarization=polarization,signalType=signalType)]=dataOut
xrDataArray.loc[dict(polarization='m')]=xrDataArray.mean(dim="polarization")    

colorDataFrame = pd.read_excel("data/all_1nm_data.xls", skiprows=63)
fCIEX = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 5], kind="cubic")
fCIEY = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 6], kind="cubic")
fCIEZ = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 7], kind="cubic")
fD65 = interp1d(colorDataFrame.values[:, 0], colorDataFrame.values[:, 2], kind="cubic")
CIEX = fCIEX(wavelengths)
CIEY = fCIEY(wavelengths)
CIEZ = fCIEZ(wavelengths)
D65 = fD65(wavelengths)
illum = D65
Yr = np.trapz(CIEY * illum, wavelengths)
nDataFrame = pd.read_csv("data/Johnson.csv", skiprows=0, nrows=49)
kDataFrame = pd.read_csv("data/Johnson.csv", skiprows=50, nrows=49)
fSplineN = interp1d(
    nDataFrame.values[:, 0] * 1000, nDataFrame.values[:, 1], kind="cubic"
)
fSplineK = interp1d(
    kDataFrame.values[:, 0] * 1000, kDataFrame.values[:, 1], kind="cubic"
)
m = fSplineN(wavelengths) + 1j * fSplineK(wavelengths)

channels=['Red','Green','Blue','Hue','Saturation','Value','CIE L*','CIE a*','CIE b*','X','Y','Z','r chromaticity','g chromaticity','b chromaticity','Red gamma','Green gamma','Blue gamma','lmax interpolated','lmax raw','lmax gauss fit']
xrDataSummary = xr.DataArray(
    dims=("folder", "ri", "channel"),
    coords={
        "folder" : folders,
        "ri": ris,
        "channel": channels,
    },
    )
for folder in folders:
    fig,ax=plt.subplots()
    for ri in ris:
        absorbance=xrDataArray.loc[dict(folder=folder,ri=ri,polarization='m',signalType='E')]
        scaleFactor=5e-14
        absorbance=absorbance/scaleFactor
        #RGB, HSV, LAB, XYZ, rgb, RGBg = absorbanceToTristim(wavelengths, absorbance, Yr, gammaFlag=True)
        tristims = absorbanceToTristim(wavelengths, absorbance, Yr, gammaFlag=True)
        channelIndex=0
        for cc in range(6):
            for pt in range (3):
                xrDataSummary.loc[dict(folder=folder,ri=ri,channel=channels[channelIndex])]=tristims[cc][pt]
                channelIndex=channelIndex+1
        baseLabel=folder+":"+"ri_"+ri   
        wavesIn=np.array(dfLumerical[baseLabel+"_p_0_l"])
        dataIn=np.mean([dfLumerical[baseLabel+"_p_0_E"],dfLumerical[baseLabel+"_p_90_E"]],axis=0)
        rangeBoolIntens = (dataIn >= 0.50 * np.max(dataIn)) 
        popt, pcov = curve_fit(gaussian, wavesIn[rangeBoolIntens], dataIn[rangeBoolIntens], p0=[wavesIn[np.argmax(dataIn)], 30, 87])
        xrDataSummary.loc[dict(folder=folder,ri=ri,channel='lmax interpolated')]=wavelengths[np.argmax(np.array(absorbance))]
        xrDataSummary.loc[dict(folder=folder,ri=ri,channel='lmax raw')]=wavesIn[np.argmax(dataIn)]
        xrDataSummary.loc[dict(folder=folder,ri=ri,channel='lmax gauss fit')]=popt[0]
        ax.plot(wavelengths,absorbance,color=tristims[5],label=ri)
    plt.legend()
dfSummary = xrDataSummary.to_dataframe('value').unstack()
dfSummary = dfSummary["value"]
dfSummary.reset_index(inplace=True)
dfSummary.to_excel(os.path.join(startpath, 'Lumerical_Data_Summary.xlsx'),sheet_name="ColorSummary")


activeColorSpectrum=['red','green','blue','cyan','grey','black','black','magenta','yellow','salmon','lightseagreen','skyblue','gold','mediumorchid','darkcyan','darkred','darkgreen','darkblue','black','black','black']

for folder in folders:
    fig,axes=plt.subplots(7,3,sharex=True, gridspec_kw=dict(wspace=0, hspace=0))
    channelIndex=0
    refIndex=np.array(ris)
    refIndex=refIndex.astype(float)

    for cc in range(7):
        for pt in range (3):
            dta=xrDataSummary.loc[dict(folder=folder,channel=channels[channelIndex])].values
            axes[cc,pt].plot(refIndex,dta,'o',color=activeColorSpectrum[channelIndex])
            axes[cc,pt].axes.get_yaxis().set_visible(False)
            channelIndex=channelIndex+1



#dfRaw = xrDataArray.loc[dict(polarization='m',signalType='E')].to_dataframe('value')
#dfRaw = dfRaw["value"]
#dfRaw.reset_index(inplace=True)

# rcParams["font.family"] = "sans-serif"
# rcParams["font.sans-serif"] = ["Arial"]

# fontSizeLarge = 14
# fontSizeLargest = 18
# fontSizeSmall = 12
# fontSizeLegend = 12

# figColor, cie = plt.subplots(1, 1, figsize=(6, 6))
# cie.plot(waves, CIEX, "r", label="CIE X")
# cie.plot(waves, CIEY, "g", label="CIE Y")
# cie.plot(waves, CIEZ, "b", label="CIE Z")
# cie.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
# cie.set_ylabel("Response", fontsize=fontSizeLarge, fontweight="bold")
# legend = cie.legend(loc="upper right", fontsize=fontSizeLegend)
# cie.set_xlim([360, 830])
# cie.set_ylim([0, 2])
# cie.set_yticks(np.array([0, 0.5, 1, 1.5, 2]))
# cie.set_yticklabels(np.array([0, 0.5, 1, 1.5, 2]), fontsize=fontSizeLarge)
# cie.set_xticks(np.array([400, 500, 600, 700, 800]))
# cie.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
# if savefigureFlag:
#     figColor.savefig("CIE_XYZ.png", dpi=figureDPI)

# figIll, ill = plt.subplots(1, 1, figsize=(6, 6))
# ill.plot(waves, illum, "k", label="D65 illuminant")
# ill.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
# ill.set_ylabel("Intensity", fontsize=fontSizeLarge, fontweight="bold")
# legend = ill.legend(loc="upper right", fontsize=fontSizeLegend)
# ill.set_xlim([360, 830])
# ill.set_xticks(np.array([400, 500, 600, 700, 800]))
# ill.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
# ill.set_ylim([0, 120])
# ill.set_yticks(np.array([0, 25, 50, 75, 100]))
# ill.set_yticklabels(np.array([0, 25, 50, 75, 100]), fontsize=fontSizeLarge)
# if savefigureFlag:
#     figIll.savefig("D65.png", dpi=figureDPI)

# figRI, opt = plt.subplots(1, 1, figsize=(6, 6))
# opt.plot(waves, fSplineN(waves), "-k", label="refractive index (n)")
# opt.plot(waves, fSplineK(waves), ":k", label="extinction coefficient (k)")
# opt.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
# opt.set_ylabel("Complex Refractive Index", fontsize=fontSizeLarge, fontweight="bold")
# legend = opt.legend(loc="upper left", fontsize=fontSizeLegend)
# opt.set_xlim([360, 830])
# opt.set_xticks(np.array([400, 500, 600, 700, 800]))
# opt.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
# opt.set_ylim([0, 5.2])
# opt.set_yticks(np.array([0, 1, 2, 3, 4, 5]))
# opt.set_yticklabels(np.array([0, 1, 2, 3, 4, 5]), fontsize=fontSizeLarge)
# if savefigureFlag:
#     figRI.savefig("RI.png", dpi=figureDPI)

# figExt, ext = plt.subplots(1, 1, figsize=(6, 6))
# ris = np.array([1.33])
# diameters = np.array([115, 100, 80, 60, 35])
# mie = "qext"
# for diameter in diameters:
#     for ri in ris:
#         if storeXarray:
#             signal = mieEfficenciesDataArray.loc[
#                 dict(ri=ri, diameter=diameter, mie=mie)
#             ].values
#         else:
#             signal = dfmieEfficencies[
#                 (dfmieEfficencies["ri"] == ri)
#                 & (dfmieEfficencies["diameter"] == diameter)
#             ][mie].values[0]
#         dfBool = (dfColorMono["diameter"] == diameter) & (dfColorMono["ri"] == ri)
#         color = dfColorMono[dfBool][["Rg", "Gg", "Bg"]].values
#         ext.plot(waves, signal, label=str(diameter) + " nm", color=color[0, :])
# ext.set_xlabel("Wavelength (nm)", fontsize=fontSizeLarge, fontweight="bold")
# ext.set_ylabel(r"Q$\mathbf{_{ext}}$", fontsize=fontSizeLarge, fontweight="bold")
# legend = ext.legend(loc="upper left", fontsize=fontSizeLegend)
# ext.set_xlim([360, 830])
# ext.set_xticks(np.array([400, 500, 600, 700, 800]))
# ext.set_xticklabels(np.array([400, 500, 600, 700, 800]), fontsize=fontSizeLarge)
# # ext.set_ylim([0, 0.5])
# # ext.set_yticks(np.array([0,.1,.2,.3,.4,.5]))
# # ext.set_yticklabels(np.array([0,.1,.2,.3,.4,.5]),fontsize = fontSizeLarge)
# if savefigureFlag:
#     figExt.savefig("Au_Ext.png", dpi=figureDPI)
