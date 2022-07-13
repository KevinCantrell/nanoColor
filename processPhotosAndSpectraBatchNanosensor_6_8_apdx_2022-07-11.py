#this script an process directories of data containing both photos and spectra
#the main folder should include the sample identifiers (type of particle, sensor chip replicate, and date)
#last 10 characters need to be the date
#e.g. Main folder name = 600uL_A_2022-06-28
#within this folder, make folders called "Spectra" and "Photos"
#within Spectra folder, first 3 characters of filename need to be the RI medium (e.g., "air")
#within Photos folder, make separate folders where the folder name is the RI medium

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



from __future__ import division
import cv2
import numpy as np
import time
try:
    import Tkinter as tk
    from tkFileDialog import askopenfilename
    from tkFileDialog import asksaveasfilename
except ImportError:
    import tkinter as tk
    from tkinter.filedialog import askopenfilename
    from tkinter.filedialog import asksaveasfilename

try:
    input=raw_input
except NameError:
    pass
import math
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.spatial.distance import cdist 
import matplotlib.pyplot as plt
import os
import datetime
import scipy.signal
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

riValues_old =           [1.0000, 1.3310, 1.3436, 1.3612, 1.3774, 1.3972, 1.4176]
riValues_G7Au_C_062719 = [1.0000, 1.3310, 1.3436, 1.3612, 1.3774, 1.3972, 1.4188]
riValues_G7Au_D_062719 = [1.0000, 1.3310, 1.3436, 1.3612, 1.3774, 1.3972, 1.4188]
riValues_G9Au_C_062719 = [1.0000, 1.3310, 1.3436, 1.3612, 1.3774, 1.3972, 1.4188]
riValues_062819 =        [1.0000, 1.3299, 1.3458, 1.3615, 1.3786, 1.3885, 1.4188]
riValues_070119 =        [1.0000, 1.3299, 1.3466, 1.3616, 1.3799, 1.3886, 1.4188]
riValues_070319 =        [1.0000, 1.3310, 1.3465, 1.3615, 1.3801, 1.3885, 1.4190]
riValues_071019 =        [1.0000, 1.3327, 1.3459, 1.3615, 1.3799, 1.3882, 1.4191]
riValues_071219 =        [1.0000, 1.3327, 1.3459, 1.3615, 1.3799, 1.3994, 1.4190]
riValues_071519 =        [1.0000, 1.3327, 1.3481, 1.3701, 1.3779, 1.3994, 1.4190]
riValues_071619 =         [1.0000, 1.3327, 1.3470, 1.3713, 1.3796, 1.3998, 1.4215]
riValues_100819 =        [1.0000, 1.3328, 1.3471, 1.3637, 1.3819, 1.4000, 1.4208]
riValues_112319 =        [1.0000, 1.3328, 1.3486, 1.3637, 1.3813, 1.3998, 1.4199]
riValues_220628 = [1.0000, 1.3325, 1.3480, 1.3617, 1.3735, 1.3844, 1.3954]
riValues_220701A = [1.0000, 1.3328, 1.3473, 1.3593, 1.3708, 1.3805, 1.3907]
riValues_220701D = [1.0000, 1.3329, 1.3475, 1.3597, 1.3716, 1.3821, 1.3892]
riValues_220705 = [1.0000, 1.3326, 1.3469, 1.3595, 1.3710, 1.3811, 1.3890]
riValues_220706 = [1.0000, 1.3329, 1.3478, 1.3587, 1.3710, 1.3808, 1.3899]


#average of G9 G-I
riValues_070519 =        [1.0000, 1.3313, 1.3464, 1.3621, 1.3805, 1.3899, 1.4206]

    
lower_lim = np.array([0,0.06,0],dtype=np.float32)
#lower_lim = np.array([0,0.04,0],dtype=np.float32)
upper_lim = np.array([360,1,1],dtype=np.float32)
lower_lim_grey = np.array([0,0,0],dtype=np.float32)
upper_lim_grey = np.array([360,0.05,1],dtype=np.float32)
scalemin=0.65
MeanWindow=1
FrameByFrameToggle=False
ColorRectOver = True
ColorRectangle = False
ColorRect = (0,0,1,1)      
GreyRectOver = False
GreyRectangle = False
rebalanceToggle=True
localRebalanceToggle=True
ListDataToggle=False
RecordFlag=False
GreyRect = (0,0,1,1) 
times=[]
clocks=[]
hues=[]
huesr=[]
HSVMeans=[]
RGBMeans=[]
LABMeans=[]
RGBGreys=[]
reds=[]
greens=[]
blues=[]
LABa=[]
LABb=[]
bWB=[]
gWB=[]
rWB=[]
dataPoint=0
interlaceStart=0
RectList=[]
#RectList=[(211, 149, 43, 132), (274, 151, 49, 129), (344, 153, 46, 127)]
#RectList=[(224, 154, 59, 153), (283, 161, 55, 145), (338, 149, 59, 155)]
#RectList=[(214, 135, 52, 153), (270, 140, 63, 151), (338, 139, 65, 148)]
#ColorRectOver=True
GreyRectOver = False
font = cv2.FONT_HERSHEY_SIMPLEX
#fourcc = cv2.VideoWriter_fourcc(*'XVID')
file_base='C:/Users/Kevin/Dropbox/Whatman/'
#RectList=[(212, 151, 75, 163)]

#defines a single Gaussian function
def gaussian(x,mu,sigma,amp):#this defines a function with the paramaters x=wavelength, mu=peak center, sigma=peak width, amp=peak height
        #Note that x data MUST be the first argument when defining a function that will be used with CurveFit
        y=(amp/np.sqrt(2*np.pi*sigma**2))*np.exp((-(x-mu)**2)/(2*sigma**2))
        return y
    
def Cross(x, y, crossPoint=0, direction='cross'):
    """
    Given a Series returns all the index values where the data values equal 
    the 'cross' value. 

    Direction can be 'rising' (for rising edge), 'falling' (for only falling 
    edge), or 'cross' for both edges
    """
    # Find if values are above or bellow yvalue crossing:
    above=y > crossPoint
    below=np.logical_not(above)
    left_shifted_above = above[1:]
    left_shifted_below = below[1:]
    x_crossings = []
    # Find indexes on left side of crossing point
    if direction == 'rising':
        idxs = (left_shifted_above & below[0:-1]).nonzero()[0]
    elif direction == 'falling':
        idxs = (left_shifted_below & above[0:-1]).nonzero()[0]
    else:
        rising = left_shifted_above & below[0:-1]
        falling = left_shifted_below & above[0:-1]
        idxs = (rising | falling).nonzero()[0]

    # Calculate x crossings with interpolation using formula for a line:
    x1 = x[idxs]
    x2 = x[idxs+1]
    y1 = y[idxs]
    y2 = y[idxs+1]
    x_crossings = (crossPoint-y1)*(x2-x1)/(y2-y1) + x1

    return x_crossings,idxs

def PeakFind(SegmentX,SegmentY,peakHeightThreshold=0):
    """
    Finds peaks in a data array based on zero crossing in the derivative (calculated with numpy's gradient function)
    
    see https://docs.scipy.org/doc/numpy/reference/generated/numpy.gradient.html
           
    Parameters
    ----------
        SegmentX: a numpy array
            The x-cordinates of the data in Segment Y,  does not have to be evenly spaced
        SegmentY: a numpy array
            The data to find peaks in
        peakHeightThreshold: number, optional
            Only peaks with whre the y-values minus the averge of two adjacent valleys is greater than the threshold will be returned
    
    Returns
    -------
        peakX: a numpy array
            the interpolated x-coordinate of the peaks in y
        peakY: a numpy array
            the y-coordinate of the peaks
        peakHeight: a numpy array
            the heights of the peaks relative to the average y-coordinate of the two adjacent valleys 
    
    Examples
    --------
    >>> SegmentY=np.array([1,0,3,2,1,4,1,0,1])
    >>> SegmentX=np.array([1,2,3,4,5,6,7,8,9])
    >>> PeakFind(SegmentX,SegmentY)
    (array([3.5, 6. ]), array([3, 4]), array([2.5, 3.5]))
    >>> peakX,peakY,peakHeight=PeakFind(SegmentX,SegmentY)

    >>> plt.plot(SegmentX,SegmentY,'-ok')
    >>> plt.plot(peakX,peakY,'^r')
    """
    
    deriv1=np.gradient(SegmentY, SegmentX, edge_order=1)
    #plt.plot(SegmentX,deriv1)
    #deriv1sm=scipy.signal.savgol_filter(SegmentY, 11, 2,deriv=1)
    #valleyX,vindex = Cross(SegmentX,deriv1,direction='rising')
    peakX,pindex = Cross(SegmentX,deriv1,direction='falling')
    #valleyY=np.min([SegmentY[vindex],SegmentY[vindex+1]],axis=0)
    peakY=np.max([SegmentY[pindex],SegmentY[pindex+1]],axis=0)
    #midPeaksBool=((peakX>=valleyX[0]) & (peakX<=valleyX[-1]))
    #peakX=peakX[midPeaksBool]
    #peakY=peakY[midPeaksBool]
    #valleyMean=valleyY[0:-1]+(np.diff(valleyY)/2)
    #peakHeightBool=(peakY-valleyMean)>=peakHeightThreshold
    #peakX=peakX[peakHeightBool]
    #peakY=peakY[peakHeightBool]
    #valleyMean=valleyMean[peakHeightBool]
#    return peakX,peakY,peakY-valleyMean
    if peakY.size==0:
        peakY=[np.max(SegmentY)]
        peakX=[SegmentX[np.argmax(SegmentY)]]
    return peakX,peakY

def absorbanceToTristim(waves,absorbance,Yr,gammaFlag=False):
    pixel=np.zeros((1,1,3),dtype=np.float32)
    pixel[0,0,0]=np.trapz(CIEX*illum*10**-absorbance, waves)/Yr
    pixel[0,0,1]=np.trapz(CIEY*illum*10**-absorbance, waves)/Yr
    pixel[0,0,2]=np.trapz(CIEZ*illum*10**-absorbance, waves)/Yr
    XYZ=pixel[0,0,:]
    RGB=cv2.cvtColor(pixel, cv2.COLOR_XYZ2RGB)
    RGBg=np.zeros((RGB.shape),dtype=np.float32)
    for cc in range(RGB.shape[2]):
        if RGB[0,0,cc]<=0.0031308:
            RGBg[0,0,cc]=12.92*RGB[0,0,cc]
        else:
            RGBg[0,0,cc]=1.055*RGB[0,0,cc]**(1/2.4)-0.055
        if RGBg[0,0,cc]>1:
            RGBg[0,0,cc]=1
        elif RGBg[0,0,cc]<0:
            RGBg[0,0,cc]=0
    if gammaFlag:
        HSV=cv2.cvtColor(RGBg, cv2.COLOR_RGB2HSV)[0,0,:]
        LAB=cv2.cvtColor(RGBg, cv2.COLOR_RGB2LAB)[0,0,:]
        RGB=RGB[0,0,:]
        RGBg=RGBg[0,0,:]
        rgb=np.zeros((RGBg.shape))
        rgb[0]=RGBg[0]/np.sum(RGBg)
        rgb[1]=RGBg[1]/np.sum(RGBg)
        rgb[2]=RGBg[2]/np.sum(RGBg)
        rat=np.zeros((RGBg.shape))
        rat[0]=RGBg[0]/RGBg[1]
        rat[1]=RGBg[0]/RGBg[2]
        rat[2]=RGBg[1]/RGBg[2]
    else:
        HSV=cv2.cvtColor(RGB, cv2.COLOR_RGB2HSV)[0,0,:]
        LAB=cv2.cvtColor(RGB, cv2.COLOR_RGB2LAB)[0,0,:]
        RGB=RGB[0,0,:]
        RGBg=RGBg[0,0,:]
        rgb=np.zeros((RGB.shape))
        rgb[0]=RGB[0]/np.sum(RGB)
        rgb[1]=RGB[1]/np.sum(RGB)
        rgb[2]=RGB[2]/np.sum(RGB)
        rat=np.zeros((RGB.shape))
        rat[0]=RGB[0]/RGB[1]
        rat[1]=RGB[0]/RGB[2]
        rat[2]=RGB[1]/RGB[2]
    HSV[0] = ShiftHOriginToValue(HSV[0],360,360.0/3,direction='ccw')
    HSV[0]=HSV[0]/360.0
    LAB[0]=LAB[0]/100.0
    LAB[1]=(LAB[1]+128)/255.0
    LAB[2]=(LAB[2]+128)/255.0
    RGBg=np.rint(RGBg*255)/255.0
    return RGB, HSV, LAB, XYZ, rgb, rat, RGBg

def interpolateResponse(OOWaves, CIEWaves, X):
    numCIEwaves=CIEWaves.shape[0]
#    Xind=np.zeros((numCIEwaves+int(np.floor(CIEWaves[0]))))
    Xind=np.zeros((int(np.ceil(CIEWaves[numCIEwaves-1])+1)))
    for index in range(numCIEwaves):
        newindex=int(CIEWaves[index])
        Xind[newindex]=X[index]
    numOOwaves=OOWaves.shape[0]
    Xout=np.zeros((numOOwaves))
    for i in range(numOOwaves):
        roundedlow=int(np.floor(OOWaves[i]))
        roundedhigh=roundedlow+1
        if roundedhigh<CIEWaves[numCIEwaves-1]:
            Xout[i]=((Xind[roundedhigh]-Xind[roundedlow])*(OOWaves[i]-roundedlow))+Xind[roundedlow]
    return (Xout)

def PolyReg(X,Y,order):
    coef,cov=np.polyfit(X,Y,order,cov=True)
    N=float(len(X))
    df=N-len(coef)
    stdErrors=np.sqrt(np.diagonal(cov)*(df-2)/df)
    p=np.poly1d(coef)
    yfit=p(X)
    res=np.array(Y)-yfit
    sy=np.sqrt( np.sum(res**2) / df )
    return {'coef':coef,'errors':stdErrors,'n':N,'sy':sy,'res':res,'poly':p}

def FormatSciUsingError(x,e,WithError=False,ExtraDigit=0):
    if abs(x)>=e:
        NonZeroErrorX=np.floor(np.log10(abs(e)))
        NonZeroX=np.floor(np.log10(abs(x)))
        formatCodeX="{0:."+str(int(NonZeroX-NonZeroErrorX+ExtraDigit))+"E}"
        formatCodeE="{0:."+str(ExtraDigit)+"E}"
    else:
        formatCodeX="{0:."+str(ExtraDigit)+"E}"
        formatCodeE="{0:."+str(ExtraDigit)+"E}"
    if WithError==True:
        return formatCodeX.format(x)+" (+/- "+formatCodeE.format(e)+")"
    else:
        return formatCodeX.format(x) 

def AnnotateFit(fit,axisHandle,annotationText='Eq',color='black',Arrow=False,xArrow=0,yArrow=0,xText=0.1,yText=0.9):
    c=fit['coef']
    e=fit['errors']
    t=len(c)
    if annotationText=='Eq':
        annotationText="y = "
        for order in range(t):
            exponent=t-order-1
            if exponent>=2:
                annotationText=annotationText+FormatSciUsingError(c[order],e[order])+"x$^{}$".format(exponent)+" + "
            elif exponent==1:
                annotationText=annotationText+FormatSciUsingError(c[order],e[order])+"x + "
            else:
                annotationText=annotationText+FormatSciUsingError(c[order],e[order])
        annotationText=annotationText+"\nsy={0:.1E}".format(fit['sy'])+", sens={0:.1E}".format(np.abs(c[0]/fit['sy']))
    if (Arrow==True):
        if (xArrow==0):
            xSpan=axisHandle.get_xlim()
            xArrow=np.mean(xSpan)
        if (yArrow==0):    
            yArrow=fit['poly'](xArrow)
        annotationObject=axisHandle.annotate(annotationText, 
                xy=(xArrow, yArrow), xycoords='data',
                xytext=(xText, yText),  textcoords='axes fraction',
                arrowprops={'color': color, 'width':1, 'headwidth':5},
                bbox={'boxstyle':'round', 'edgecolor':color,'facecolor':'0.8'}
                )
    else:
        xSpan=axisHandle.get_xlim()
        xArrow=np.mean(xSpan)
        ySpan=axisHandle.get_ylim()
        yArrow=np.mean(ySpan)
        annotationObject=axisHandle.annotate(annotationText, 
                xy=(xArrow, yArrow), xycoords='data',
                xytext=(xText, yText),  textcoords='axes fraction',
                ha="left", va="center",
                bbox={'boxstyle':'round', 'edgecolor':color,'facecolor':'0.8'}
                )
    annotationObject.draggable()
    
def onmouse(event,x,y,flags,params):
    global ColorRectangle,ColorRect,ix,iy,ixg,iyg,ColorRectOver,GreyRectangle,GreyRect,GreyRectOver,rebalanceToggle,RectList
    if event == cv2.EVENT_LBUTTONDOWN:
            ColorRectangle = True
            ColorRectOver = False
            ix,iy = x,y
    elif event == cv2.EVENT_RBUTTONDOWN:
        if GreyRectOver==True:
            GreyRectOver=False
        else:
            GreyRectangle = True
            GreyRectOver = False
            ixg,iyg = x,y
    elif event == cv2.EVENT_MOUSEMOVE:
        if ColorRectangle == True:
            cv2.rectangle(img,(ix,iy),(x,y),(255,255,255),2)
            ColorRect = (min(ix,x),min(iy,y),abs(ix-x),abs(iy-y))
            cv2.imshow('Result',img)
            #cv2.waitKey(1)
        if GreyRectangle == True:
            cv2.rectangle(img,(ixg,iyg),(x,y),(0,0,0),2)
            GreyRect = (min(ixg,x),min(iyg,y),abs(ixg-x),abs(iyg-y))
            cv2.imshow('Result',img)
            #cv2.waitKey(1)
    elif event == cv2.EVENT_LBUTTONUP:
        ColorRectangle = False
        if ix!=x & iy!=y:
            ColorRectOver = True
            cv2.rectangle(img,(ix,iy),(x,y),(255,255,255),2)
            ColorRect = (min(ix,x),min(iy,y),abs(ix-x),abs(iy-y))
            x1,y1,w,h = ColorRect
            RectList.append(ColorRect)
            cv2.imshow('Result',img)
        else:
            ColorRectOver = False
    elif event == cv2.EVENT_RBUTTONUP:
        GreyRectangle = False
        if ixg!=x & iyg!=y:
            GreyRectOver = True
            cv2.rectangle(img,(ixg,iyg),(x,y),(0,0,0),2)
            GreyRect = (min(ixg,x),min(iyg,y),abs(ixg-x),abs(iyg-y))
            #x1,y1,w,h = GreyRect        
            #cv2.imshow('Result',img)
            rebalanceToggle=True
        else:
            GreyRectOver = False
            
def OpenCVDisplayedHistogram(image,channel,mask,NumBins,DataMin,DataMax,x,y,w,h,DisplayImage,color,integrationWindow,labelFlag,labelText=""):
        avgVal=cv2.meanStdDev(image,mask=mask)
        histdata = cv2.calcHist([image],[channel],mask,[NumBins],[DataMin,DataMax])
        domValue=np.argmax(histdata)
        #domCount=np.max(histdata)/np.sum(histdata) 
        #sortArg=np.argsort(histdata,axis=0)
        #domValue=np.sum(histdata[sortArg[-5:][:,0]][:,0]*sortArg[-5:][:,0])/np.sum(histdata[sortArg[-5:][:,0]][:,0])
        #domCount=np.sum(histdata[sortArg[-5:][:,0]][:,0])/np.sum(histdata)
        numpixels=np.sum(np.array(mask))/np.max(np.array(mask))
        cv2.normalize(histdata, histdata, 0, h, cv2.NORM_MINMAX);
        if w>NumBins:
            binWidth = w/NumBins
        else:
            binWidth=1
        #img = np.zeros((h, NumBins*binWidth, 3), np.uint8)
        for i in range(NumBins):
            freq = int(histdata[i])
            cv2.rectangle(DisplayImage, ((i*binWidth)+x, y+h), (((i+1)*binWidth)+x, y+h-freq), color)
        if labelFlag:
            cv2.putText(img," "+labelText+" m="+'{0:.3f}'.format(domValue/float(NumBins-1)*(DataMax-DataMin))+" n="+'{0:.2E}'.format(numpixels)+" a="+'{0:.3f}'.format(avgVal[0][channel][0])+" s="+'{0:.4f}'.format(avgVal[1][channel][0]),(x,y+h+12), font, 0.4,color,1,cv2.LINE_AA)
        return (avgVal[0][channel][0],avgVal[1][channel][0],domValue/float(NumBins-1)*(DataMax-DataMin))
    
    
#    def OpenCVDisplayedHistogram(image,channel,mask,NumBins,DataMin,DataMax,x,y,w,h,DisplayImage,color,integrationWindow,labelFlag):
#        avgVal=cv2.mean(image,mask=mask)[channel]
#        pix=np.sum(mask)/255
#        histdata = cv2.calcHist([image],[channel],mask,[NumBins],[DataMin,DataMax])
#        domValue=np.argmax(histdata)
#        numpixels=sum(np.array(histdata[domValue-integrationWindow:domValue+integrationWindow+1]))
#        cv2.normalize(histdata, histdata, 0, h, cv2.NORM_MINMAX);
#        binWidth = w/NumBins
#        for i in xrange(NumBins):
#            freq = int(histdata[i])
#            cv2.rectangle(DisplayImage, ((i*binWidth)+x, y+h), (((i+1)*binWidth)+x, y+h-freq), color)
#        if labelFlag:
#            cv2.putText(img,"avg="+'{0:.3f}'.format(avgVal)+" pix="+str(pix).zfill(4),(x,y+h+20), font, 0.6,color,1,cv2.LINE_AA)
#        return (domValue/float(NumBins),numpixels)

def boxcarSmooth(data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(data, window, 'same')
        
def OpenCVDisplayedScatter(xdata,ydata,x,y,w,h,size,color,ydataRangemin=None, ydataRangemax=None,xdataRangemin=None, xdataRangemax=None):      
        if xdataRangemin==None: 
             xdataRangemin=min(xdata)       
        if xdataRangemax==None: 
             xdataRangemax=max(xdata) 
        if ydataRangemin==None: 
             ydataRangemin=min(ydata) 
        if ydataRangemax==None: 
             ydataRangemax=max(ydata)
        xdataRange=xdataRangemax-xdataRangemin
        ydataRange=ydataRangemax-ydataRangemin
        xscale=w/xdataRange
        yscale=h/ydataRange
        xdata=xdata-xdataRangemin
        ydata=ydata-ydataRangemin
        xdata=(xdata*xscale).astype(int)
        ydata=(ydata*yscale).astype(int)
        for i in range(xdata.size):
            cv2.circle(img,( xdata[i]+x,h-ydata[i]+y ), size , color, -1)
#        for i in range(xdata.size-1):
#            cv2.line(img,( xdata[i]+x,h-ydata[i]+y ), ( xdata[i+1]+x,h-ydata[i+1]+y ), color, size )
        cv2.rectangle(img,(x,y),(x+w,y+h),color,1)
        cv2.putText(img,str(round(xdataRangemax,0)),(x+w-15,y+h+15), font, 0.5,color,1,cv2.LINE_AA)
        cv2.putText(img,str(round(xdataRangemin,0)),(x-5,y+h+15), font, 0.5,color,1,cv2.LINE_AA)
        cv2.putText(img,str(round(ydataRangemax,2)),(x-40,y+10), font, 0.5,color,1,cv2.LINE_AA)
        cv2.putText(img,str(round(ydataRangemin,2)),(x-40,y+h-5), font, 0.5,color,1,cv2.LINE_AA)
#x=array([0,1.0,2.0,3.0,4.0,5.0])
#y=array([.1,11,23,32,39,50])
#pts=np.transpose(np.vstack((x, y)))
#vx, vy, cx, cy = cv2.fitLine(pts, cv2.DIST_L2, 0, .01, .01)
#slope=(cy-vy)/(cx-vx)
#intercept=cy-(slope*cx)
#src1=np.transpose(np.vstack((ones_like(x), x)))
#src2=y
#coeffs=cv2.solve(src1,src2,flags=cv2.DECOMP_SVD)
#slp=coeffs[1][1]
#inter=coeffs[1][0]
        
def RebalanceImageCV(frame,rfactor,gfactor,bfactor):
    #could I scale up to a 16-bit image here?
    if frame.dtype=='float32':
        offset=np.zeros(frame[:,:,0].shape,dtype=np.float32)
    else:
        offset=np.zeros(frame[:,:,0].shape,dtype="uint8")
    frame[:,:,0]=cv2.scaleAdd(frame[:,:,0], bfactor, offset)
    frame[:,:,1]=cv2.scaleAdd(frame[:,:,1], gfactor, offset)
    frame[:,:,2]=cv2.scaleAdd(frame[:,:,2], rfactor, offset)
    return frame

def SubsampleData(arr, n):
    end =  n * int(len(arr)/n)
    return np.mean(arr[:end].reshape(-1, n), 1)
#    cv.Reshape(arr, newCn, newRows=0)

#def ShiftHOriginToGreen(hsvROI):
#    shifthsv=np.copy(hsvROI[:,:,0])
#    shifthsv[hsvROI[:,:,0]>=120]=shifthsv[hsvROI[:,:,0]>=120]-120
#    shifthsv[hsvROI[:,:,0]<120]=shifthsv[hsvROI[:,:,0]<120]+240
#    hsvROI[:,:,0]=shifthsv
#    return hsvROI

def ShiftHOriginToGreen(hue,maxHue):
    shifthsv=np.copy(hue)
    shifthsv[hue>=maxHue/3.0]=shifthsv[hue>=maxHue/3.0]-maxHue/3.0
    shifthsv[hue<maxHue/3.0]=shifthsv[hue<maxHue/3.0]+maxHue*2/3.0
    hue=shifthsv
    return hue

def ShiftHOriginToValue(hue,maxHue,newOrigin,direction='cw'):
    shifthsv=np.copy(hue).astype('float')
    shiftAmount=maxHue-newOrigin
    shifthsv[hue<newOrigin]=shifthsv[hue<newOrigin]+shiftAmount
    shifthsv[hue>=newOrigin]=shifthsv[hue>=newOrigin]-newOrigin
    hue=shifthsv
    if direction=='ccw':
        hue=maxHue-hue
    return hue

def FindLargestContour(binaryImage, estimated_center=[0,0],x_tolerance=0,y_tolerance=0):
    if float(float(cv2.__version__[0])+float(cv2.__version__[2])/10)>=4:
        contours,hierarchy = cv2.findContours(binaryImage,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    else:
        image,contours,hierarchy = cv2.findContours(binaryImage,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    MaxArea=0
    ContourIndex=0
    LargestContour=0
    if len(contours)>=1:
        for contour in contours:
            area=cv2.contourArea(contour)
            M = cv2.moments(contour)
            cx=0
            cy=0
            if M['m00']>0:
                cx = int(M['m10']/M['m00'])
                cy = int(M['m01']/M['m00'])           
            if (area>MaxArea) & ((estimated_center==[0,0]) | ((cx>=estimated_center[0]-x_tolerance) & (cx<= estimated_center[0]+x_tolerance) & (cy>= estimated_center[1]-y_tolerance) & (cy<= estimated_center[1]+y_tolerance))   ):
                MaxArea=area
                LargestContour=ContourIndex
                centx=cx
                centy=cy
            ContourIndex=ContourIndex+1
    if MaxArea!=0:
        return contours[LargestContour],centx,centy,MaxArea
    else:
        return False,0,0,0
#    cv2.drawContours(frame,[contours[LargestContour]],0,(0,0,255),2)
#        M = cv2.moments(contours[LargestContour])
#        if M['m00']>0:
#            cx = int(M['m10']/M['m00'])
#            cy = int(M['m01']/M['m00'])
#            cv2.circle(frame,(cx,cy), 2, (0,0,255), -1)
#            ptsFound[0,0]=cx
#            ptsFound[0,1]=cy
#            waypoints=waypoints+1
#            ptsArea[0]=MaxArea

def scaleAndResize(chanToScale,lowScale,highScale):
    chanScaled=(chanToScale-lowScale)*(255/(highScale-lowScale))
    chanScaled=chanScaled.astype('uint8')
    chanScaledFrame = cv2.bitwise_and(chanScaled,chanScaled, mask= maskROI)   
    chanScaledResizeFrame = cv2.resize(chanScaledFrame, (int(chanScaledFrame.shape[1]*resFrameScale/2),int(chanScaledFrame.shape[0]*resFrameScale/2)), interpolation = cv2.INTER_AREA)
    return chanScaledResizeFrame
            
cv2.namedWindow('Result',cv2.WINDOW_GUI_NORMAL)
cv2.setMouseCallback('Result',onmouse)

linLUTfloat=np.zeros((256),dtype='float32')
linLUTint=np.zeros((256),dtype='uint8')
linLUTabs=np.zeros((256),dtype='float32')
for chan in range(256):
    val=chan/255.0
    if (val<=0.04045):
        val=val/12.92
    else:
        val=((val+0.055)/1.055)**2.4
    linLUTfloat[chan]=val
    linLUTint[chan]=int(round(val*255))
    if val==0:
        linLUTabs[chan]=255/64.0
    else:
        linLUTabs[chan]=-np.log10(val)
        

#getting rid of tkinter blank window
root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)

print('Please select a file')
    
#allows user to pick file
file_path = askopenfilename() 

#variable for start path
startPath = os.path.dirname(file_path)

#stores files to be processed
filesToProcessPhoto = []
filesToProcessSpectrum = []
transmissionPhotoFilesToProcess = []
calibrationToProcess=[]
#prompting user for what they want
mode = input('Root is ' + os.path.dirname(file_path) + '\n\n' + 'Process one file (1), all files in this directory (2), or all files in this and all subdirectories (3)' + '\n')
fileTypes = input('Process photos only (1), spectra only (2), or both photos and spectra (3)' + '\n')
processPhotos=False
processSpectra=False
if (fileTypes == '1') | (fileTypes == 1)| (fileTypes == '3') | (fileTypes == 3):
    processPhotos=True
if (fileTypes == '2') | (fileTypes == 2)| (fileTypes == '3') | (fileTypes == 3):
    processSpectra=True
    
if (mode == '1') | (mode == 1):
    if (file_path[-4:]=='.jpg') & (file_path[:3]!='pro') & (processPhotos):
            filesToProcessPhoto.append(file_path)
    if (file_path[-3:]== "csv") & (processSpectra): #if last 3 characters in filename are csv, then do the next line
            filesToProcessSpectrum.append(file_path) #joins the directory to the filename to create a filepath, then adds that filepath to the "files to process" list
    
if (mode == '2') | (mode == 2):
    for filename in os.listdir(startPath):
        if (filename[-4:]=='.jpg') & (filename[:3]!='pro') & (processPhotos):
                filesToProcessPhoto.append(os.path.join(startPath, filename))
        if (filename[-3:]== "csv") & (processSpectra): #if last 3 characters in filename are csv, then do the next line
                filesToProcessSpectrum.append(os.path.join(startPath,filename)) #joins the directory to the filename to create a filepath, then adds that filepath to the "files to process" list

if (mode == '3') | (mode == 3):
    for (dirpath, dirnames, filenames) in os.walk(startPath):
        for filename in filenames:
            #print(filename)
            if (filename[-4:]=='.jpg') & (filename[:3]!='pro') & (processPhotos):
                dir2=os.path.dirname(os.path.normpath(dirpath))
                folder=dir2[dir2.rfind('\\')+1:]
                if filename.find('Checker')==-1:
                    filesToProcessPhoto.append(os.path.join(os.path.normpath(dirpath), filename))
                else:
                    calibrationToProcess.append(os.path.join(os.path.normpath(dirpath), filename))                
            if (filename[-3:]=='csv') & (processSpectra):
                dir2=os.path.normpath(dirpath)
                folder=dir2[dir2.rfind('\\')+1:]
#                if (folder=='Transmission Spectra') | (folder=='Spectra')
                filesToProcessSpectrum.append(os.path.join(os.path.normpath(dirpath), filename))     

spectraData=[]
colorDataMean=[]
    
if processSpectra:                    
    dfHead=pd.read_csv(filesToProcessSpectrum[0],skiprows=1,skipfooter=513,delimiter=',',index_col=0,header=None,engine='python')
    dfType=pd.read_csv(filesToProcessSpectrum[0],skiprows=64,skipfooter=523,delimiter=',',header=None,engine='python')
    if dfType.iloc[0,0]=='sel_pixel_start':
        skiprows=89
    else:
        skiprows=80
    dfSpectrum=pd.read_csv(filesToProcessSpectrum[0],skiprows=skiprows,delimiter=',',index_col=0,header=None,names=['Waves','Wavenumber','Raman_Shift','Dark','Reference','Raw','Dark_Subtracted','Transmission','Absorbance','Irradiance','blk'],engine='python')
    dfSpectrum['Wavelength']=float(dfHead.loc['coefs_a0'].values[0])+float(dfHead.loc['coefs_a1'].values[0])*dfSpectrum.index+float(dfHead.loc['coefs_a2'].values[0])*dfSpectrum.index**2+float(dfHead.loc['coefs_a3'].values[0])*dfSpectrum.index**3
    waves=np.array(dfSpectrum['Wavelength'])
    rangeBoolCIE = (dfSpectrum['Wavelength'] >= 380) & (dfSpectrum['Wavelength'] <= 780) #do not change these wavelengths
    waves=waves[rangeBoolCIE]
    
    XLFileDF = pd.read_excel('all_1nm_data.xls',skiprows=63)
    CIEX=interpolateResponse(waves, XLFileDF.values[:,0],XLFileDF.values[:,5])
    CIEY=interpolateResponse(waves, XLFileDF.values[:,0],XLFileDF.values[:,6])
    CIEZ=interpolateResponse(waves, XLFileDF.values[:,0],XLFileDF.values[:,7])
    D65=interpolateResponse(waves, XLFileDF.values[:,0],XLFileDF.values[:,2])
    A=interpolateResponse(waves, XLFileDF.values[:,0],XLFileDF.values[:,1])
    XLFileDF_illum = pd.read_excel('fluorescent_illuminants.xls',skiprows=1)
    fSplineF1 = interp1d(XLFileDF_illum.values[:,0],XLFileDF_illum.values[:,1], kind='cubic')
    F1=fSplineF1(waves)
    #F1=interpolateResponse(waves, XLFileDF_illum.values[:,0],XLFileDF_illum.values[:,1])
    illum=D65
    Yr=np.trapz(CIEY*illum, waves)

    largestAbsorbance=-1           
    spectraData=[]
    rawSpectra=[]    
    rangeBool=(dfSpectrum['Wavelength']>=480)&(dfSpectrum['Wavelength']<=840)
    rawSpectra.append({'medium':'all','filename': 'all','type':'all','ri':'all','slide':'all','prep':'all','spectra':dfSpectrum['Wavelength'][rangeBool].values})
    for filename in filesToProcessSpectrum:
        baseName=os.path.basename(os.path.basename(filename))
        #reads data from array containing list of files
        info = os.stat(filename)
        cTime=datetime.datetime.fromtimestamp(info.st_ctime)
        dfHead=pd.read_csv(filename,skiprows=1,skipfooter=513,delimiter=',',index_col=0,header=None,engine='python')
        dfType=pd.read_csv(filename,skiprows=64,skipfooter=523,delimiter=',',header=None,engine='python')
        if dfType.iloc[0,0]=='sel_pixel_start':
            skiprows=89
        else:
            skiprows=80
        dfSpectrum=pd.read_csv(filename,skiprows=skiprows,delimiter=',',index_col=0,header=None,names=['Waves','Wavenumber','Raman_Shift','Dark','Reference','Raw','Dark_Subtracted','Transmission','Absorbance','Irradiance','blk'],engine='python')
        dfSpectrum['Wavelength']=float(dfHead.loc['coefs_a0'].values[0])+float(dfHead.loc['coefs_a1'].values[0])*dfSpectrum.index+float(dfHead.loc['coefs_a2'].values[0])*dfSpectrum.index**2+float(dfHead.loc['coefs_a3'].values[0])*dfSpectrum.index**3
    
        dfSpectrum['Absorbance']=dfSpectrum.iloc[:, 8]
        
        #Cantrell
        RGB, HSV, LAB, XYZ, rgb, rat, RGBg=absorbanceToTristim(dfSpectrum['Wavelength'][rangeBoolCIE],dfSpectrum['Absorbance'][rangeBoolCIE],Yr,gammaFlag=True)
        
        #truncates data to region of interest (the plasmonic peak region) based on wavelength and then plots it
        rangeBool=(dfSpectrum['Wavelength']>=480)&(dfSpectrum['Wavelength']<=840)#This takes the x data from dfSpectrum and makes a Boolean array with values "true" between 500 and 600 and "false" everywhere else
        #ax.plot(dfSpectrum['Wavelength'][rangeBool],dfSpectrum['Absorbance'][rangeBool],color=RGBg)
        
       
        #returns full row of data for row where absorbance is a maximum within the wavelength range defined in rangeBool
        maxVals = dfSpectrum.loc[dfSpectrum['Absorbance'][rangeBool].idxmax()]
        if maxVals['Absorbance']>largestAbsorbance:
            largestAbsorbance=maxVals['Absorbance']
        
        #truncates absorbance array to only the upper 25% of absorbance values
        rangeBoolIntens=(dfSpectrum['Absorbance']>=0.75*maxVals['Absorbance'])&(dfSpectrum['Wavelength']>=480)&(dfSpectrum['Wavelength']<=840)
        #y = dfSpectrum['Absorbance'][rangeBoolIntens]
        #x = np.array(dfSpectrum['Wavelength'][rangeBoolIntens])
        
        #calculates Gaussian Fit to upper 25% of spectrum
        #popt does parameter optimization for the gaussian function and specifies the fit will be for the x and y data from rangeBoolIntens
        #with initial parameter guesses for mu, sigma, amp of 520, 30, and 87
        popt, pcov = curve_fit(gaussian,dfSpectrum['Wavelength'][rangeBoolIntens],dfSpectrum['Absorbance'][rangeBoolIntens],p0=[dfSpectrum['Wavelength'][rangeBoolIntens].iloc[np.argmax(dfSpectrum['Absorbance'][rangeBoolIntens])],30,87])
        #this creates a gaussian curve (y data) using the optimized parameters determined using popt above
        GaussData=gaussian(dfSpectrum['Wavelength'],popt[0],popt[1],popt[2])
        #axGauss.plot(dfSpectrum['Wavelength'][rangeBool],GaussData[rangeBool],color=RGBg)
        #axGauss.plot(dfSpectrum['Wavelength'],dfSpectrum['Absorbance'],color=RGBg*0.5)
        
        #smoothing data, 2nd argument is the pixels to smooth, 3rd argument is polynomial order)
        smoothedData = scipy.signal.savgol_filter(np.array(dfSpectrum['Absorbance']), 31, 5)
        
        #Cantrell
        peakdata=PeakFind(np.array(dfSpectrum['Wavelength'][rangeBool]),np.array(dfSpectrum['Absorbance'][rangeBool]),peakHeightThreshold=0)
            
        derMaxAbsIndex = np.argmax(peakdata[1])
        derDataLambda = peakdata[0][derMaxAbsIndex]
        derDataAbsorbance = peakdata[1][derMaxAbsIndex]
        
        deriv1sm=scipy.signal.savgol_filter(np.array(dfSpectrum['Absorbance']), 31, 2,deriv=1)
        deriv2sm=scipy.signal.savgol_filter(np.array(dfSpectrum['Absorbance']), 31, 2,deriv=2)
        #axDer1.plot(dfSpectrum['Wavelength'][rangeBool],deriv1sm[rangeBool],color=RGBg)
        #axDer1.set_title('1st Derivative')
        #axDer2.plot(dfSpectrum['Wavelength'][rangeBool],deriv2sm[rangeBool],color=RGBg)
        #axDer2.set_title('2nd Derivative')
    
        gaussLambda=popt[0]
        gaussSigma=popt[1]
        targetIndex = np.argmin(np.abs(np.array(dfSpectrum['Wavelength'])-gaussLambda))
        targetX = dfSpectrum['Wavelength'][targetIndex]
        targetY =deriv2sm[targetIndex]
    
    
        #calculating smoothed max
        y = smoothedData[rangeBool] #get smoothed absorbance spectrum values
        x = np.array(dfSpectrum['Wavelength'][rangeBool].values) #gets wavelengths for absorbance values
        smoothedMaxAbsIndex = np.argmax(y) #finds index where absorbance is a maximum
        smoothedDataLambda = x[smoothedMaxAbsIndex] #gives wavelength at that index
        smoothedDataAbsorbance = y[smoothedMaxAbsIndex] #gives absorbance at that index
    
        ri = 0
        mainfolder=os.path.split(os.path.split(os.path.split(filename)[0])[0])[1]
        gen=mainfolder[0:3]
        #first 3 characters of folder name that contains photos needs to be the RI medium
        sensorMedium=baseName[0:3]
        slide=mainfolder[6]
        #last 6 characters of folder name containing ALL data needs to be the date
        prep=mainfolder[-6:]
        sam=mainfolder[8:14]
        surface=mainfolder[15:18]
    
        
        riValues=riValues_old
        if prep=='220706':
            riValues=riValues_220706
        if prep=='220705':
            riValues=riValues_220705
        if prep=='220701' and slide=='A' or 'B' or 'C':
            riValues=riValues_220701A
        if prep=='220701' and slide=='D':
            riValues=riValues_220701D
        if prep=='220628':
            riValues=riValues_220628
        if prep=='062719':
            if gen=='G7Au':
                if slide=='C':
                    riValues=riValues_G7Au_C_062719
                if slide=='D':
                    riValues=riValues_G7Au_D_062719
            if gen=='G9Au':
                if slide=='C':
                    riValues=riValues_G9Au_C_062719
        if prep=='062819':        
            riValues=riValues_062819
        if prep=='070119':        
            riValues=riValues_070119
        if prep=='070319':        
            riValues=riValues_070319
        if prep=='070519':        
            riValues=riValues_070519        
        if prep == '071019':
            riValues = riValues_071019
        if prep == '071219':
            riValues = riValues_071219
        if prep == '071519':
            riValues = riValues_071519
        if prep == '071619':
            riValues = riValues_071619
        if prep == '100819':
            riValues = riValues_100819
        if prep == '112319':
            riValues = riValues_112319
            
        if sensorMedium  == 'air':
            ri = riValues[0]
        if (sensorMedium  == 'h20') | (sensorMedium  == 'h2o'):
            ri = riValues[1]
        if sensorMedium  == '10%':
            ri = riValues[2]
        if sensorMedium  == '20%':
            ri = riValues[3]
        if sensorMedium  == '30%':
            ri = riValues[4]
        if sensorMedium  == '40%':
            ri = riValues[5]
        if sensorMedium  == '50%':
            ri = riValues[6]
        spectraData.append({'medium':sensorMedium,'filename': filename,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'R':RGB[0],'G':RGB[1],'B':RGB[2],'H':HSV[0],'S':HSV[1],'V':HSV[2],'L*':LAB[0],'a*':LAB[1],'b*':LAB[2],'r':rgb[0],'g':rgb[1],'b':rgb[2],'gaussLambda':gaussLambda,'gaussSigma':gaussSigma,'Max Abs':smoothedDataAbsorbance,'smoothedDataLambda':smoothedDataLambda,'derDataLambda':derDataLambda,'type':gen,'ri':ri,'slide':slide,'prep':prep,'SAM':sam,'Surface':surface})
        rawSpectra.append({'medium':sensorMedium,'filename': filename,'type':gen,'SAM':sam,'Surface':surface,'ri':ri,'slide':slide,'prep':prep,'spectra':dfSpectrum['Absorbance'][rangeBool].values})
if processPhotos:
    colorDataMean=[]
    colorDataMost=[]
    colorDataMeanOriginal=[]
    colorDataMostOriginal=[]
    ParameterStats=np.zeros((16,9,len(filesToProcessPhoto)+1))
    counter=0
    for fileName in filesToProcessPhoto:
            recNum=0
            counter=counter+1
            print('Processing file '+str(counter)+' of '+str(len(filesToProcessPhoto)))
            info = os.stat(fileName)
            cTime=datetime.datetime.fromtimestamp(info.st_ctime)
            frame=cv2.imread(fileName,cv2.IMREAD_ANYDEPTH|cv2.IMREAD_COLOR)
            originalFrame=np.copy(frame)
    
            #frameFloat=np.zeros((frame.shape),dtype=np.float32)
            frameFloat=np.divide(frame,255,dtype=np.float32)
            frame=frameFloat
            imgScale=50.0/max(frame.shape[1]/16.0,frame.shape[0]/9.0)
            img = np.zeros((900, 1700, 3), np.uint8)
            CamFrame = cv2.resize(frame, (int(frame.shape[1]*imgScale),int(frame.shape[0]*imgScale)), interpolation = cv2.INTER_AREA)
            img[0:int(frame.shape[0]*imgScale),0:int(frame.shape[1]*imgScale),:]=CamFrame*255
            ColorRect=[int(frame.shape[1]*imgScale*0.1),int(frame.shape[0]*imgScale*0.1/2),int(frame.shape[1]*imgScale*0.8),int(frame.shape[0]*imgScale*1.00/2)]
    #            ColorRect=[int(frame.shape[1]*imgScale*0.1),int(frame.shape[0]*imgScale*0.1/2),int(frame.shape[1]*imgScale*0.8),int(frame.shape[0]*imgScale/2)]
            RectList=[ColorRect]
            x1,y1,w,h = ColorRect   
            #cv2.rectangle(img,(x1,y1),(x1+w,y1+h),(255,255,255),2)
            wbMask=np.zeros((frame.shape[0],frame.shape[1],1), np.uint8)
    #        wbMask[int(frame.shape[0]/3):int(frame.shape[0]*2/3),:]=255
            wbMask[int(frame.shape[0]/2):int(frame.shape[0]*2/3),int(frame.shape[1]*0.1):int(frame.shape[1]*0.9)]=255
            #wbMask[0:int(frame.shape[0]/2),0:int(frame.shape[1]*.1)]=255
            #wbMask[int(frame.shape[0]*.9/2):int(frame.shape[0]/2),:]=255
            #wbMask[0:int(frame.shape[0]/2),int(frame.shape[1]*.9):int(frame.shape[1])]=255
            wbMaskScale = cv2.resize(wbMask, (int(wbMask.shape[1]*imgScale),int(wbMask.shape[0]*imgScale)), interpolation = cv2.INTER_AREA)
            RGBGreyROI=cv2.mean(frame,wbMask)
            bscale=RGBGreyROI[0]
            gscale=RGBGreyROI[1]
            rscale=RGBGreyROI[2]
            #scalemax=max(rscale,gscale,bscale)
            if min(rscale,gscale,bscale)!=0:
                if scalemin!=0:
                    rfactor=float(scalemin)/float(rscale)
                    gfactor=float(scalemin)/float(gscale)
                    bfactor=float(scalemin)/float(bscale)
                else:
                    rfactor=float(min(rscale,gscale,bscale))/float(rscale)
                    gfactor=float(min(rscale,gscale,bscale))/float(gscale)
                    bfactor=float(min(rscale,gscale,bscale))/float(bscale)
            else:
                rfactor=float(1)
                gfactor=float(1)
                bfactor=float(1)
            if rebalanceToggle:
                frame=RebalanceImageCV(frame,rfactor,gfactor,bfactor)
                CamFrame = cv2.resize(frame, (int(frame.shape[1]*imgScale),int(frame.shape[0]*imgScale)), interpolation = cv2.INTER_AREA)
                img[0:int(frame.shape[0]*imgScale),int(frame.shape[1]*imgScale):int(frame.shape[1]*imgScale*2),:]=CamFrame*255
    
            if ColorRectOver:
                x1,y1,w,h = ColorRect
                #cv2.rectangle(img,(int(frame.shape[1]*imgScale)+(recNum*260),0),(int(frame.shape[1]*imgScale)+(recNum*260)+260,900),(255-(recNum*75),255-(recNum*75),255-(recNum*75)),2)
                #cv2.rectangle(img,(x1,y1),(x1+w,y1+h),(255,255,255),2)
                #cv2.rectangle(img,(x1,y1+int(frame.shape[0]*imgScale)),(x1+w,y1+h+int(frame.shape[0]*imgScale)),(255,255,255),2)
                #img[y1+450:y1+450+h,x1:x1+w,:]= CamFrame[y1:y1+h,x1:x1+w,:]
                x1f=int(x1/imgScale)
                y1f=int(y1/imgScale)
                wf=int(w/imgScale)
                hf=int(h/imgScale)
                hsvFrame = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)
                labFrame = cv2.cvtColor(frame, cv2.COLOR_BGR2LAB)
                hsvOriginalFrame = cv2.cvtColor(originalFrame, cv2.COLOR_BGR2HSV)
                labOriginalFrame = cv2.cvtColor(originalFrame, cv2.COLOR_BGR2LAB)
                maskFrame = cv2.inRange(hsvFrame, lower_lim, upper_lim)
                maskFrameScale = cv2.resize(maskFrame, (int(frame.shape[1]*imgScale),int(frame.shape[0]*imgScale)), interpolation = cv2.INTER_AREA)
    
                largestRegion,cntX,cntY,area=FindLargestContour(maskFrame,[maskFrame.shape[1]/2,maskFrame.shape[0]/4],maskFrame.shape[1]*0.2,maskFrame.shape[0]*0.2)
                if (cntX!=0)|(cntY!=0):
                    xc,yc,wc,hc = cv2.boundingRect(largestRegion)
                    rgbROI = frame[yc:yc+hc, xc:xc+wc]
                    rgbROIoriginal = originalFrame[yc:yc+hc, xc:xc+wc]
                    (x,y),radius = cv2.minEnclosingCircle(largestRegion)
                    center = (int(x)-xc,int(y)-yc)
                    radius = int(radius*0.8)
                    centerMask=np.zeros((rgbROI.shape[0],rgbROI.shape[1]), dtype=np.uint8)
                    cv2.circle(centerMask,center,radius,(255),-1)
                    #cv2.rectangle(img,(x,y),(x+w,y+h),(0,255,0),2)
                    #cv2.drawContours(frame,[largestRegion],0,(0,0,255),2)
                else:
                    rgbROI = frame[y1f:y1f+hf, x1f:x1f+wf]
                    rgbROIoriginal = originalFrame[y1f:y1f+hf, x1f:x1f+wf]
                    centerMask=np.zeros((rgbROI.shape[0],rgbROI.shape[1]), dtype=np.uint8)
              
                balanceFrame = cv2.bitwise_and(frame,frame, mask= wbMask)
                balanceFrameScale = cv2.resize(balanceFrame, (int(balanceFrame.shape[1]*imgScale),int(balanceFrame.shape[0]*imgScale)), interpolation = cv2.INTER_AREA)
    
                hsvFrame[:,:,0]=hsvFrame[:,:,0]/360.0
                labFrame[:,:,0]=labFrame[:,:,0]/100.0
                labFrame[:,:,1]=(labFrame[:,:,1]+127)/255.0
                labFrame[:,:,2]=(labFrame[:,:,2]+127)/255.0
                
                #hsvOriginalFrame[:,:,0]=hsvOriginalFrame[:,:,0]/360.0
                #labOriginalFrame[:,:,0]=labOriginalFrame[:,:,0]/100.0
                #labOriginalFrame[:,:,1]=(labOriginalFrame[:,:,1]+127)/255.0
                #labOriginalFrame[:,:,2]=(labOriginalFrame[:,:,2]+127)/255.0
                
                wbOffset=380
                localWBoffset=260+380
                row=0                
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(frame,2,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(0,0,255),5,True,"R")
                #Rw,numRw=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,0,128),5,False)
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(frame,1,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(0,255,0),5,True,"G")
                #Gw,numGw=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,128,0),5,False)
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(frame,0,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(255,50,50),5,True,"B")
                #Bw,numBw=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(recNum*258),5+(row*55),256,40,img,(128,0,0),5,False)
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(hsvFrame,0,wbMask,360,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),180,40,img,(255,255,0),1,True,"H")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(hsvFrame,1,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(200,200,200),1,True,"S")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(hsvFrame,2,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(128,128,128),1,True,"V")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(labFrame,0,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(255,255,255),1,True,"L")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(labFrame,1,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(255,0,255),1,True,"a")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(labFrame,2,wbMask,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(0,255,255),1,True,"b")
                row+=1
    
                row=0                
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(originalFrame,2,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(0,0,128),5,False,"R")
                #Rw,numRw=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,0,128),5,False)
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(originalFrame,1,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(0,128,0),5,False,"G")
                #Gw,numGw=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,128,0),5,False)
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(originalFrame,0,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(128,25,25),5,False,"B")
                #Bw,numBw=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(recNum*258),5+(row*55),256,40,img,(128,0,0),5,False)
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(hsvOriginalFrame,0,wbMask,180,0,179,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),180,40,img,(128,128,0),1,False,"H")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(hsvOriginalFrame,1,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(100,100,100),1,False,"S")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(hsvOriginalFrame,2,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(64,64,64),1,False,"V")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(labOriginalFrame,0,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(128,128,128),1,False,"L")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(labOriginalFrame,1,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(128,0,128),1,False,"a")
                row+=1
                ParameterStats[row,6,recNum],ParameterStats[row,7,recNum],ParameterStats[row,8,recNum]=OpenCVDisplayedHistogram(labOriginalFrame,2,wbMask,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset+wbOffset),5+(row*58),256,40,img,(0,128,128),1,False,"b")
                row+=1 
            
                   
                if localRebalanceToggle:
                    hsvROI = cv2.cvtColor(rgbROI, cv2.COLOR_BGR2HSV)
                    maskGreyROI = cv2.inRange(hsvROI, lower_lim_grey, upper_lim_grey)
                    #maskGreyROI = cv2.bitwise_not(maskROI)
                    #radius = int(radius*0.8)
                    #centerMask=np.zeros((rgbROI.shape[0],rgbROI.shape[1]), dtype=np.uint8)
                    cv2.circle(maskGreyROI,center,int(radius/0.8*1.1),(0),-1)
                    
                    RGBGrey=cv2.mean(rgbROI, mask=maskGreyROI)
                    Bw=RGBGrey[0]
                    Gw=RGBGrey[1]
                    Rw=RGBGrey[2]
                    if min(Rw,Gw,Bw)!=0:
                        if scalemin!=0:
                            Rr=float(scalemin)/float(Rw)
                            Gr=float(scalemin)/float(Gw)
                            Br=float(scalemin)/float(Bw)
                        else:
                            Rr=float(min([Bw,Gw,Rw]))/float(Rw)
                            Gr=float(min([Bw,Gw,Rw]))/float(Gw)
                            Br=float(min([Bw,Gw,Rw]))/float(Bw)
                    else:
                        Rr=1
                        Gr=1
                        Br=1
                    rgbROI=RebalanceImageCV(rgbROI,Rr,Gr,Br)
                else:
                    Rr=1
                    Gr=1
                    Br=1                  
                    maskGreyROI = cv2.inRange(hsvROI, lower_lim_grey, upper_lim_grey)
                
                hsvROI = cv2.cvtColor(rgbROI, cv2.COLOR_BGR2HSV)
                hsvROIoriginal = cv2.cvtColor(rgbROIoriginal, cv2.COLOR_BGR2HSV)
                
                hsvROI[:,:,0] = ShiftHOriginToValue(hsvROI[:,:,0],360,360.0/3,direction='ccw')
                #hsvROI = ShiftHOriginToGreen(hsvROI)
                labROI = cv2.cvtColor(rgbROI, cv2.COLOR_BGR2Lab)
                labROIoriginal = cv2.cvtColor(rgbROIoriginal, cv2.COLOR_BGR2LAB)
                logsrgbROI=cv2.LUT(np.rint(rgbROI*255).astype('uint8'), linLUTabs)*64
                logsrgbROIoriginal=cv2.LUT(rgbROIoriginal, linLUTabs)*64
    
                #absROI=np.float32(-np.log10(rgbROI/255.0))
                #absROI[absROI>2]=2
                maskROIsquare = cv2.inRange(hsvROI, lower_lim, upper_lim)
                maskROI=cv2.bitwise_and(maskROIsquare,centerMask)
                resFrame = cv2.bitwise_and(rgbROI,rgbROI, mask= maskROI)

                
                largestSide=max(resFrame.shape[0],resFrame.shape[1])
                resFrameScale=(img.shape[0]/4)/largestSide

                CamMaskFrame = cv2.resize(resFrame, (int(resFrame.shape[1]*resFrameScale),int(resFrame.shape[0]*resFrameScale)), interpolation = cv2.INTER_AREA)
                img[460:460+CamMaskFrame.shape[0],0:0+CamMaskFrame.shape[1],:]= CamMaskFrame*255
                wbFrame = cv2.bitwise_and(rgbROI,rgbROI, mask= maskGreyROI)
                CamWBFrame = cv2.resize(wbFrame, (int(wbFrame.shape[1]*resFrameScale),int(wbFrame.shape[0]*resFrameScale)), interpolation = cv2.INTER_AREA)
                img[460:460+CamWBFrame.shape[0],CamMaskFrame.shape[1]:CamMaskFrame.shape[1]+CamWBFrame.shape[1],:]= CamWBFrame*255
                hsvROI[:,:,0]=hsvROI[:,:,0]/360.0
                labROI[:,:,0]=labROI[:,:,0]/100.0
                labROI[:,:,1]=(labROI[:,:,1]+127)/255.0
                labROI[:,:,2]=(labROI[:,:,2]+127)/255.0
                RGBGrey=cv2.mean(rgbROI, mask=maskGreyROI)
                HSVMean=cv2.mean(hsvROI, mask=maskROI)
                RGBMean=cv2.mean(rgbROI, mask=maskROI)
                LABMean=cv2.mean(labROI, mask=maskROI)
                
                row=0                
                #OpenCVDisplayedHistogram(image,channel,mask,NumBins,DataMin,DataMax,x,y,w,h,DisplayImage,color,integrationWindow,labelFlag,labelText="")
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(rgbROI,2,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,0,255),5,True,"R")
                #Rw,numRw=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,0,128),5,False)
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(rgbROI,1,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,255,0),5,True,"G")
                #Gw,numGw=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,128,0),5,False)
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(rgbROI,0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(255,50,50),5,True,"B")
                #Bw,numBw=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(recNum*258),5+(row*55),256,40,img,(128,0,0),5,False)
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(hsvROI,0,maskROI,360,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),180,40,img,(255,255,0),1,True,"H")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(hsvROI,1,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(200,200,200),1,True,"S")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(hsvROI,2,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(128,128,128),1,True,"V")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(labROI,0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(255,255,255),1,True,"L")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(labROI,1,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(255,0,255),1,True,"a")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(labROI,2,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,255,255),1,True,"b")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(np.float32(rgbROI[:,:,2]/np.sum(rgbROI,axis=2).astype(float)),0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(100,100,255),5,True,"r")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(np.float32(rgbROI[:,:,1]/np.sum(rgbROI,axis=2).astype(float)),0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(100,255,100),5,True,"g")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(np.float32(rgbROI[:,:,0]/np.sum(rgbROI,axis=2).astype(float)),0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(255,150,150),5,True,"b")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(logsrgbROI,2,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,255,150),5,True,"Ra")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(logsrgbROI,1,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(255,0,150),5,True,"Ga")
                row+=1
                ParameterStats[row,0,recNum],ParameterStats[row,1,recNum],ParameterStats[row,2,recNum]=OpenCVDisplayedHistogram(logsrgbROI,0,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(150,255,0),5,True,"Ba")
                row+=1
    
                row=0                
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,0,255),5,True,"R")
                #Rw,numRw=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,0,128),5,False)
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,255,0),5,True,"G")
                #Gw,numGw=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,128,0),5,False)
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(255,50,50),5,True,"B")
                #Bw,numBw=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(recNum*258),5+(row*55),256,40,img,(128,0,0),5,False)
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(hsvROI,0,maskGreyROI,360,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),180,40,img,(255,255,0),1,True,"H")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(hsvROI,1,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(200,200,200),1,True,"S")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(hsvROI,2,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(128,128,128),1,True,"V")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(labROI,0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(255,255,255),1,True,"L")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(labROI,1,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(255,0,255),1,True,"a")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(labROI,2,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,255,255),1,True,"b")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(np.float32(rgbROI[:,:,2]/np.sum(rgbROI,axis=2).astype(float)),0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(100,100,255),5,True,"r")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(np.float32(rgbROI[:,:,1]/np.sum(rgbROI,axis=2).astype(float)),0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(100,255,100),5,True,"g")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(np.float32(rgbROI[:,:,0]/np.sum(rgbROI,axis=2).astype(float)),0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(255,150,150),5,True,"b")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(logsrgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,255,150),5,True,"Ra")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(logsrgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(255,0,150),5,True,"Ga")
                row+=1
                ParameterStats[row,3,recNum],ParameterStats[row,4,recNum],ParameterStats[row,5,recNum]=OpenCVDisplayedHistogram(logsrgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(150,255,0),5,True,"Ba")
                row+=1
    
    
                row=0                
                #OpenCVDisplayedHistogram(image,channel,mask,NumBins,DataMin,DataMax,x,y,w,h,DisplayImage,color,integrationWindow,labelFlag,labelText="")
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(rgbROIoriginal,2,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,0,128),5,False,"R")
                #Rw,numRw=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,0,128),5,False)
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(rgbROIoriginal,1,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,128,0),5,False,"G")
                #Gw,numGw=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,128,0),5,False)
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(rgbROIoriginal,0,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(128,25,25),5,False,"B")
                #Bw,numBw=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(recNum*258),5+(row*55),256,40,img,(128,0,0),5,False)
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(hsvROIoriginal,0,maskROI,180,0,179,int(frame.shape[1]*imgScale)+(260),5+(row*58),180,40,img,(128,128,0),1,False,"H")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(hsvROIoriginal,1,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(100,100,100),1,False,"S")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(hsvROIoriginal,2,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(64,64,64),1,False,"V")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(labROIoriginal,0,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(128,128,128),1,False,"L")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(labROIoriginal,1,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(128,0,128),1,False,"a")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(labROIoriginal,2,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,128,128),1,False,"b")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(np.float32(rgbROIoriginal[:,:,2]/np.sum(rgbROIoriginal,axis=2).astype(float)),0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(50,50,128),5,False,"r")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(np.float32(rgbROIoriginal[:,:,1]/np.sum(rgbROIoriginal,axis=2).astype(float)),0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(50,128,50),5,False,"g")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(np.float32(rgbROIoriginal[:,:,0]/np.sum(rgbROIoriginal,axis=2).astype(float)),0,maskROI,256,0,1,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(128,75,75),5,False,"b")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(logsrgbROIoriginal,2,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(0,128,75),5,False,"Ra")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(logsrgbROIoriginal,1,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(128,0,75),5,False,"Ga")
                row+=1
                ParameterStats[row,0,recNum+1],ParameterStats[row,1,recNum+1],ParameterStats[row,2,recNum+1]=OpenCVDisplayedHistogram(logsrgbROIoriginal,0,maskROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*58),256,40,img,(75,128,0),5,False,"Ba")
                row+=1
    
                row=0                
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(rgbROIoriginal,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,0,128),5,False,"R")
                #Rw,numRw=OpenCVDisplayedHistogram(rgbROI,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,0,128),5,False)
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(rgbROIoriginal,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,128,0),5,False,"G")
                #Gw,numGw=OpenCVDisplayedHistogram(rgbROI,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(260),5+(row*55),256,40,img,(0,128,0),5,False)
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(rgbROIoriginal,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(128,25,25),5,False,"B")
                #Bw,numBw=OpenCVDisplayedHistogram(rgbROI,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(recNum+1*258),5+(row*55),256,40,img,(128,0,0),5,False)
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(hsvROIoriginal,0,maskGreyROI,180,0,179,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),180,40,img,(128,128,0),1,False,"H")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(hsvROIoriginal,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(100,100,100),1,False,"S")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(hsvROIoriginal,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(64,64,64),1,False,"V")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(labROIoriginal,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(128,128,128),1,False,"L")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(labROIoriginal,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(128,0,128),1,False,"a")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(labROIoriginal,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,128,128),1,False,"b")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(np.float32(rgbROIoriginal[:,:,2]/np.sum(rgbROIoriginal,axis=2).astype(float)),0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(50,50,128),5,False,"r")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(np.float32(rgbROIoriginal[:,:,1]/np.sum(rgbROIoriginal,axis=2).astype(float)),0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(50,128,50),5,False,"g")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(np.float32(rgbROIoriginal[:,:,0]/np.sum(rgbROIoriginal,axis=2).astype(float)),0,maskGreyROI,256,0,1,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(128,75,75),5,False,"b")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(logsrgbROIoriginal,2,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(0,128,75),5,False,"Ra")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(logsrgbROIoriginal,1,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(128,0,75),5,False,"Ga")
                row+=1
                ParameterStats[row,3,recNum+1],ParameterStats[row,4,recNum+1],ParameterStats[row,5,recNum+1]=OpenCVDisplayedHistogram(logsrgbROIoriginal,0,maskGreyROI,256,0,255,int(frame.shape[1]*imgScale)+(localWBoffset),5+(row*58),256,40,img,(75,128,0),5,False,"Ba")
                row+=1

                ParameterStats[15,:,recNum]=area

#'R':ParameterStats[0,0,recNum],'G':ParameterStats[1,0,recNum],'B':ParameterStats[2,0,recNum],'H':ParameterStats[3,0,recNum]
                chanToScale=hsvROI[:,:,0]
#                chanMask = cv2.bitwise_and(chanToScale,chanToScale, mask= maskROI)
#                if chanMask[chanMask!=0].size!=0:
#                    lowScale=np.min(chanMask[chanMask!=0])
#                else:
#                    lowScale=0
#                highScale=np.max(chanMask)
                lowScale=ParameterStats[3,0,0]-(ParameterStats[3,1,0]*3)
                highScale=ParameterStats[3,0,0]+(ParameterStats[3,1,0]*3)
                chanScaledResizeFrame=scaleAndResize(chanToScale,lowScale,highScale)
                img[470+int(img.shape[0]/4):470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0],10:10+chanScaledResizeFrame.shape[1],0]=chanScaledResizeFrame
                img[470+int(img.shape[0]/4):470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0],10:10+chanScaledResizeFrame.shape[1],1]=chanScaledResizeFrame
                img[470+int(img.shape[0]/4):470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0],10:10+chanScaledResizeFrame.shape[1],2]=chanScaledResizeFrame
                cv2.putText(img,'l={0:.2f}'.format(lowScale)+' h={0:.2f}'.format(highScale),(10,470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0]), font, 0.4,(255,255,255),1,cv2.LINE_AA)
                
                chanToScale=rgbROI[:,:,0]
#                chanMask = cv2.bitwise_and(chanToScale,chanToScale, mask= maskROI)
#                if chanMask[chanMask!=0].size!=0:
#                    lowScale=np.min(chanMask[chanMask!=0])
#                else:
#                    lowScale=0
#                highScale=np.max(chanMask)
                lowScale=ParameterStats[2,0,0]-(ParameterStats[2,1,0]*3)
                highScale=ParameterStats[2,0,0]+(ParameterStats[2,1,0]*3)
                chanScaledResizeFrame=scaleAndResize(chanToScale,lowScale,highScale)
                img[470+int(img.shape[0]/4):470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0],10+chanScaledResizeFrame.shape[1]:10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1],0]= chanScaledResizeFrame
                cv2.putText(img,'l={0:.2f}'.format(lowScale)+' h={0:.2f}'.format(highScale),(10+chanScaledResizeFrame.shape[1],470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0]), font, 0.4,(255,0,0),1,cv2.LINE_AA)

                chanToScale=rgbROI[:,:,1]
#                chanMask = cv2.bitwise_and(chanToScale,chanToScale, mask= maskROI)
#                if chanMask[chanMask!=0].size!=0:
#                    lowScale=np.min(chanMask[chanMask!=0])
#                else:
#                    lowScale=0
#                highScale=np.max(chanMask)
                lowScale=ParameterStats[1,0,0]-(ParameterStats[1,1,0]*3)
                highScale=ParameterStats[1,0,0]+(ParameterStats[1,1,0]*3)
                chanScaledResizeFrame=scaleAndResize(chanToScale,lowScale,highScale)
                img[470+int(img.shape[0]/4):470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0],10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]:10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1],1]= chanScaledResizeFrame
                cv2.putText(img,'l={0:.2f}'.format(lowScale)+' h={0:.2f}'.format(highScale),(10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1],470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0]), font, 0.4,(0,255,0),1,cv2.LINE_AA)

                chanToScale=rgbROI[:,:,2]
#                chanMask = cv2.bitwise_and(chanToScale,chanToScale, mask= maskROI)
#                if chanMask[chanMask!=0].size!=0:
#                    lowScale=np.min(chanMask[chanMask!=0])
#                else:
#                    lowScale=0
#                highScale=np.max(chanMask)
                lowScale=ParameterStats[0,0,0]-(ParameterStats[0,1,0]*3)
                highScale=ParameterStats[0,0,0]+(ParameterStats[0,1,0]*3)
                chanScaledResizeFrame=scaleAndResize(chanToScale,lowScale,highScale)
                img[470+int(img.shape[0]/4):470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0],10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]:10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1],2]= chanScaledResizeFrame
                cv2.putText(img,'l={0:.2f}'.format(lowScale)+' h={0:.2f}'.format(highScale),(10+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1]+chanScaledResizeFrame.shape[1],470+int(img.shape[0]/4)+chanScaledResizeFrame.shape[0]), font, 0.4,(0,0,255),1,cv2.LINE_AA)

                keypress=cv2.waitKey(1) & 0xFF
                if keypress == ord('q'):
                    break
    #            gen=os.path.basename(fileName)[:2]
    #            ri=os.path.basename(fileName)[-7:-4]
    #            slide=os.path.basename(fileName)[os.path.basename(fileName).lower().find('sample')+6:os.path.basename(fileName).lower().find('sample')+7]
                ri = 0
                mainfolder=os.path.split(os.path.split(os.path.split(filename)[0])[0])[1]
                gen=mainfolder[0:3]
                sensorMedium=os.path.dirname(fileName)[-3:]
                slide=mainfolder[6]
                #last 10 characters of folder name should be the date of experiment
                prep=mainfolder[-6:]
                sam=mainfolder[8:14]
                surface=mainfolder[15:18]
    
                dir2=os.path.dirname(fileName)
                folder=dir2[dir2.rfind('\\')+1:]
                
                riValues=riValues_old
                if prep=='220706':
                    riValues=riValues_220706
                if prep=='220705':
                    riValues=riValues_220705
                if prep=='220701' and slide=='A' or 'B' or 'C':
                    riValues=riValues_220701A
                if prep=='220701' and slide=='D':
                    riValues=riValues_220701D
                if prep=='220628':
                    riValues=riValues_220628
                if prep=='062719':
                    if gen=='G7Au':
                        if slide=='C':
                            riValues=riValues_G7Au_C_062719
                        if slide=='D':
                            riValues=riValues_G7Au_D_062719
                    if gen=='G9Au':
                        if slide=='C':
                            riValues=riValues_G9Au_C_062719
                if prep=='062819':        
                    riValues=riValues_062819
                if prep=='070119':        
                    riValues=riValues_070119
                if prep=='070319':        
                    riValues=riValues_070319
                if prep=='070519':        
                    riValues=riValues_070519  
                if prep == '071019':
                    riValues = riValues_071019
                if prep == '071219':
                    riValues = riValues_071219
                if prep == '071519':
                    riValues = riValues_071519
                if prep == '071619':
                    riValues = riValues_071619
                if prep == '100819':
                    riValues = riValues_100819 
                if prep == '112319':
                    riValues = riValues_112319
            
                if sensorMedium  == 'air':
                    ri = riValues[0]
                if (sensorMedium  == 'h20') | (sensorMedium  == 'h2o'):
                    ri = riValues[1]
                if sensorMedium  == '10%':
                    ri = riValues[2]
                if sensorMedium  == '20%':
                    ri = riValues[3]
                if sensorMedium  == '30%':
                    ri = riValues[4]
                if sensorMedium  == '40%':
                    ri = riValues[5]
                if sensorMedium  == '50%':
                    ri = riValues[6]
    #            colorDataMean.append({'filename': fileName,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'R':ParameterStats[0,0,recNum],'G':ParameterStats[1,0,recNum],'B':ParameterStats[2,0,recNum],'H':ParameterStats[3,0,recNum],'S':ParameterStats[4,0,recNum],'V':ParameterStats[5,0,recNum],'L*':ParameterStats[6,0,recNum],'a*':ParameterStats[7,0,recNum],'b*':ParameterStats[8,0,recNum],'r':ParameterStats[9,0,recNum],'g':ParameterStats[10,0,recNum],'b':ParameterStats[11,0,recNum],'G/R':ParameterStats[12,0,recNum],'B/R':ParameterStats[13,0,recNum],'G/B':ParameterStats[14,0,recNum],'Rw':rfactor,'Gw':gfactor,'Bw':bfactor,'type':os.path.dirname(fileName)[os.path.dirname(fileName).find('Photos')+7:os.path.dirname(fileName).find('Photos')+10],'ri':os.path.dirname(fileName)[-3:],'slide':os.path.dirname(fileName)[os.path.dirname(fileName).lower().find('slide')+6:os.path.dirname(fileName).lower().find('slide')+7]})
    #            colorDataMost.append({'filename': fileName,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'R':ParameterStats[0,2,recNum],'G':ParameterStats[1,2,recNum],'B':ParameterStats[2,2,recNum],'H':ParameterStats[3,2,recNum],'S':ParameterStats[4,2,recNum],'V':ParameterStats[5,2,recNum],'L*':ParameterStats[6,2,recNum],'a*':ParameterStats[7,2,recNum],'b*':ParameterStats[8,2,recNum],'r':ParameterStats[9,2,recNum],'g':ParameterStats[10,2,recNum],'b':ParameterStats[11,2,recNum],'G/R':ParameterStats[12,2,recNum],'B/R':ParameterStats[13,2,recNum],'G/B':ParameterStats[14,2,recNum],'Rw':rfactor,'Gw':gfactor,'Bw':bfactor,'type':os.path.dirname(fileName)[os.path.dirname(fileName).find('Photos')+7:os.path.dirname(fileName).find('Photos')+10],'ri':os.path.dirname(fileName)[-3:],'slide':os.path.dirname(fileName)[os.path.dirname(fileName).lower().find('slide')+6:os.path.dirname(fileName).lower().find('slide')+7]})
                colorDataMean.append({'medium':sensorMedium,'filename': fileName,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'area':ParameterStats[15,0,recNum],'R':ParameterStats[0,0,recNum],'G':ParameterStats[1,0,recNum],'B':ParameterStats[2,0,recNum],'H':ParameterStats[3,0,recNum],'S':ParameterStats[4,0,recNum],'V':ParameterStats[5,0,recNum],'L*':ParameterStats[6,0,recNum],'a*':ParameterStats[7,0,recNum],'b*':ParameterStats[8,0,recNum],'r':ParameterStats[9,0,recNum],'g':ParameterStats[10,0,recNum],'b':ParameterStats[11,0,recNum],'Ra':ParameterStats[12,0,recNum],'Ga':ParameterStats[13,0,recNum],'Ba':ParameterStats[14,0,recNum],'Rw':rfactor,'Gw':gfactor,'Bw':bfactor,'Rl':Rr,'Gl':Gr,'Bl':Br,'type':gen,'ri':ri,'slide':slide,'prep':prep,'SAM':sam,'Surface':surface})
                colorDataMost.append({'medium':sensorMedium,'filename': fileName,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'area':ParameterStats[15,0,recNum],'R':ParameterStats[0,2,recNum],'G':ParameterStats[1,2,recNum],'B':ParameterStats[2,2,recNum],'H':ParameterStats[3,2,recNum],'S':ParameterStats[4,2,recNum],'V':ParameterStats[5,2,recNum],'L*':ParameterStats[6,2,recNum],'a*':ParameterStats[7,2,recNum],'b*':ParameterStats[8,2,recNum],'r':ParameterStats[9,2,recNum],'g':ParameterStats[10,2,recNum],'b':ParameterStats[11,2,recNum],'Ra':ParameterStats[12,2,recNum],'Ga':ParameterStats[13,2,recNum],'Ba':ParameterStats[14,2,recNum],'Rw':rfactor,'Gw':gfactor,'Bw':bfactor,'Rl':Rr,'Gl':Gr,'Bl':Br,'type':gen,'ri':ri,'slide':slide,'prep':prep,'SAM':sam,'Surface':surface})
                colorDataMeanOriginal.append({'medium':sensorMedium,'filename': fileName,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'area':ParameterStats[15,0,recNum],'R':ParameterStats[0,0,recNum+1],'G':ParameterStats[1,0,recNum+1],'B':ParameterStats[2,0,recNum+1],'H':ParameterStats[3,0,recNum+1],'S':ParameterStats[4,0,recNum+1],'V':ParameterStats[5,0,recNum+1],'L*':ParameterStats[6,0,recNum+1],'a*':ParameterStats[7,0,recNum+1],'b*':ParameterStats[8,0,recNum+1],'r':ParameterStats[9,0,recNum+1],'g':ParameterStats[10,0,recNum+1],'b':ParameterStats[11,0,recNum+1],'Ra':ParameterStats[12,0,recNum+1],'Ga':ParameterStats[13,0,recNum+1],'Ba':ParameterStats[14,0,recNum+1],'Rw':rfactor,'Gw':gfactor,'Bw':bfactor,'Rl':Rr,'Gl':Gr,'Bl':Br,'type':gen,'ri':ri,'slide':slide,'prep':prep,'SAM':sam,'Surface':surface})
                colorDataMostOriginal.append({'medium':sensorMedium,'filename': fileName,'cTime': cTime.strftime("%Y-%m-%d %H:%M:%S"),'area':ParameterStats[15,0,recNum],'R':ParameterStats[0,2,recNum+1],'G':ParameterStats[1,2,recNum+1],'B':ParameterStats[2,2,recNum+1],'H':ParameterStats[3,2,recNum+1],'S':ParameterStats[4,2,recNum+1],'V':ParameterStats[5,2,recNum+1],'L*':ParameterStats[6,2,recNum+1],'a*':ParameterStats[7,2,recNum+1],'b*':ParameterStats[8,2,recNum+1],'r':ParameterStats[9,2,recNum+1],'g':ParameterStats[10,2,recNum+1],'b':ParameterStats[11,2,recNum+1],'Ra':ParameterStats[12,2,recNum+1],'Ga':ParameterStats[13,2,recNum+1],'Ba':ParameterStats[14,2,recNum+1],'Rw':rfactor,'Gw':gfactor,'Bw':bfactor,'Rl':Rr,'Gl':Gr,'Bl':Br,'type':gen,'ri':ri,'slide':slide,'prep':prep,'SAM':sam,'Surface':surface})
                cv2.imwrite(os.path.dirname(fileName)+'/processed'+os.path.basename(fileName), img)
                cv2.imshow('Result', img)

                #scaleFrame=resFrame*255
                #scaleFrame=scaleFrame.astype(np.uint8)
                #cv2.imwrite(gen+'_'+slide+'_'+prep+'_'+sensorMedium+'.jpg', scaleFrame)

                recNum=recNum+1
        #cv2.destroyAllWindows()
if (mode != '1') | (mode != 1):
    now = datetime.datetime.now()
    writer = pd.ExcelWriter(startPath+'/DataSummary'+now.strftime("%m-%d-%Yat%H-%M")+'.xlsx')
    workbook  = writer.book

    if processSpectra: 
        dfSpectraData=pd.DataFrame(spectraData)
        dfSpectraData=dfSpectraData[['filename','cTime','prep','type','SAM','Surface','slide','medium','ri','gaussLambda','gaussSigma','Max Abs','smoothedDataLambda','derDataLambda','R','G','B','H','S','V','L*','a*','b*','r','g','b']]
        dfSpectraData.to_excel(writer,'Spectra_Values')
        dfSummarySpectraData=dfSpectraData.groupby(['type','ri','slide']).describe()
        #dfColorAllMean.to_excel(writer,'Color_Values_Mean')
        dfSummarySpectraData.to_excel(writer,'Spectra_Stats')
        dfSummaryMeanSpectraData=dfSpectraData.groupby(['type','ri','slide']).mean()
        #dfColorAllMean.to_excel(writer,'Color_Values_Mean')
        dfSummaryMeanSpectraData.to_excel(writer,'Spectra_Means')

    if processPhotos:
        dfColorMean=pd.DataFrame(colorDataMean)
        dfColorMean=dfColorMean[['filename','cTime','prep','type','SAM','Surface','slide','area','medium','ri','R','G','B','H','S','V','L*','a*','b*','r','g','b','Ra','Ga','Ba','Rw','Gw','Bw','Rl','Gl','Bl']]
        dfColorMost=pd.DataFrame(colorDataMost)
        dfColorMost=dfColorMost[['filename','cTime','prep','type','SAM','Surface','slide','area','medium','ri','R','G','B','H','S','V','L*','a*','b*','r','g','b','Ra','Ga','Ba','Rw','Gw','Bw','Rl','Gl','Bl']]
        colorDataMeanOriginal=pd.DataFrame(colorDataMeanOriginal)
        colorDataMeanOriginal=colorDataMeanOriginal[['filename','cTime','prep','type','SAM','Surface','slide','area','medium','ri','R','G','B','H','S','V','L*','a*','b*','r','g','b','Ra','Ga','Ba','Rw','Gw','Bw','Rl','Gl','Bl']]
        colorDataMostOriginal=pd.DataFrame(colorDataMostOriginal)
        colorDataMostOriginal=colorDataMostOriginal[['filename','cTime','prep','type','SAM','Surface','slide','area','medium','ri','R','G','B','H','S','V','L*','a*','b*','r','g','b','Ra','Ga','Ba','Rw','Gw','Bw','Rl','Gl','Bl']]
        dfColorMean.to_excel(writer,'Color_Values_Mean')
        dfColorMost.to_excel(writer,'Color_Values_Most')
        colorDataMeanOriginal.to_excel(writer,'Color_Values_Mean_Original')
        colorDataMostOriginal.to_excel(writer,'Color_Values_Most_Original')
    #writer.save()
        #dfColorAllMean = dfColorAllMean.append(dfColorMean)
        #dfColorAllMost = dfColorAllMost.append(dfColorMost)
    
        #writer = pd.ExcelWriter('StillPhotoSummaryAll.xlsx')
        dfSummaryMean=dfColorMean.groupby(['type','ri','slide']).describe()
        #dfColorAllMean.to_excel(writer,'Color_Values_Mean')
        dfSummaryMean.to_excel(writer,'Color_Stats_Mean')
        dfSummaryMost=dfColorMost.groupby(['type','ri','slide']).describe()
        #dfColorAllMost.to_excel(writer,'Color_Values_Most')
        dfSummaryMost.to_excel(writer,'Color_Stats_Most')
        dfSummaryMeanMean=dfColorMean.groupby(['type','ri','slide']).mean()
        #dfColorAllMean.to_excel(writer,'Color_Values_Mean')
        dfSummaryMeanMean.to_excel(writer,'Color_Mean_of_Means')

#if (mode == '3') | (mode == 3):
#    typeList=dfColorMean['type'].unique()
#    fitDataSepPhoto=[]
#    fitDataComboPhoto=[]
#    fitDataSepSpectrum=[]
#    fitDataComboSpectrum=[]
#    for generation in typeList:
#        figSlidesSepPhoto, axesSepPhoto = plt.subplots(5,3,sharey=False,sharex=True,figsize=(16,16))
#        figSlidesSepPhoto.set_facecolor((0.7,0.7,0.7))
#        activeColor=['black','red','green','blue','cyan','grey','black','white','magenta','yellow','salmon','lightseagreen','skyblue','gold','mediumorchid','darkcyan','darkred','darkgreen','darkblue']
#        activeChannel=['LambdaMax','Red','Green','Blue','Hue','Saturation','Value','CIE L*','CIE a*','CIE b*','r chromaticity','g chromaticity','b chromaticity','Ra','Ga','Ba','Red gamma','Green gamma','Blue gamma']
#        dfGenPhoto=dfColorMean[(dfColorMean['type']==generation) & (dfColorMean['ri']>1.3)]
#        slideList=dfGenPhoto['slide'].unique()
#        riList=dfGenPhoto['ri'].unique()
#        Xrange=np.arange(min(riList),max(riList)+0.001,0.001)
#        annotation_summary=[['' for i in range(3)] for j in range(5)]
#        if processPhotos: 
#            for slide in slideList:
#                if slide=='A':
#                    slideColor='xkcd:red'
#                elif slide=='B':
#                    slideColor='xkcd:orange'
#                elif slide=='C':
#                    slideColor='xkcd:yellow'
#                elif slide=='D':
#                    slideColor='xkcd:green'
#                elif slide=='E':
#                    slideColor='xkcd:blue'
#                elif slide=='F':
#                    slideColor='xkcd:violet'
#                else: 
#                    slideColor='xkcd:black'            
#                activeColumn=1
#                dfSlide=dfGenPhoto[dfGenPhoto['slide']==slide]
#                for cc in range(5):
#                    for feature in range(3):
#                        axesSepPhoto[cc,feature].set_facecolor((0.7,0.7,0.7))
#                        axesSepPhoto[cc,feature].scatter(dfSlide['ri'],dfSlide.iloc[:,activeColumn+7],alpha=0.2,color=slideColor)
#                        fit=PolyReg(dfSlide['ri'],dfSlide.iloc[:,activeColumn+7],1)
#                        axesSepPhoto[cc,feature].plot(Xrange,fit['poly'](Xrange),'-',color=slideColor)
#                        axesSepPhoto[cc,feature].set(ylabel=activeChannel[activeColumn])
#                        fitDataSepPhoto.append({'gen':generation,'slide':slide,'cc':activeChannel[activeColumn],'slope': fit['coef'][0],'intercept': fit['coef'][1],'error in slope':fit['errors'][0],'error in intercept':fit['errors'][1],'n':fit['n'],'standard error in y':fit['sy'],'sensitivity':np.abs(fit['coef'][0]/fit['sy'])})
#                        activeColumn=activeColumn+1
#                        annotation_summary[cc][feature]=annotation_summary[cc][feature]+'slide'+slide+' m={0:.1E}'.format(fit['coef'][0])+' sy={0:.1E}'.format(fit['sy'])+' sens={0:.1E}'.format(np.abs(fit['coef'][0]/fit['sy']))+'\n'
#    #        for cc in range(5):
#    #            for feature in range(3):
#    #                AnnotateFit(fit,axesSepPhoto[cc,feature],annotation_summary[cc][feature][:-1])
#            figSlidesSepPhoto.tight_layout()
#            figSlidesSepPhoto.savefig(os.path.join(startPath,generation+'_SepPhoto.png'))
#            
#            activeColumn=1
#            figSlidesComboPhoto, axesComboPhoto = plt.subplots(5,3,sharey=False,sharex=True,figsize=(16,16))
#            figSlidesComboPhoto.set_facecolor((0.7,0.7,0.7))
#            for cc in range(5):
#                for feature in range(3):
#                    axesComboPhoto[cc,feature].set_facecolor((0.7,0.7,0.7))
#                    axesComboPhoto[cc,feature].scatter(dfGenPhoto['ri'],dfGenPhoto.iloc[:,activeColumn+7],alpha=0.2)
#                    fitCombo=PolyReg(dfGenPhoto['ri'],dfGenPhoto.iloc[:,activeColumn+7],1)
#                    axesComboPhoto[cc,feature].plot(Xrange,fitCombo['poly'](Xrange),'-')
#                    AnnotateFit(fitCombo,axesComboPhoto[cc,feature],Arrow=False,xText=0.05,yText=0.50)
#                    axesComboPhoto[cc,feature].set(ylabel=activeChannel[activeColumn])
#                    fitDataComboPhoto.append({'gen':generation,'cc':activeChannel[activeColumn],'slope': fitCombo['coef'][0],'intercept': fitCombo['coef'][1],'error in slope':fitCombo['errors'][0],'error in intercept':fitCombo['errors'][1],'n':fitCombo['n'],'standard error in y':fitCombo['sy'],'sensitivity':np.abs(fitCombo['coef'][0]/fitCombo['sy'])})
#                    activeColumn=activeColumn+1
#            figSlidesComboPhoto.tight_layout()
#            figSlidesComboPhoto.savefig(os.path.join(startPath,generation+'_ComboPhoto.png'))
#        
#        if processSpectra: 
#            figSlidesSepSpectrum, axesSepSpectrum = plt.subplots(5,3,sharey=False,sharex=True,figsize=(16,16))
#            figSlidesSepSpectrum.set_facecolor((0.7,0.7,0.7))
#            activeColorSpectrum=['black','red','green','blue','cyan','grey','black','white','magenta','yellow','salmon','lightseagreen','skyblue','gold','mediumorchid','darkcyan','darkred','darkgreen','darkblue']
#            activeChannelSpectrum=['Gauss','Max','Derivative','Red','Green','Blue','Hue','Saturation','Value','CIE L*','CIE a*','CIE b*','r chromaticity','g chromaticity','b chromaticity','Ra','Ga','Ba','Red gamma','Green gamma','Blue gamma']
#            dfGenSpectrum=dfSpectraData[(dfSpectraData['type']==generation) & (dfSpectraData['ri']>1.3)]
#            slideList=dfGenSpectrum['slide'].unique()
#            riList=dfGenSpectrum['ri'].unique()
#            Xrange=np.arange(min(riList),max(riList)+0.001,0.001)
#            annotation_summary=[['' for i in range(3)] for j in range(5)]
#            for slide in slideList:
#                if slide=='A':
#                    slideColor='xkcd:red'
#                elif slide=='B':
#                    slideColor='xkcd:orange'
#                elif slide=='C':
#                    slideColor='xkcd:yellow'
#                elif slide=='D':
#                    slideColor='xkcd:green'
#                elif slide=='E':
#                    slideColor='xkcd:blue'
#                elif slide=='F':
#                    slideColor='xkcd:violet'
#                else: 
#                    slideColor='xkcd:black'   
#                activeColumn=0
#                dfSlide=dfGenSpectrum[dfGenSpectrum['slide']==slide]
#                waterLmax=dfSlide['gaussLambda'][dfSlide['medium']=='h20']
#                for cc in range(5):
#                    for feature in range(3):
#                        axesSepSpectrum[cc,feature].set_facecolor((0.7,0.7,0.7))
#                        axesSepSpectrum[cc,feature].scatter(dfSlide['ri'],dfSlide.iloc[:,activeColumn+7],alpha=0.2,color=slideColor)
#                        fit=PolyReg(dfSlide['ri'],dfSlide.iloc[:,activeColumn+7],1)
#                        axesSepSpectrum[cc,feature].plot(Xrange,fit['poly'](Xrange),'-',color=slideColor)
#                        axesSepSpectrum[cc,feature].set(ylabel=activeChannelSpectrum[activeColumn])
#                        fitDataSepSpectrum.append({'gen':generation,'slide':slide,'cc':activeChannelSpectrum[activeColumn],'slope': fit['coef'][0],'intercept': fit['coef'][1],'error in slope':fit['errors'][0],'error in intercept':fit['errors'][1],'n':fit['n'],'standard error in y':fit['sy'],'sensitivity':np.abs(fit['coef'][0]/fit['sy']),'lmax_water':waterLmax.values[0]})
#                        activeColumn=activeColumn+1
#                        annotation_summary[cc][feature]=annotation_summary[cc][feature]+'slide'+slide+' m={0:.1E}'.format(fit['coef'][0])+' sy={0:.1E}'.format(fit['sy'])+' sens={0:.1E}'.format(np.abs(fit['coef'][0]/fit['sy']))+'\n'
#    #        for cc in range(5):
#    #            for feature in range(3):
#    #                AnnotateFit(fit,axesSepSpectrum[cc,feature],annotation_summary[cc][feature][:-1])
#            figSlidesSepSpectrum.tight_layout()
#            figSlidesSepSpectrum.savefig(os.path.join(startPath,generation+'_SepSpectrum.png'))
#              
#
#        
#            activeColumn=0
#            figSlidesComboSpectrum, axesComboSpectrum = plt.subplots(5,3,sharey=False,sharex=True,figsize=(16,16))
#            figSlidesComboSpectrum.set_facecolor((0.7,0.7,0.7))
#            for cc in range(5):
#                for feature in range(3):
#                    axesComboSpectrum[cc,feature].set_facecolor((0.7,0.7,0.7))
#                    axesComboSpectrum[cc,feature].scatter(dfGenSpectrum['ri'],dfGenSpectrum.iloc[:,activeColumn+7],alpha=0.2)
#                    fitCombo=PolyReg(dfGenSpectrum['ri'],dfGenSpectrum.iloc[:,activeColumn+7],1)
#                    axesComboSpectrum[cc,feature].plot(Xrange,fitCombo['poly'](Xrange),'-')
#                    AnnotateFit(fitCombo,axesComboSpectrum[cc,feature],Arrow=False,xText=0.05,yText=0.50)
#                    axesComboSpectrum[cc,feature].set(ylabel=activeChannelSpectrum[activeColumn])
#                    fitDataComboSpectrum.append({'gen':generation,'cc':activeChannelSpectrum[activeColumn],'slope': fitCombo['coef'][0],'intercept': fitCombo['coef'][1],'error in slope':fitCombo['errors'][0],'error in intercept':fitCombo['errors'][1],'n':fitCombo['n'],'standard error in y':fitCombo['sy'],'sensitivity':np.abs(fitCombo['coef'][0]/fitCombo['sy'])})
#                    activeColumn=activeColumn+1
#            figSlidesComboSpectrum.tight_layout()    
#            figSlidesComboSpectrum.savefig(os.path.join(startPath,generation+'_ComboSpectrum.png'))
#    
#    if processPhotos: 
#        dfFitsSepPhoto=pd.DataFrame.from_dict(fitDataSepPhoto)
#        dfFitsComboPhoto=pd.DataFrame.from_dict(fitDataComboPhoto)  
#        dfFitsSepPhoto=dfFitsSepPhoto.sort_values(by=['cc', 'gen','slide'])
#        dfFitsComboPhoto=dfFitsComboPhoto.sort_values(by=['cc', 'gen'])
#        dfFitsSepPhoto.to_excel(writer,'Photo_Fits_Separated')
#        dfFitsComboPhoto.to_excel(writer,'Photo_Fits_Combined')
#    if processSpectra: 
#        dfFitsSepSpectrum=pd.DataFrame.from_dict(fitDataSepSpectrum)
#        dfFitsComboSpectrum=pd.DataFrame.from_dict(fitDataComboSpectrum)
#        dfFitsSepSpectrum=dfFitsSepSpectrum.sort_values(by=['cc', 'gen','slide'])
#        dfFitsComboSpectrum=dfFitsComboSpectrum.sort_values(by=['cc', 'gen'])
#        dfFitsSepSpectrum.to_excel(writer,'Spectrum_Fits_Separated')
#        dfFitsComboSpectrum.to_excel(writer,'Spectrum_Fits_Combined')
    
writer.save()
