"""
This utility module contains all common imports needed, and all basic functions used by the other modules.
"""

# basic:
import re
import math
import os
import numpy as np
import numpy.ma as ma

# other:
from scipy.ndimage import rotate
from astroquery.mast import Tesscut
import io
import sys
from contextlib import contextmanager

# matplotlib:
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["savefig.dpi"] = 150
rcParams["figure.dpi"] = 150
rcParams["font.size"] = 8
rcParams['savefig.bbox'] = 'tight'
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

# astropy
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.stats import sigma_clip # function
from astropy.stats import SigmaClip # class

from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.visualization import ZScaleInterval

from astropy.nddata.utils import Cutout2D
from astropy.modeling.models import Ellipse2D

# photutils:
from photutils.background import Background2D, SExtractorBackground
from photutils import Background2D, MedianBackground
from photutils.isophote import build_ellipse_model
from photutils.aperture import EllipticalAperture
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse

""" Warning suppression: prevents warnings that arise from JWST images, to clean up printed output.
Here is an example warning that is suppressed:
WARNING: FITSFixedWarning: 'datfix' made the change 'Set DATE-BEG to '2022-07-02T14:38:52.612' from MJD-BEG.
Set DATE-AVG to '2022-07-02T14:48:10.988' from MJD-AVG.
Set DATE-END to '2022-07-02T14:57:29.338' from MJD-END'. [astropy.wcs.wcs] 
"""
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


@contextmanager
def suppress_stdout():
    """function to supress print ( I don't think I ended up using this?)"""
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout



def quick_open(data, sci=0, plot=True):
    """
    Quickly imports and plots a .fits file to view the data.
    
    Parameters:
        data (string): Path to .fits file.
    Returns:
        hdulist0, image_data0, prihdr0: hdu object, data, and header of opened file.
    """
    # open and get data and header info
    hdulist0 = fits.open(data)
    print(hdulist0.info())
    prihdr0 = hdulist0[sci].header
    image_data0 = hdulist0[sci].data

    # find if the data is 2D or 3D. If 3D, work with the zero index of the image data
    if len(image_data0.shape) > 2:
        print('image data is 3D, using 2D axis')
        image_data0 = image_data0[0]
    else:
        print('image is 2D.')

    # uses the file name from the path for the plot title
    name = str(data).split('/')[-1]

    if plot==True:
        quick_plot(image_data0, name=name)

    return hdulist0, image_data0, prihdr0


def quick_plot(data, name='Plot'):
    """
    Quickly plots a data array with both the raw and ZScaled data.

    Parameters:
        data (numpy array): image data array (called image_data0 or image_data elsewhere).
        name (string): the name of the data given to the plot.
    """
    # plot
    fig = plt.figure(figsize=(12,6), dpi=150)  

    imin, imax = np.min(data), np.max(data)
    plt.subplot(1,3,1)
    plt.imshow(data, origin="lower", vmin=imin, vmax=imax, cmap="gray")
    plt.title(name+': raw')
    plt.colorbar(shrink=.45)

    plt.subplot(1,3,2)
    plt.imshow(data, origin="lower", vmin=imin, vmax=imax, cmap="gray", norm='asinh')
    plt.title(name+': asinh')
    plt.colorbar(shrink=.45)

    imin,imax = apply_zscale(data)
    plt.subplot(1,3,3)
    plt.imshow(data, origin="lower", vmin=imin, vmax=imax, cmap='viridis')
    plt.title(name+': ZScale')
    plt.colorbar(shrink=.45)
    
    plt.tight_layout()


def create_working_file(path, file, pos, size, cutout=False, sci=0):
    """
    Creates a working file for convenience and modification. If needed, it will also create a smaller cutout of the data.

    Parameters:
        path (string): Path from the project to original .fits file.
        file (string): Path from the project to the working and save location.
    Return:
    """
    # opens data
    hdulist0, image_data0, prihdr0 = quick_open(path, sci, plot=False)

    wcs = WCS(prihdr0)

    # creates cutout if needed
    try:
        if cutout == True:
            cutout = Cutout2D(image_data0, pos, size, wcs=wcs)
            save_cutout(path, cutout, file+'.fits', sci) # writes to a new file
        else:
            cutout = image_data0
            hdulist = hdulist0
            hdulist.writeto(file+'.fits', overwrite=True)
    except:
        print('Cannot make a smaller cutout. Perhaps it already exists, and is being used by another process.')
        print('If a new cutout is desired, close ds9 and restart kernal to update.')
        print('Continuing to next step, using previous cutout if available...')
    
    # opens newly created file, and closes old one.
    # redefine variables: variables with 0 at the end indicate the original file. All others will refer to the new (cutout) file.
    hdulist0.close()
    hdulist = fits.open(file+'.fits')
    prihdr = hdulist[sci].header
    image_data = hdulist[sci].data
    print(f'data shape: {image_data.shape}')

    return hdulist, image_data, prihdr



def apply_zscale(data):
    """
    A short cut to find ZScale interval of data.

    Paramters:
        data (numpy array): Two-dimensional image data.
    Returns:
        imin, imax : minimum and maximum values of ZScaled data.
    """
    interval = ZScaleInterval()
    (imin,imax) = interval.get_limits(data)
    return imin,imax


def fix_zeros(data):
    """Fixes divide by zero errors by replacing zero with 1e-5."""
    data[data == 0.] = 1e-5
    return data


def drawregions(fitsfile, science=0):
    """
    Open .fits file in DS9 and upload.

    Parameters:
        science (int): the index of the science image of the fits file, default is 0.
    """
    tmp = fits.open(fitsfile)[science]
    print(tmp)
    print(f'data shape:\t {np.shape(tmp.data)}\n')

    if len(np.shape(tmp.data)) == 3:
        ymax, xmax, nimages = np.shape(tmp.data)
    elif len(np.shape(tmp.data)) == 2:
        ymax, xmax = np.shape(tmp.data)

    xc = math.floor(xmax/2)
    yc = math.floor(ymax/2)

    regfile = open(fitsfile.replace('fits','reg'), mode='w')
    regfile.write('# Region file format: DS9 version 4.1\n')
    regfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    regfile.write('image\n')
    regfile.write('point(%8.4f,%8.4f) # point=x 15 color=red text={BG center}\n' % (xc+100,yc+100))
    regfile.write('point(%8.4f,%8.4f) # point=cross 15 color=blue text={FG center}\n' % (xc-100,yc-100))
    regfile.write('ellipse(%8.4f,%8.4f,300,150,90.0) # color=blue\n' % (xc+100,yc+100))
    regfile.write('ellipse(%8.4f,%8.4f,300,150,0.0) # color=red\n' % (xc-100,yc-100))
    regfile.close()

    print('Please open DS9 with the newly generated .reg file and move')
    print('the x to the forground galaxy center, ')
    print('the + to the background galaxy center, ')
    print('the blue ellipse to encompass the foreground galaxy, and')    
    print(f'the red ellipse to encompass the background galaxy.\n')
    print('WHEN SAVING, OPT FOR IMAGE so everything is in PIXELS for TheOG to deal with.')
    print('Best practice is to make a separate .reg file in WCS.')
    

    #stream = os.popen('pwd')
    #pwd = stream.read()
    #print(output)

    #os.system('ls -l')
    #os.system('pwd')
    #os.system('open ds9 %s -regionfile %s & ' % ( pwd+'/'+fitsfile, pwd+'/'+fitsfile.replace('fits','reg') ) )



def maskobjects(fitsfile):
    """
    open fitsfile in DS9 and upload 
    """
    tmp = fits.open(fitsfile)[1]
    print(tmp)
    print(np.shape(tmp.data))

    if len(np.shape(tmp.data)) == 3:
        ymax, xmax, nimages = np.shape(tmp.data)
    elif len(np.shape(tmp.data)) == 2:
        ymax, xmax = np.shape(tmp.data)

    xc = math.floor(xmax/2)
    yc = math.floor(ymax/2)

    regfile = open(fitsfile.replace('.fits','_mask_objects.reg'), mode='w')
    regfile.write('# Region file format: DS9 version 4.1\n')
    regfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    regfile.write('image\n')
    regfile.write('ellipse(%8.4f,%8.4f,300,150,90.0) # color=magenta\n' % (xc+100,yc+100))
    regfile.close()
    print('Please open DS9 with the newly generated .reg file')
    print('Copy and paste the ellipse as often as needed')
    print('Move and adjust each ellipse onto objects that need to be masked.')

    print('WHEN SAVING, OPT FOR IMAGE so everything is in PIXELS for TheOG to deal with.')
    print('Best practice is to make a separate .reg file in WCS.')
    

def readmaskregions(regfile):
    masks = []
    reg = open(regfile, mode='r')
    lines = reg.readlines()
    reg.close()
    for line in lines:
        if line.count('ellipse') > 0 and  line.count('orange') > 0:
            tmp2 = tmp[1].split(',')
            x0,y0,a,b,pa = float(tmp2[0]),float(tmp2[1]),float(tmp2[2]),float(tmp2[3]),float(tmp2[4])
            masks.append([x0,y0,a,b,pa])
    return masks
    
def readregion(regfile):
    """
    read in the .ref file provided
    """

    x0_fg,y0_fg,a_fg,b_fg,pa_fg, x0_bg,y0_bg,a_bg,b_bg,pa_bg = 0,0,0,0,0, 0,0,0,0,0

    # re.split('; |, ',str)
    
    reg = open(regfile, mode='r')
    lines = reg.readlines()
    reg.close()
    for line in lines:
        print(line)
        if line.count('ellipse') > 0 and  line.count('blue') > 0:
            print('Found the foreground galaxy ellipse!')
            tmp = re.split(r"\(|\) ",line) 
            tmp2 = tmp[1].split(',')
            x0_fg_ell,y0_fg_ell,a_fg,b_fg,pa_fg = float(tmp2[0]),float(tmp2[1]),float(tmp2[2]),float(tmp2[3]),float(tmp2[4])
            print(x0_fg_ell,y0_fg_ell,a_fg,b_fg,pa_fg,'\n')

        if line.count('ellipse') > 0 and  line.count('red') > 0:
            print('Found the background galaxy ellipse!')
            tmp = re.split(r"\(|\) ",line) 
            tmp2 = tmp[1].split(',')
            x0_bg_ell,y0_bg_ell,a_bg,b_bg,pa_bg = float(tmp2[0]),float(tmp2[1]),float(tmp2[2]),float(tmp2[3]),float(tmp2[4])
            print(x0_bg_ell,y0_bg_ell,a_bg,b_bg,pa_bg,'\n')

        if line.count('point') > 0 and  line.count('red') > 0:
            print('Found the background galaxy center marker!')
            tmp = re.split(r"\(|\) ",line) 
            tmp2 = tmp[1].split(',')
            x0_bg,y0_bg = float(tmp2[0]),float(tmp2[1])
            print(x0_bg,y0_bg)

        if line.count('point') > 0 and  line.count('blue') > 0:
            print('Found the foreground galaxy center marker!')
            tmp = re.split(r"\(|\) ",line) 
            tmp2 = tmp[1].split(',')
            x0_fg,y0_fg = float(tmp2[0]),float(tmp2[1])
            print(x0_fg,y0_fg)
    if x0_fg == 0:
        x0_fg = x0_fg_ell
    if y0_fg == 0:
        y0_fg = y0_fg_ell
    if x0_bg == 0:
        x0_bg = x0_bg_ell
    if y0_bg == 0:
        y0_bg = y0_bg_ell

    return x0_fg,y0_fg,a_fg,b_fg,pa_fg, x0_bg,y0_bg,a_bg,b_bg,pa_bg 




def save_cutout(originalfile, cutout, outfile, science=0):
    # Note: There is also a save_cutout_2 function

    # Load the image and the WCS of the original image
    hdu = fits.open(originalfile)[science] #input science as the index of the image to look at
    wcs = WCS(hdu.header)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    hdu.writeto(outfile, overwrite=True)


def normalize(arr, t_min=0, t_max=1):
    """
    Normalizes an array. Used in apply_mask.
    """
    norm_arr = []
    diff = t_max - t_min
    imin,imax = np.min(arr), np.max(arr)
    diff_arr = imax - imin    
    for i in arr:
        temp = (((i - imin)*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return norm_arr


def make_ellipse_mask(data,aper):
    """
    Make a mask of an ellipse in the data (from original code)
    """
    (x0,y0,a,b,pa) = aper

    theta = Angle(pa, 'deg')
    ell_mask = Ellipse2D(amplitude=1., x_0=x0, y_0=y0, a=a, b=b,
                    theta=theta.radian)
    ymax, xmax  = np.shape(data)
    y, x = np.mgrid[0:ymax, 0:xmax]

    return ell_mask(x,y)


def apply_mask(data,aper):
    """
    Makes and applies the mask to a data array. Get only the data within an ellipse.
    """
    data_norm = normalize(data,0,1)         # normalize
    e = make_ellipse_mask(data_norm,aper)   # make mask of ellipse
    e = np.abs(e-1)                         # invert mask
    masked = data_norm - e                  # apply mask
    masked[masked < 0.] = 0                 # remove less than 0
    
    return masked


def plot_ellipses(data,header,fgaper,bgaper,figname,vmin=0.0001,vmax=1, norm='linear', ZScale=False, label=False):

    if ZScale is True:
        vmin,vmax = apply_zscale(data)

    # unpack the ellipse information for foreground and background apertures
    (x0_fg,y0_fg,a_fg,b_fg,pa_fg) = fgaper
    (x0_bg,y0_bg,a_bg,b_bg,pa_bg) = bgaper

    geometry_fg = EllipseGeometry(x0=x0_fg, y0=y0_fg, sma=a_fg, eps=(1. - b_fg/a_fg), pa=pa_fg * np.pi / 180.)
    aper_fg = EllipticalAperture((geometry_fg.x0, geometry_fg.y0), geometry_fg.sma, geometry_fg.sma * (1 - geometry_fg.eps), geometry_fg.pa)

    geometry_bg = EllipseGeometry(x0=x0_bg, y0=y0_bg, sma=a_bg, eps=(1. - b_bg/a_bg),pa=pa_bg * np.pi / 180.)
    aper_bg = EllipticalAperture((geometry_bg.x0, geometry_bg.y0), geometry_bg.sma,geometry_bg.sma * (1 - geometry_bg.eps),geometry_bg.pa)

    fig = plt.figure(dpi=150)
    ax1 = plt.subplot(1,1,1, projection=WCS(header).celestial)
    cbar = ax1.imshow(data, origin='lower', vmin=vmin,vmax=vmax,norm=norm,cmap='viridis')
    aper_fg.plot(color='blue')
    aper_bg.plot(color='red')

    if label==True:
        # adds labels to each region
        ax1.text(x0_bg-a_bg/3,y0_bg-a_bg/3,"B'", color='black')
        ax1.text(x0_fg+a_fg/3,y0_fg+a_fg/3,"F'", color='black')
        ax1.text(x0_fg-a_fg/2,y0_fg-a_fg/2,"I", color='black')
        figname = figname[0:-4]+'_label.png'

    ax1.coords['ra'].set_axislabel('Right Ascension')
    ax1.coords['dec'].set_axislabel('Declination')
    plt.title('Ellipse Apertures')
    # fig.colorbar(cbar)
    fig.savefig(figname)