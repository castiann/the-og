"""
This module contains all functions needed to analyze overlapping galaxies using the rotation technique.
The module utils contains most other imports, and access to utility functions.
"""
from . utils import *


# image rotation around an RA DEC
def imrotate(data,header, x0, y0, dimen=u.Quantity((2, 2), u.arcmin), rotangle=180.,  figname='', logplot=False, vmin=0.0001,vmax=1, ZScale=False):
    """
    Starting from an image and x,y positions of one of the galaxies,
    rotate around that position and return the rotated image
    Plot the original, rotated imate and difference in threeplot if figname is not empty.
    """
    if ZScale is True:
        vmin,vmax = apply_zscale(data)

    # rotate does not like nan as an input. Replacing nan with 0
    data[~np.isfinite(data)] = 0.

    # calculate the central position of the galaxy in pixels in RA/DEC
    c = pixel_to_skycoord(x0, y0, WCS(header).celestial)
    print() #just for spacing
    print(f'input position (x,y) in px:\t {x0}, {y0}')
    print(f'input position (ra,dec) in deg:\t {c.ra.deg}, {c.dec.deg}\n')

    # cut out size of the postage stamp (bigger than both galaxies, not too big but can be bigger than the image). Fills the rest to 0.
    cutout = Cutout2D(data, position=c, size=dimen, mode='partial', wcs=WCS(header).celestial,fill_value=0.)
    
    # print out some useful info
    print(f'original position:\t {cutout.position_original}')
    print(f'original shape:\t\t {data.shape}')
    # print(f'cutout position:\t {cutout.position_cutout}')
    print(f'cutout shape:\t\t {cutout.shape}')

    # use scipy to rotate the cutout
    rotated = rotate(cutout.data, rotangle, axes=(1,0), cval=0.0)

    # get ymax and xmax of the original data
    ymax, xmax = np.shape(data)
    dimen=(ymax, xmax) # now back in pixels and for the size of the data that was given in.

    # get corners of cutout
    xc0,yc0 = cutout.position_original
    dy, dx = cutout.shape
    dx = np.floor(dx/2)         # use floor to work in full pixels
    dy = np.floor(dy/2)
    # print(f'cutout position (center):\t\t ({xc0}, {yc0})')

    dymax, dxmax = np.shape(rotated.data)
    
    xcmin = int(xc0-dx)
    xcmax = int(xc0+dx)
    ycmin = int(yc0-dy)
    ycmax = int(yc0+dy)
    print(f'cutout corners:\t\t x:({xcmin}, {xcmax}), y:({ycmin}, {ycmax}) \n')

    # fixes potential rounding error caused by floor. Adds the extra pixel back if necessary
    if xcmax - xcmin < dxmax:
        xcmax +=1
    if ycmax - ycmin < dymax:
        ycmax +=1
    
    # Get the center of the data array, and create a recut of the rotated data, but within the bounds of the original data
    cdata = pixel_to_skycoord(xmax/2, ymax/2, WCS(header).celestial)
    recut = Cutout2D(rotated.data, position=cdata, size=dimen, mode='partial', wcs=cutout.wcs,fill_value=0.)
    recut = recut.data
 
    # open a new hdu for the rotated image. Note: fits.PrimaryHDU converts data and header into HDU
    hdu = fits.PrimaryHDU(recut)
    hdu.header = header
    
    # puts the new cutout data into plots
    if figname !='' and logplot==False:
        PlotThree(im1=data,im2=recut.data,im3=data-hdu.data,header=header,
              title1='Original Image', title2='Rotated Image', title3='Difference Image', figname=figname,
              vmin=vmin,vmax=vmax,cmap='viridis')
    elif figname !='' and logplot==True:
        PlotThree(im1=np.log10(data),im2=np.log10(hdu.data),im3=np.log10(data-hdu.data),header=header,
              title1='Original Image', title2='Rotated Image', title3='Difference Image', figname=figname,
              vmin=vmin,vmax=vmax,cmap='viridis')
    return hdu



def elip_bbox(fgaper,bgaper):
    """
    Draws a bounding box around an ellipse. Used for jiggle and rotate to locate ideal center.
    """
    # unpack the ellipse information for foreground and background apertures
    (x0_fg,y0_fg,a_fg,b_fg,pa_fg) = fgaper
    (x0_bg,y0_bg,a_bg,b_bg,pa_bg) = bgaper

    # geometry
    x1 = a_fg*np.cos(pa_fg)
    y1 = a_fg*np.sin(pa_fg)
    x2 = b_fg*np.cos(pa_fg+np.pi/2)
    y2 = -a_fg*np.sin(pa_fg+np.pi/2)
    width_fg = (x1**2+x2**2)**(1/2)*2
    height_fg = (y1**2+y2**2)**(1/2)*2
    len_fg = max(width_fg,height_fg) # to make the bounding box a square, for simplification

    x1 = a_bg*np.cos(pa_bg)
    y1 = a_bg*np.sin(pa_bg)
    x2 = b_bg*np.cos(pa_bg+np.pi/2)
    y2 = a_bg*np.sin(pa_bg+np.pi/2)
    width_bg = (x1**2+x2**2)**(1/2)*2
    height_bg = (y1**2+y2**2)**(1/2)*2
    len_bg = max(width_bg, height_bg)

    return len_fg, len_bg


def jiggle_and_rotate(data,header, x0, y0, pixrange=10, pixelstep=2, figname='', regfile=None, foreground=None, finetune=False):
    """
    Move the center of the rotation around and minimize the difference between original and rotated image. 
    Return the best center to rotate around and the rotated image.

    NOTE: We should be minimizing the difference within the ellipse, not the entire image

    note: best done within the ellipse area
    pixrange: range where the best center is located
    pixelstep: how many steps within the range that should be checked

    If you wish to use the automatic mode using ds9 regions, specify the regfile. MAYBE ADD MORE INFO HERE LATER
    """
    
    # IF USING AUTOMATIC, CODE STARTS FROM HERE
    if regfile != None:
        with suppress_stdout(): # supresses print in readregion
            x0_fg,y0_fg,a_fg,b_fg,pa_fg, x0_bg,y0_bg,a_bg,b_bg,pa_bg = readregion(regfile)
            # pack the ellipse information for foreground and background apertures
            fgaper = (x0_fg,y0_fg,a_fg,b_fg,pa_fg)
            bgaper = (x0_bg,y0_bg,a_bg,b_bg,pa_bg)
        len_fg, len_bg = elip_bbox(fgaper,bgaper)
        # len_fg, len_bg = x0_fg*2, y0_fg*2

        if foreground == True:
            pixrange = int(len_fg)
            aper = fgaper
        elif foreground == False:
            pixrange = int(len_bg)
            aper = bgaper
        else:
            raise Exception('Please indicate if this is the foreground or the background using foreground=True or foreground=False')
        pixelstep = int(pixrange/5)

    if finetune == True:
        pixrange = int(pixrange/5)
        pixelstep = int(pixrange/5)

    if pixelstep == 0:  # force a minimum pixelstep of 1
        pixelstep = 1

    print(f'pixrange = {pixrange}')
    print(f'pixelstep = {pixelstep}')

    # IF USING MANUAL, CODE STARTS FROM HERE

    masked_data = apply_mask(data,aper)

    diff_min = np.sum(np.abs(masked_data))
    best_dx = 0.
    best_dy = 0.

    for dx in np.arange(-pixrange,pixrange, pixelstep):
        for dy in np.arange(-pixrange,pixrange, pixelstep):
            with suppress_stdout(): # supresses print in imrotate
                rotim = imrotate(data,header, x0+dx, y0+dy, rotangle=180.)

            rotim = apply_mask(rotim.data,aper)
            diff = np.sum(np.abs(masked_data - rotim.data))
            print(f'checking {dx},{dy}')
            if diff < diff_min:
                diff_min = diff
                print(f'new best center: difference {diff}')
                best_dx = dx
                best_dy = dy

    imrotate(data,header,x0+best_dx,y0+best_dy,rotangle=180.,figname=figname) # maybe label the best center on plot

    return best_dx, best_dy





def mkobjectsmask(data,regions):
    """
    mask all the object in the regions list. 


    """

    ymax, xmax  = np.shape(data)
    y, x = np.mgrid[0:ymax, 0:xmax]
    [x0,y0,a,b,pa] = regions[0]
    theta = Angle(pa, 'deg')
    ell_mask = Ellipse2D(amplitude=1., x_0=x0, y_0=y0, a=a, b=b,theta=theta.radian)
    mask = ell_mask(x,y)

    for region in regions:
    
        [x0,y0,a,b,pa] = region

        theta = Angle(pa, 'deg')
        ell_mask = Ellipse2D(amplitude=1., x_0=x0, y_0=y0, a=a, b=b,theta=theta.radian)

    return ell_mask(x,y)



# model galaxy with ellipse  task.


def sigmaclip_mask(data, aper, sigmalim, figname = ''):
    """
    Sigma-clip the data

    sigma clipping the residual image, making sure one just counts the area within the aperture
    Clipping agressively as the tail end of the residual contains the signal of attenuation.
    Return a mask with the pixels that still can be used. 


    Parameters
    ----------
    data : numpy array
         Two-dimensional image data

    aper : elliptical aperture structure (x0,y0,a,b,pa) 
         parameters used to define an elliptical apertutre: central pixel coordinates (X0,Y0), 
         major axis (a), minor axis (b) and position angle (pa) on the sky for galaxy aperture.

    sigmalim : float
         The threshold for clipping used in the sigma-clipping. 

    figname : string (optional)
         The path and filename for the figure output f the histogram of pixel value pre and post
         sigma clipping.


    Returns
    -------
    mask : numpy array
        Two-dimensional array of the same dimensions as data with 1 for pixels within the aperture 
        that satisfy the sigma clipped threshold and 0 for pixels outside the aperture or iteratively 
        clipped from the pixel collection.
    
    """
    
    (x0,y0,a,b,pa) = aper

    
    ymax, xmax  = np.shape(data)
    y, x = np.mgrid[0:ymax, 0:xmax]
    theta = Angle(pa, 'deg')

    # deep copy the data to save the mask in
    mask = np.copy(data)
    # defining the galaxy area
    ell = Ellipse2D(amplitude=1., x_0=x0, y_0=y0, a=a, b=b,theta=theta.radian)
    # mask
    # ell_bg(x, y)
    filtered_data = sigma_clip(data[np.where(ell(x, y)>0)], sigma=sigmalim, maxiters=None,
                            cenfunc='mean', masked=True, copy=True)
    print(filtered_data)
    mask[np.where(ell(x, y)==0)] = 0.
    mask[np.where(ell(x, y)>0)][filtered_data==True] = 1.
    mask[np.where(ell(x, y)>0)][filtered_data==False] = 0.
    
    if figname !='':
        plt.hist(np.ravel(data[np.where(ell_bg(x, y)>0)]), bins=100,range=[-0.1,0.1],alpha=0.4,density=True,label='data');
        plt.hist(np.ravel(filtered_data), bins=100,range=[-0.1,0.1],alpha=0.4,density=True,label='sigma clipped');
        plt.xlabel('')
        plt.legend()
    
    # return np.std(np.ravel(data[np.where(ell(x, y)>0)])),np.std(np.ravel(filtered_data))
    return mask
    
def sky_estimate(data, mask=''):
    """
    estimate the sky in the unmasked area
    """

    sigma_clip = SigmaClip(sigma=3., maxiters=10)
    bkg_estimator = MedianBackground()
    if mask !='':
       bkg = Background2D(data, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    else: 
       bkg = Background2D(data[mask], (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return bkg

# cut out ellipse
# lookslike this can be done with regions package.


def inellipse(x, y, h, k, a, b):
    z = ((x-h)**2)/a**2 + ((y-k)**2)/b**2
    if z < 1: #the equation is ^ this expression == 1
        return True
    else:
        return False



    
def PlotThree(im1,im2,im3,header,
              title1='', title2='', title3='', figname='plotthree.png',
              vmin=0.0001,vmax=0.1,cmap='viridis'):
    """
    Three plots side-by-side
    """
    
    use_wcs = False

    # Figure showing the results
    fig = plt.figure(figsize=(12, 6), dpi=150)

    # if use_wcs == True:
    #     ax1 = plt.subplot(1,3,1, projection=WCS(header).celestial)
    #     ax1.coords['ra'].set_axislabel('Right Ascension')
    #     ax1.coords['dec'].set_axislabel('Declination')
    # else:
    #     ax1 = plt.subplot(1,3,1)
    # ax1.imshow(im1, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
    # ax1.set_title(title1)

    ax1 = plt.subplot(1,3,1, 
                    #   projection=WCS(header).celestial
                      )
    ax1.imshow(im1, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
    # ax1.coords['ra'].set_axislabel('Right Ascension')
    # ax1.coords['dec'].set_axislabel('Declination')
    ax1.set_title(title1)

    ax2 = plt.subplot(1,3,2, 
                    #   projection=WCS(header).celestial
                      )
    ax2.imshow(im2, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
    # ax2.coords['ra'].set_axislabel('Right Ascension')
    # ax2.coords['dec'].set_axislabel('Declination')
    ax2.set_title(title2)

    ax3 = plt.subplot(1,3,3, 
                    #   projection=WCS(header).celestial
                      )
    ax3.imshow(im3, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
    # ax3.coords['ra'].set_axislabel('Right Ascension')
    # ax3.coords['dec'].set_axislabel('Declination')
    ax3.set_title(title3)

    # fig.subplots_adjust(wspace=0.9)
    fig.savefig(figname)



def PlotSix(im1,im2,im3,im4,im5,im6,header,
              title1='', title2='', title3='', title4='', title5='', title6='', fgaper='',bgaper='', figname='plotsix.png',
              vmin='',vmax='',norm='linear',ZScale=False, celestial=False):
    """
    Organize six plots into two rows of three plots side-by-side. 

    Convenience function to plot six images in a 2x3 grid. For example, original data, rotated around the foreground galaxy 
    center, rotated around the background galaxy center, the uncertainty image in transmission, the transmission image and 
    the signal-to-noise ratio image (Edit: changing function to plot specifically the data mentioned above)

    Parameters
    ----------
    im1,im2,im3,im4,im5,im6 : numpy array
        The data in each image image to be displayed from top left to bottom right in order.

    title1,title2,title3,title4,title5,title6 : string (optional)
        The title of each image to be displayed above each.

    figname : string (default "plotsix")
        Path and filename of the output file for the figure. 

    vmin : float
        The minimum value in the color stretch of each image.
    vmax : 

    cmap
    

    Returns
    -------
    None, file is written to figname file path.

    """

    # unpack the ellipse information for foreground and background apertures
    if fgaper !='':
        (x0_fg,y0_fg,a_fg,b_fg,pa_fg) = fgaper
        geometry_fg = EllipseGeometry(x0=x0_fg, y0=y0_fg, sma=a_fg, eps=(1. - b_fg/a_fg), pa=pa_fg * np.pi / 180.)
        aper_fg = EllipticalAperture((geometry_fg.x0, geometry_fg.y0), geometry_fg.sma,geometry_fg.sma * (1 - geometry_fg.eps),geometry_fg.pa)
    if bgaper !='':
        (x0_bg,y0_bg,a_bg,b_bg,pa_bg) = bgaper
        geometry_bg = EllipseGeometry(x0=x0_bg, y0=y0_bg, sma=a_bg, eps=(1. - b_bg/a_bg),pa=pa_bg * np.pi / 180.)
        aper_bg = EllipticalAperture((geometry_bg.x0, geometry_bg.y0), geometry_bg.sma,geometry_bg.sma * (1 - geometry_bg.eps),geometry_bg.pa)

    # uncertainty should always use zscale, or it won't be visible
    vmin5,vmax5 = apply_zscale(im5)

    if ZScale is True:
        #uses the same scale from im1 to be consistent
        vmin1,vmax1 = apply_zscale(im1)
        vmin2,vmax2 = apply_zscale(im1)
        vmin3,vmax3 = apply_zscale(im1)
    else:
        vmin1 = np.min(im1)
        vmin2 = np.min(im2)
        vmin3 = np.min(im3)
        vmax1 = np.max(im1)
        vmax2 = np.max(im2)
        vmax3 = np.max(im3)

    if vmin =='':
        vmin4,vmin6=0,0
    else:
        vmin4,vmin6=vmin,vmin
    if vmax =='':
        vmax4,vmax6=1,1
    else:
        vmax4,vmax6=vmax,vmax

        
    # Figure showing the results
    fig = plt.figure(figsize=(12,6), dpi=150)

    if celestial == True:
        cel = WCS(header).celestial
    else:
        cel = None

    sh=1

    ax = plt.subplot(2,3,1, projection=cel)
    im = ax.imshow(im1, origin='lower', vmin=vmin1,vmax=vmax1,norm=norm)
    if celestial == True:
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
    else:
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')
    ax.set_title(title1)
    fig.colorbar(im, ax=ax, shrink=sh)

    ax = plt.subplot(2,3,2, projection=cel)
    im = ax.imshow(im2, origin='lower', vmin=vmin2,vmax=vmax2,norm=norm)
    if celestial == True:
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
    else:
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')
    ax.set_title(title2)
    fig.colorbar(im, ax=ax, shrink=sh)

    ax = plt.subplot(2,3,3, projection=cel)
    im = ax.imshow(im3, origin='lower', vmin=vmin3,vmax=vmax3,norm=norm)
    if celestial == True:
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
    else:
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')
    ax.set_title(title3)
    fig.colorbar(im, ax=ax, shrink=sh)
    
    ax = plt.subplot(2,3,4, projection=cel)
    im = ax.imshow(im4, origin='lower', vmin=vmin4,vmax=vmax4, norm=norm, interpolation='none', cmap='gray')
    if celestial == True:
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
    else:
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')
    ax.set_title(title4)
    fig.colorbar(im, ax=ax, shrink=sh)

    ax = plt.subplot(2,3,5, projection=cel)
    im = ax.imshow(im5, origin='lower', vmin=vmin5,vmax=vmax5, interpolation='none', cmap='gray')
    if celestial == True:
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
    else:
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')
    ax.set_title(title5, pad=15)
    cbar = fig.colorbar(im, ax=ax, shrink=sh)
    cbar.ax.yaxis.get_major_formatter().set_useOffset(True)
    cbar.ax.yaxis.get_major_formatter().set_scientific(True)
    cbar.ax.yaxis.get_offset_text().set_position((5, -20))
    cbar.update_ticks()

    ax = plt.subplot(2,3,6, projection=cel)
    im = ax.imshow(im6, origin='lower', vmin=vmin6,vmax=vmax6, norm=norm, interpolation='none', cmap='gray')
    if celestial == True:
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')
    else:
        plt.xlabel('Right Ascension')
        plt.ylabel('Declination')
    ax.set_title(title6)
    fig.colorbar(im, ax=ax, shrink=sh)

    # fig.subplots_adjust(wspace=0.5)
    # fig.subplots_adjust(hspace=0.1)
    # fig.subplots_adjust(hspace=.01)
    fig.tight_layout()
    fig.savefig(figname)





   

def background(data):
    # from photutils import Background2D, MedianBackground
    sigma_clip = SigmaClip(sigma=3., maxiters=10)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return bkg



def transmission(data,header,fg,bg,fgaper,bgaper,factor=2.,figname='',norm='linear',vmin='',vmax='',ZScale=False,celestial=False):
    """
    Make a transmission map with noise and s/n map. Optional to plot them with PlotSix routine.
    
    Plot in PlotSix 

    Parameters
    ----------
    data,header,fg,bg,fgaper,bgaper,factor=2.,figname=''

    data : numpy array
         Two-dimensional image data

    header : fits header instance
         The FITS header information to be used for WCS information for plot.

    fgaper, bgaper : elliptical aperture structure (x0,y0,a,b,pa) 
         parameters used to define an elliptical apertutre: central pixel coordinates (X0,Y0), 
         major axis (a), minor axis (b) and position angle (pa) on the sky for galaxy aperture.
         fgaper : the foreground galaxy aperture. 
         bgaper : the background galaxy aperture. 

    factor : float
         Optional growth factor for the elliptical apertures for foreground and background galaxy both. 

    figname : string (optional)
         The path and filename for the figure output f the histogram of pixel value pre and post
         sigma clipping.

    Returns
    -------
    trans : numpy array
         Transmission map based on data, foreground and background galaxy images: T= (Data - FG)/ BG
         Same dimension as data. 

    sn : numpy array
         The signal-to-noise map based on the data

    DTmap : numpy array
         The uncertainty in the transmission map derived from the difference maps between foreground model 
         and data and background model and data as well as a sky value determined outside both galaxies.
         Values are quadratically added. Difference between galaxy and model is scaled wi

    """
    bg = fix_zeros(bg)

    # transmission data
    trans = (data-fg)/bg

    # generate the x,y pixels for the data
    ymax, xmax  = np.shape(data)
    y, x = np.mgrid[0:ymax, 0:xmax]

    # unpack the ellipse information for foreground and background apertures
    (x0_fg,y0_fg,a_fg,b_fg,pa_fg) = fgaper
    (x0_bg,y0_bg,a_bg,b_bg,pa_bg) = bgaper
    
    # print some useful info
    print(f'Aperture Growth Factor:\t {factor}')
    print(f'Foreground\t major axis:\t {a_fg}\t minor axis:\t{b_fg}')
    print(f'Background\t major axis:\t {a_bg}\t minor axis:\t{b_bg}')

    # Foreground Galaxy Mask
    theta_fg = Angle(pa_fg, 'deg')
    ell_fg = Ellipse2D(amplitude=1., x_0=x0_fg, y_0=y0_fg, a=factor*a_fg, b=factor*b_fg, theta=theta_fg.radian)
    aperture_fg = ell_fg(x, y)
        
    # Background Galaxy Mask
    theta_bg = Angle(pa_fg, 'deg')
    ell_bg = Ellipse2D(amplitude=1., x_0=x0_bg, y_0=y0_bg, a=factor*a_bg, b=factor*b_bg, theta=theta_bg.radian)
    aperture_fg = ell_fg(x, y)

    # Get the background
    bkgnd = background(data)

    # the variance due to the sky 
    DeltaI = np.std(bkgnd.data)

    # compute the residual image of the foreground image
    residual_fg = data-fg
    filtered_residual_fg = sigma_clip(residual_fg[np.where(ell_fg(x, y)>0)], sigma=2., maxiters=None, cenfunc='mean', masked=False, copy=False)

    # compute the residual image of the background image
    residual_bg = data-bg
    filtered_residual_bg = sigma_clip(residual_bg[np.where(ell_bg(x, y)>0)], sigma=2., maxiters=None,cenfunc='mean', masked=False, copy=False)

    # compute the uncertainties from the foreground/background asymmetry
    DeltaB = np.std(np.ravel(filtered_residual_bg))
    DeltaF = np.std(np.ravel(filtered_residual_fg))

    # DTmap = np.sqrt( (DeltaI/bg)**2 + (DeltaF/bg)**2 + (DeltaB/bg**2)**2 )

    # Transmission uncertainty
    DTmap = np.sqrt( (residual_fg/np.max(np.ravel(residual_fg))*DeltaF )**2. + (residual_bg/np.max(np.ravel(residual_bg))*DeltaB)**2 + DeltaI**2 )
   
    # Signal to noise
    sn = trans/DTmap

    # plt.imshow(trans, vmin=0, vmax=1)

    if figname != '':
        PlotSix(data,fg,bg,trans,DTmap,sn,header,
              title1='Data', title2='Foreground Galaxy Model', title3='Background Galaxy Model', title4='Transmission Map', title5='Transmission Uncertainty', title6='S/N Map', 
              figname=figname, 
              vmin=vmin, vmax=vmax, 
              norm=norm, ZScale=ZScale, celestial=celestial)

    return trans,sn,DTmap

    
def attenuation(data,vmin=0,vmax=3,xlim='',ylim='',figname=''):
    if xlim=='' and ylim=='':
        xlim=(0,data.shape[0])
        ylim=(0,data.shape[1])

    data = fix_zeros(data)
    a = -1.086*np.log(data)
    plt.figure(dpi=300)
    plt.imshow(a,vmin=vmin,vmax=vmax,interpolation='none')
    plt.title('Attenuation Map')
    plt.colorbar(label='$A_V$')
    plt.xlim(xlim)
    plt.ylim(ylim)

    if figname != '':
        plt.savefig(figname)



def RadialTransmission (trans,DTmap,header,aper,z,figname=''):
    """
    Compute the radial transmission curve for the foreground galaxy in the aperture.

    Compute 

    Parameters
    ----------
    trans : numpy array
        The transmission map computed elsewhere. 

    DTmap : numpy array
        The uncertainty map of the transmission map. Same dimension as the transmission map.

    header : fits header instance
         The FITS header information to be used for WCS information for plot.

    aper : elliptical aperture structure (x0,y0,a,b,pa) 
         parameters used to define an elliptical apertutre: central pixel coordinates (X0,Y0), 
         major axis (a), minor axis (b) and position angle (pa) on the sky for galaxy aperture.

    z : float
         Redshift of the target galaxy. Used to compute the distance using the Planck 2018 cosmology.
         This distance is used to convert angular radii into physical ones expressed in kiloparsec. 
 
    figname : string (optional)
         The path and filename for the figure output of the radial plot. 

    Returns
    -------
    out : type
        Explanation of `out`.

    """
    
    
    (x0,y0,a,b,pa) = aper

    if ('CDELT1' in header) == True:
        pixscale1 = header['CDELT1']
    elif ('CD1_1' in header) == True:
        pixscale1 = header['CD1_1']
         
    if ('CDELT2' in header) == True:
        pixscale2 = header['CDELT2']
    elif ('CD2_2' in header) == True:
        pixscale2 = header['CD2_2']
    #
    ra0  = header['CRVAL1'] + pixscale1*( x0-header['CRPIX1'] )
    dec0 = header['CRVAL2'] + pixscale2*( y0-header['CRPIX2'] )

    ymax, xmax  = np.shape(trans)
    y, x = np.mgrid[0:ymax, 0:xmax]

    ra  = header['CRVAL1'] + pixscale1*( x-header['CRPIX1'] )
    dec = header['CRVAL2'] + pixscale2*( y-header['CRPIX2'] )

    #### for generating a sequence of mu and phi values
    (mu,phi) = proper_motion(ra0,dec0,ra,dec)

    ### Inclination
    bamin = 0.1
    cosinc = np.sqrt( ( (b/a)**2 - (bamin)**2 )/(1-(bamin)**2)) # definition from Hubble (1926)
    
    # Make a radius map (in pixels)

    #### for converting all PHI values to radians betweeon -pi and +pi
    phi[np.where(phi>np.pi)] -= 2*np.pi

    R = mu*np.sqrt( (np.sin(phi - pa)/cosinc)**2 + np.cos(phi - pa)**2)

    # sort R, trans and DTmap according to radius
    
    
    # return to pixel scale
    R_pix = R/pixscale1 # returning this to pixels

    from astropy.cosmology import Planck18 as cosmo
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z)/60. 
    
    R_kpc = R*kpc_per_arcsec
    
    incl = (180.0/np.pi)*np.arccos(cosinc)
    if np.isnan(incl) == True:
        incl = 90. # edge-on or poorly constrained inclination systems

    fig = plt.figure(figsize=(12, 12), dpi=150)
    plt.plot(np.ravel(R),np.ravel(trans))
    plt.xlabel('Radius (")')
    plt.ylabel('Transmission')
    fig.savefig(figname)
     
    

    return np.ravel(trans), np.ravel(R), np.ravel(R_pix), np.ravel(R_kpc)




def proper_motion(ra0,dec0,ra1,dec1):
    mu_dec = dec1-dec0                                   # difference in declinations
    mu_ra = ra1-ra0                                      # difference in right ascension
    mu = np.sqrt( (mu_ra*np.cos(dec1))**2 + mu_dec**2 )     # approximate mu for very nearby points
    alpha = mu_ra*np.cos(dec0)/mu                                             # sin(phi) = alpha
    beta = mu_dec/mu                                                            # cos(phi) = beta

    phi = np.arccos(beta)
    phi[np.where(alpha<0.0)] = 2*np.pi - np.arccos(beta[np.where(alpha<0.0)])

    return (mu,phi)



