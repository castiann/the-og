"""
This module contains all functions needed to analyze overlapping galaxies using the isophote technique.
The module utils contains most other imports, and access to utility functions.
"""
from . utils import *

def fit_isophote(data, aper,  # removed header from
                #  mask=np.zeros([2,2]),
                 use_mask=False, 
                 firstgood=True, isotablefile=''):

    """
    Fit an isophote within the elliptical aperture provided. 
    Use the mask to include/exclude pixels in that aperture.
    Firstgood means that the routine will go with the first modified aperture 
    that provides a meaningful isophotal fit. 

    Follows this tutorial: https://photutils.readthedocs.io/en/latest/user_guide/isophote.html

    Parameters:
        data : numpy array
            Two-dimensional image data

        header : fits header
            The header information to go with the data.

        aperture : elliptical aperture
            aperture used in the isophote fit. The aperture is sized by a scale factor (alternating up and down) until a fit is achieved 

        mask: numpy array
            numpy array with 0 and 1 for pixels to be used or not. 
            For masking unrelated sources and other pixels not to be used in an isophote fit.

        firstgood : boolean
            Whether to go with the first acceptable isophote fit

        isotablefile : string (optional)
            The path and filename for the output of the isophote table. 

    Returns
        isolist : photutils.isophote.isophote.IsophoteList

        geometry : photutils.isophote.geometry.EllipseGeometry 
    
        aper : aperture for ellipse
            The aperture scaled to which one was used for the isofit.
    
    """
    # unpack aperture
    (x0,y0,a,b,pa) = aper

    # If there is a mask provided, make a masked version of the data.
    # if np.all(mask==0): # if using the default array
    #     print('No mask applied')
    # else:
    #     data = ma.masked_equal(data, mask) # note, this doesn't seem to work
    #     print(f'Mask applied: {data.shape}, {type(data)}')


    # fix any potential missing data
    # "A masked array is an array that has an associated mask to indicate which elements are valid and which are not. 
    # The concept of masking comes in handy when you want to ignore missing or invalid entries during computations."
    
    # clean up data by removing bad values
    # invalids = np.isfinite(data)
    # invalids = invalids * 1       # sets True to be 1, and False to be 0
    # invalids = np.abs(invalids-1) # inverts values
    # print(np.max(invalids))
    # print(np.mean(invalids))
    
    # data = ma.masked_invalid(data)
    # data = ma.masked_where(data < 0, data)
    
    if use_mask == True:
        mask = make_ellipse_mask(data, aper)
        mask = np.abs(mask-1)
        print(mask)
        print(f'Mask: {mask.shape}, {type(mask)}')
        data = ma.masked_equal(data, mask) # I don't think this is what we want
        # data = ma.array(data, mask=invalids)
        data = ma.array(data, mask=mask)
        # data = apply_mask(data, aper)
        print(f'Masked data; {data.shape}, {type(data)}')
        # return data
    else:
        print('No mask applied')



    if firstgood==False:
        diff_min = np.sum(np.abs(data)) # there is no model, max value of residual
        best_sizefactor = 1.
        
    for sizefactor in [1.0, 1.1,
                        # 0.9, 1.2, 0.8, 1.3, 0.7, 1.4, 0.6, 1.5, 0.5, 1.6, 0.4, 1.7, 0.3, 1.8, 0.2, 1.9, 0.1, 2.
                        ]: 
        print('Using size factor %4.1f on the aperture' % sizefactor)
        # define the geometry
        geometry = EllipseGeometry(x0=x0, y0=y0, sma=a*sizefactor, eps=(1. - b/a),
                            pa= (pa*np.pi/180.) )
        # the ellipse aperture first guess
        ellipse = Ellipse(data, geometry)
        print(type(ellipse))
        # quick_plot(ellipse)

        
        # fit isophotes 
        isolist = ellipse.fit_image()
        if firstgood==True and len(isolist.pa) > 0: # ok there is a solution
            # right now, this goes with the first solution, is there a way to identify the *best* solution or "good enough" solution?
            print('ok a size size factor of %4.1f on the aperture seems to work.' % sizefactor)

            # the aperture (also to be returned for later plotting)
            aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                                          geometry.sma * (1 - geometry.eps), geometry.pa)

            
            # save the isophotal model to table
            if isotablefile !='':
                iso_model_table = isolist.to_table()
                iso_model_table.write(isotablefile+'.csv', format='csv', overwrite=True)
                print('Isophote fit saved to %s.' % isotablefile)
            else:
                print('Isophote fit not saved to a file.')
            # return isolist, geometry, aper

        elif firstgood==False and len(isolist.pa) > 0:
            model_image = build_ellipse_model(data.shape, isolist)
            diff = np.sum(np.abs(data-model_image)) # I may need to apply the masking thing, same as jiggle and rotate
            if diff < diff_min:
                diff_min = diff
                best_sizefactor = sizefactor

    if firstgood==False:
        # define the geometry
        geometry = EllipseGeometry(x0=x0, y0=y0, sma=a*best_sizefactor, eps=(1. - b/a),
                            pa= (pa*np.pi/180.) )
        # the ellipse aperture first guess
        ellipse = Ellipse(data, geometry)
        aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                                          geometry.sma * (1 - geometry.eps), geometry.pa)

        
    return isolist, geometry, aper, data


def plot_isophote(data,isolist,aperture,figfile): # note, removed header
    """
    Plot the isophotal model

    plot a three-panel figure with the orginal data, the ellipse model, and the residual (data-model) shown.


    Parameters
    ----------
    data : numpy array
         Two-dimensional image data

    header : fits header
         The header information to go with the data.

    isolist : isophote list from fit_isophote
         the list of isophote values (center x- and y-position, major axis, position angle and ellipticity)

    aperture : 
         aperture used in the isophote fit. Over-plot over the residual image. 

    figname : string (optional)
         The path and filename for the figure output


    Returns
    -------
    none
    
    """
    model_image = build_ellipse_model(data.shape, isolist)
    residual = data - model_image

    # plot the data and aperture op top
    # plt.imshow(data, origin='lower')
    # aperture.plot(color='blue')

    imin, imax = apply_zscale(data)

    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows=1, ncols=3)
    fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)
    ax1.imshow(data, origin='lower', vmin=imin, vmax=imax)
    ax1.set_title('Data')

    imin, imax = apply_zscale(model_image)

    smas = np.linspace(10, 50, 5)
    for sma in smas:
        iso = isolist.get_closest(sma)
        x, y, = iso.sampled_coordinates()
        ax1.plot(x, y, color='white')
    aperture.plot(color='white')

    ax2.imshow(model_image, origin='lower', vmin=imin, vmax=imax)
    ax2.set_title('Ellipse Model')

    imin, imax = apply_zscale(residual)

    ax3.imshow(residual, origin='lower', vmin=imin, vmax=imax)
    ax3.set_title('Residual')
    plt.show()
    plt.savefig(figfile)

    return model_image
    

def plot_isofit(isolist,figname):
    """
    Plot in four panels the parameters of the isophotal fit.
    Adapted from the example here: https://photutils.readthedocs.io/en/stable/isophote-3.py (Link 404) 
    https://photutils.readthedocs.io/en/latest/user_guide/isophote.html ?
    """
    
    plt.figure(figsize=(8, 8))
    plt.subplots_adjust(hspace=0.35, wspace=0.35)

    plt.subplot(2, 2, 1)
    plt.errorbar(isolist.sma, isolist.eps, yerr=isolist.ellip_err,
             fmt='o', markersize=4)
    plt.xlabel('Semimajor Axis Length (pix)')
    plt.ylabel('Ellipticity')

    plt.subplot(2, 2, 2)
    plt.errorbar(isolist.sma, isolist.pa / np.pi * 180.,
                yerr=isolist.pa_err / np.pi * 80., fmt='o', markersize=4)
    plt.xlabel('Semimajor Axis Length (pix)')
    plt.ylabel('PA (deg)')

    plt.subplot(2, 2, 3)
    plt.errorbar(isolist.sma, isolist.x0, yerr=isolist.x0_err, fmt='o',
                markersize=4)
    plt.xlabel('Semimajor Axis Length (pix)')
    plt.ylabel('x0')

    plt.subplot(2, 2, 4)
    plt.errorbar(isolist.sma, isolist.y0, yerr=isolist.y0_err, fmt='o',
                     markersize=4)
    plt.xlabel('Semimajor Axis Length (pix)')
    plt.ylabel('y0')
    plt.savefig(figname)
    
def plot_isophote_fit(data,isolist,figname):
    """

    Parameters
    ----------
    data : numpy array
         Two-dimensional image data

    isolist : isophote list from fit_isophote
         The list of isophote values (center x- and y-position, major axis, position angle and ellipticity)

    figname : string (optional)
         The path and filename for the figure output


    Returns
    -------
    model_image : numpy array
         The elliptical isophotes as values in a numpy array the same size as data.

    residual : numpy array
         The resiudal of the data after subtraction of the isophotal model.

    """

    
    model_image = build_ellipse_model(data.shape, isolist)

    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows=1, ncols=3)
    fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)
    ax1.imshow(data, origin='lower')
    ax1.set_title('Galaxy Pair')

    smas = np.linspace(10, 50, 5)
    for sma in smas:
        iso = isolist.get_closest(sma)
        x, y, = iso.sampled_coordinates()
        ax1.plot(x, y, color='white')

    mod = ax2.imshow(model_image, origin='lower')

    ax2.set_title('Ellipse Model')
    ax3.imshow(residual_fg, origin='lower')
    ax3.set_title('Residual')
    fig.savefig(figdir+figname)

    return model_image, residual


def iterfit_isophote(data,header, x0, y0, a, b, pa, mask='', sigmalim=3, niter = 100, breaklim=0.1, isotablefile=''):
    """
    Fit an isophote within the ellipse provided. Fit until  
    Firstgood means that the routine will go with the first modified aperture 
    that provides a meaningful isophotal fit. 

    """
  
    # first iteration of the mask is the whole target object ellipse
    mask = mkmask(data,x0,y0,a,b,pa)
    ref_std = np.std(np.ravel( data[np.where(mask>0.)] )) # reference standard deviation (to start) is that of the data in the ellipse
    resid_std = 10000000.0

    for i in range(niter):
        print(i,ref_std, resid_std)
        # define the geometry
        geometry = EllipseGeometry(x0=x0, y0=y0, sma=a, eps=(1. - b/a),
                            pa= (pa*np.pi/180.) )
        # the ellipse aperture first guess
        ellipse = Ellipse(data*mask, geometry)

        # fit isophotes 
        # isolist = ellipse.fit_image()
        # print(len(isolist))
        try:
            isolist, geometry, aper = fit_isophote(data,header, x0, y0, a, b, pa, mask=mask, firstgood=False, isotablefile='')
        except:
            print('IsoFit did not work.')
        # the aperture (also to be returned for later plotting)
        # aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
        #                            geometry.sma * (1 - geometry.eps), geometry.pa)

        if len(isolist) > 0:
            # make a model image
            model_image = build_ellipse_model(np.shape(data), isolist)

            # make a mask based on the sigma clipping of the residual.
            mask = sigmaclip_mask(data-model_image, x0, y0, a, b, pa, sigmalim=sigmalim)

            # compute the residual of the ellipse fit in the aperture
            resid_std = np.std(np.ravel( (data-model_image)[np.where(mask>0.)]))
            # compare to the reference standard deviation
            if abs(resid_std - ref_std)  < breaklim*ref_std: # yes the std have converged to within breaklim
                print('Final iteration ',i,ref_std, resid_std)
                # save the isophotal model to table
                if isotablefile !='':
                    iso_model_table = isolist.to_table()
                    iso_model_table.write(isotablefile+'.csv', format='csv', overwrite=True)
                    print('Isophote fit saved to %s.' % isotablefile)
                else:
                    print('Isophote fit not saved to a file.')
                return isolist, geometry, aper
        else:
            ref_std = resid_std # the new limit to beat is that of the current residual
            print('Iteration %d did not converge to within %f of standard deviation, currently: %f' % (i,breaklim,resid_std))

    print('No solution was found.')
    return 0 # no solution was found
    

# sigma clip a pixel collection
# I *know* there is something for that in the photometry package.
