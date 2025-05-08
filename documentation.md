# Documentation for the-og
All additional documentation regarding the-og is located here.

## Set Up
For these instructions, the overlapping galaxy pair IC 720 is used as an example. When using the code for another overlapping pair, simply replace `ic720` with the desired name.

These steps will lead you through several of the convenience functions created for the set up of this program.
1. Import the directory containing data (in the form of .fits or .fit files) into the root the-og folder called `ic720_data`.
2. Create an empty directory in the-og called `ic720_output`. This is the directory that will contain a copy of the data (or a cutout of the data) and all files that are created with the code.
3. Create a Jupyter Notebook in the-og and call it something like `ic720_rotate`, `ic720_isophote`, or `ic720_analysis`.
4. In the first cell, paste the following code to import the needed modules:
    ```python
    from TheOG.utils import *
    import TheOG.rotate_tools as rot
    import TheOG.isophote_tools as iso
    %load_ext autoreload
    %autoreload 2
    ```
5. In the next cell, we specify the needed paths.
    - `datapath` is the relative path from the notebook to the .fits file.
    - `nick` is the nickname given to the file; this is used in the titles of plots.
    - `file` is the relative path from the notebook to the output directory *plus* the nickname, which will be included in the file names.
    - It's recommended to specify which index contains the image data; the default value of `sci` is 0, but it is sometimes 1.
    ```python
    datapath = 'ic720_data/ic720.fits'
    nick = 'ic720'
    file = 'ic720_output/ic720'
    sci = 0
    ```
6. Using the function `quick_open()` from the utils module shows a peak of the .fits header and plots a few images of the data. It has optional returns of the `hdulist0`, `image_data0`, and `prihdr0` (the 0's indicate these are for the original, unaltered file).
    ```python
    hdulist0, image_data0, prihdr0 = quick_open(datapath, sci)
    ```
7. At this step, choose whether or not to make a smaller cutout of the data. If a smaller cutout is needed, specify its center position and size. Then, you can quickly plot it using the `quick_plot` convenience function:
    ```python
    pos,size = (2210,3398),(2000, 2000)
    hdulist, image_data, prihdr = create_working_file(path=datapath, sci=sci, file=file, cutout=True, pos=pos, size=size)
    quick_plot(image_data, name='IC 720 Cutout')
    ```
    If no cutout is needed, then simply create the working file. By default `cutout=False`:
    ```python
    hdulist, image_data, prihdr = create_working_file(path=datapath, sci=sci, file=file)
    ```

## Apertures
Apertures are drawn around each of the galaxies.
1. Run the following code to draw apertures in DS9:
    ```python
    rot.drawregions(fitsfile=file+'.fits', science=sci)
    ```
2. Open the working file located at `ic720_output/ic720.fits`. Open the apertures by `region>open`, navigate to `ic720_output/ic720.reg`, and open the file.
3. Hit ctrl+r to enter region editing mode. Move the blue + to the center of the foreground galaxy, and the red x to the center of the background galaxy. Then adjust the blue ellipse to encompass the entire foreground galaxy, and the red ellipse to encompass the entire background galaxy. Here are some tips for this step:
    - Ensure that you create the ellipses on the working file, and not the original data file.
    - Use a scale where it's easy to identify the centers of the galaxy to place the center markers.
    - The closer the markers are to the center, the better the final transmission plots will be.
    - Use the information panels to set the ellipses to the same position as the center markers.
    - Use different scales to ensure that the ellipses encompass the entire galaxy, and most importantly all of the dust lanes.
    - Hold shift to rotate ellipses.
    - The ellipse regions should be overlapping.
    - Save the region with `region>save`, and specify the IMAGE coordinate system.
    - Specify the file extension .reg to avoid errors when saving.
Name this file `ic720_regions.reg`. This ensures that rerunning the code does not overwrite this file.
4. Read the regions into the notebook using `readregion`:
    ```python
    x0_fg,y0_fg,a_fg,b_fg,pa_fg, x0_bg,y0_bg,a_bg,b_bg,pa_bg  = rot.readregion(file+'_regions.reg')
    ```
5. Define the apertures, and plot the ellipses. Some images look better with `ZScale=True`, but if not it can be set to false. There are also parameters vmin, vmax, and norm to adjust the image.
    ```python
    aper_fg = (x0_fg,y0_fg,a_fg,b_fg,pa_fg)
    aper_bg = (x0_bg,y0_bg,a_bg,b_bg,pa_bg)

    rot.plot_ellipses(image_data, prihdr, aper_fg, aper_bg, figname=file+'_elip_aper.png', ZScale=True)
    ```

## Rotation Method
The rotation method:
1. The function `imrotate` will rotate the image 180 degrees (by default), and find the difference image (image data - rotated data).
    ```python
    # foreground
    fg_rot_model = rot.imrotate(data=image_data,header=prihdr, x0=x0_fg, y0=y0_fg, dimen=u.Quantity((5, 5), u.arcmin), rotangle=180., figname=file+'_foreground_rotated.png',ZScale=True)

    # background
    bg_rot_model = rot.imrotate(data=image_data,header=prihdr, x0=x0_bg, y0=y0_bg, dimen=u.Quantity((5, 5), u.arcmin), rotangle=180., figname=file+'_background_rotated.png',ZScale=True)
    ```
    The goal here is to *minimize* the difference i.e. make the target galaxy disappear as much as possible.
2. To optimize the difference image, we can try using different center points of the galaxy using `jiggle_and_rotate`. This will search the area around the center, calculate the difference image at this point, and choose the center that yields the smallest difference. This function can take anywhere from 2 to 10 minutes, so choose the area to search carefully. Smaller areas will let the code run faster.

    To aid in this step, there are two different modes: manual and automatic.
    - Manual:
        - Look at your image and choose a distance away from your center point where you think the best center is located. (`pixrange`).
        - Then choose how many steps within that range you want to check for a better center (`pixelstep`).
        ```python
        # foreground, manual rough (no fine tuning)
        fg_dx, fg_dy = rot.jiggle_and_rotate(data=image_data, header=prihdr, x0=x0_fg, y0=y0_fg, pixrange=20, pixelstep=4, figname=file+'_fg_jiggle.png', foreground=True)
        
        print(fg_dx,fg_dy)
        ```
    - Automatic:
        - Open DS9 and open the `ic720_regions.reg` file created earlier. Adjust ellipses to surround the area where you believe the best center is located.
        - Save the file as `ic720_center_regions.reg`.
        - Do not change the centers of the regions, as the code won't recognize this change.
        - Adding the .reg file as a parameter will enable automatic mode (`pixrange` and `pixelstep`, if entered, will be ignored).
        ```python
        # foreground, automatic rough (no fine tuning)
        fg_dx, fg_dy = rot.jiggle_and_rotate(data=image_data, header=prihdr, x0=x0_fg, y0=y0_fg, figname=file+'_fg_jiggle.png', regfile=file+'_center_regions.reg', foreground=True)
        
        print(fg_dx,fg_dy)
        ```
    The code will return the best *change* in x and y position to get to the best center. The centers chosen by `jiggle_and_rotate` can then be defined:
    ```python
    x0_fg_jiggle = x0_fg + fg_dx
    y0_fg_jiggle = y0_fg + fg_dy

    x0_bg_jiggle = x0_bg + bg_dx
    y0_bg_jiggle = y0_bg + bg_dy

    # If jiggle_and_rotate gave desirable results, add the change here:
    aper_fg = (x0_fg_jiggle,y0_fg_jiggle,a_fg,b_fg,pa_fg)
    aper_bg = (x0_bg,y0_bg,a_bg,b_bg,pa_bg)
    ```
    From my experience, this function *does* improve the center point for the foreground galaxy, and in turn gives better results for the transmission plot. However, it tends to give a worse center for the background galaxy. Use the new center values only if the difference image seems to have improved from the original.

3. Use the `transmission` function to calculate and plot the transmission, as well as uncertainty and S/N.
    ```python
    trans,sn_trans,DTmap = rot.transmission(image_data, prihdr, fg_rot_model.data, bg_rot_model.data, aper_fg, aper_bg,figname=file+'_transmission', ZScale=True, celestial=False, norm='linear')
    ```

## Isophote Method

