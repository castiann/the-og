# [the-og](https://github.com/castiann/the-og)
Created by [Ian Castellano](https://orcid.org/my-orcid?orcid=0009-0006-0104-6045) and [Benne Holwerda](https://orcid.org/0000-0002-4884-6756).

This program is used to map dust attenuation of galaxies using the overlapping galaxy technique. The methods used are the rotation method and the isophotal method, as demonstrated by [Benne Holwerda 2009](https://ui.adsabs.harvard.edu/abs/2009AJ....137.3000H/abstract). This project is shared under an attribution, noncommerical, share-alike license ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)).

## Using the Program
This project was created using Python 3.12.3.

To start the program, I recommend first looking at the example notebook provided using the rotation method on IC 720 (ic_720_rotate.ipynb). The file guides you through the steps of opening and analyzing the overlapping galaxy pair. The file can then be duplicated to use as a template. Note: as of the current version of this project, the isophotal method on IC 720 (ic_720_isophote.ipynb) is currently still a work in progress, and will be updated with a newer release.

In the first python cell, you will specify the path to the data location and the export location, and give the data set nickname.

## File Structure
- Project notebooks must be in the main directory the-og. For example, the-og/ic_720_rotate.ipynb.
- Data should be kept in subdirectories. For example, the-og/mw_compare_data/ic720.fits.
- Results and exports will be sent to subdirectories. For example, the-og/ic_720_output/ic_720_transmission_single.png. Create this directory before running the code.
