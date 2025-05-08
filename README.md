# [the-og](https://github.com/castiann/the-og)
Created by [Ian Castellano](https://orcid.org/0009-0006-0104-6045) and [Benne Holwerda](https://orcid.org/0000-0002-4884-6756).

This program is used to map dust attenuation of galaxies using the overlapping galaxy technique. The methods used are the rotation method and the isophotal method, as demonstrated in [Benne Holwerda 2009](https://ui.adsabs.harvard.edu/abs/2009AJ....137.3000H/abstract).

## Using the Program
This project was created using Python 3.12.3.

To start the program, I recommend first looking at the example notebook provided using the rotation method on IC 720 (`ic_720_rotate.ipynb`) and following along with the [documentation file](./documentation). These files guide you through the steps of opening and analyzing the overlapping galaxy pair. The files can then be duplicated to use as a template. Note: as of the current version of this project, the isophotal method on IC 720 (`ic720_isophote.ipynb`) is currently still a work in progress, and will be updated with a newer release.

In the first python cell, you will specify the path to the data location and the export location, and give the data set nickname.

## File Structure
- Project notebooks must be in the main directory the-og. For example, `the-og/ic720_rotate.ipynb`.
- Data should be kept in subdirectories. For example, `the-og/ic720_data/ic720.fits`.
- Results and exports will be sent to subdirectories. For example, `the-og/ic720_output`. Create this directory before running the code.
