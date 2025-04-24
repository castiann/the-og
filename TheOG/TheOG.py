# this is the main TheOG file
# import everything we need here;

# NOTE: Currently unused, due to import issues





import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams["savefig.dpi"] = 150
rcParams["figure.dpi"] = 150
rcParams["font.size"] = 10
rcParams['savefig.bbox'] = 'tight'

from astroquery.mast import Tesscut

# astropy
from astropy.table import Table
from astropy.io import fits
from astropy.stats import sigma_clip # function
from astropy.stats import SigmaClip # class

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord

from astropy.nddata.utils import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import rotate




from photutils.background import Background2D, SExtractorBackground
from photutils import Background2D, MedianBackground
from photutils.isophote import build_ellipse_model

from Image_tools import *
# import Image_tools
# from . import Image_tools
# from TheOG import Image_tools
# from .TheOG import Image_tools

def TheOG():
    print("TheOG has been fully loaded.")


    return 0


def og(filename, regfile):
    """
    This is the fully automated version of overlapping galaxy analysis.
    Once a fitsfile and a DS9 .reg file (in pixels) are available,  
    """

    
    
    
    
    return 0


def close_og():
    """
    The same as above but iterative as we model a close pair.
    """
    
