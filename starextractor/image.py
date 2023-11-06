import numpy as np
from photutils.detection import DAOStarFinder
from photutils.profiles import RadialProfile
from photutils.centroids.gaussian import centroid_1dg
from astropy.io import fits
from astropy.nddata import extract_array
from PIL import Image

def _make_cutouts(data,xypos,cutout_shape):
    """
    Extract smaller arrays of the given shape and positions from a larger array.
    """
    cutouts_data = []
    for xpos, ypos in xypos:
        cutouts_data.append(extract_array(data,cutout_shape, (ypos, xpos),fill_value=0))
    return np.array(cutouts_data)

def read_image(imagefile):
    """
    Read an astronomical image file. 
    Currently, supported image formats include .fits, generic image format(such as .bmp), and .npy.

    Usage:
        >>> imagefile = 'obs/fits/img_00000.fits' 
        >>> # imagefile = 'obs/bmp/img_00000.bmp'
        >>> # imagefile = 'obs/npy/img_00000.npy'
        >>> image_raw = read_image(imagefile)
    Inputs:
        imagefile -> [str] path of image file
    Outputs:
        image_raw -> [2d array of float] Original grayscale image
    """
    if imagefile.split('.')[1] in ['fits','fit']: # Load an astronomical image in format of fits
        unit_list = fits.open(imagefile)
        image_raw = unit_list[0].data.astype(float) 
        # convert origin from top left to bottom left
        image_raw = image_raw[::-1] 
    elif imagefile.split('.')[1] == 'npy':    
        image_raw = np.load(imagefile)
    else:
        image_raw = Image.open(imagefile).convert('L') # Load an astronomical image in generic image format
        image_raw = np.asarray(image_raw)
        # convert origin from top left to bottom left
        image_raw = image_raw[::-1]

    return image_raw  

def fwhm_estimate(image,bkg_sigma,mask_region,fwhm=10,sigma_radius=2.5,max_control_points=200):
    """
    Estimate the Full width at half maximum (FWHM) of the gaussian kernel for a given image.

    Usage:
        >>> fwhm = fwhm_estimate(image,bkg_sigma,mask_region)
    Inputs:
        image -> [2d array of float] Image of background subtracted
        bkg_sigma -> [float] Median of the background noise level
        mask_region -> [2d array of bool] Masked array for the case that edge area of the image is masked.
        fwhm -> [float,optional,default=10.0] Full width at half maximum (FWHM) of the gaussian kernel for the image.
        sigma_radius -> [float,optional,default=2.5] Radius of the gaussian kernel in unit of sigma.
        max_control_points -> [int,optional,default=200] Maximum number of sources to extract.
    """
    # Detect stars in an image and extract their centroids using the DAOFIND algorithm.
    daofind = DAOStarFinder(5*bkg_sigma,fwhm,sigma_radius=sigma_radius,brightest=max_control_points,sharplo=0.2,sharphi=1,roundlo=-0.05,roundhi=0.05) 
    star_spots = daofind(image,mask=mask_region)  
    xy_centroids = np.transpose((star_spots['xcentroid'], star_spots['ycentroid']))
    xypos = star_spots['xcentroid','ycentroid']
    cutout_shape = daofind.kernel.shape
    cutouts_image = _make_cutouts(image,xypos,cutout_shape)
    # Calculate the radial profiles using concentric circular apertures.
    edge_radii = np.arange(cutout_shape[0])
    rps = []
    for i in range(len(xypos)):
        xycen = centroid_1dg(cutouts_image[i])
        rp = RadialProfile(cutouts_image[i], xycen, edge_radii)
        rps.append(rp.gaussian_fwhm)
    fwhm = np.mean(rps)
    return fwhm    

def source_extract(bkg_sigma,fwhm,threshold,max_control_points,sharplo=-3,sharphi=3,roundlo=-2,roundhi=2):
    """
    Detect stars in an image and extract their centroids using the DAOFIND algorithm.

    Usage: 
        >>> daofind = source_extract(bkg_sigma,fwhm,threshold,max_control_points)
    Inputs:
        bkg_sigma -> [float] Median of the background noise level
        fwhm -> [float] Full width at half maximum (FWHM) of the gaussian kernel for the image.
        threshold -> [float] The detection threshold at the n-sigma noise level, above which sources are selected.
        max_control_points -> [int] Maximum number of sources to extract.
        sharplo -> [float,optional,default=-3] The lower bound on sharpness for object detection.
        sharphi -> [float,optional,default=3] The upper bound on sharpness for object detection.
        roundlo -> [float,optional,default=-2] The lower bound on roundness for object detection.
        roundhi -> [float,optional,default=2] The upper bound on roundness for object detection.
    Outputs:
        daofind -> DAOStarFinder object
    Note: 
        1. Roundness and sharpness measure the deviation of the energy distribution of the star spot from that of the Gaussian function.
        2. The roundness measures the bilateral (2-fold) to four-fold symmetry of the energy distribution of the star spot. 
        As the value approaches zero, the energy distribution approximates a pattern of decreasing concentric rings.
        3. The sharpness measures the radial profile of the energy distribution of the star spot. 
        The closer the value is to one, the closer the energy distribution is to a one-dimensional Gaussian function.   
    """  
    # DAOStarFinder searches convolved background-subtracted data image for local density maxima that have a peak amplitude greater than threshold.
    # DAOStarFinder finds the centroids of sources by fitting the marginal x and y 1D distributions of the Gaussian kernel to the marginal x and y distributions of the unconvolved background-subtracted data image.
    daofind = DAOStarFinder(threshold*bkg_sigma,fwhm,sigma_radius=2.5,brightest=max_control_points,sharplo=sharplo,sharphi=sharphi,roundlo=roundlo,roundhi=roundhi) 
    # the sigma_radius determines the size of a gaussian kernel, and the default value is approximately equal to FWHM;

    return daofind