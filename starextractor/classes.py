import numpy as np
from photutils.aperture import CircularAperture,aperture_photometry
from photutils.psf import IntegratedGaussianPRF,PSFPhotometry,IterativePSFPhotometry
from photutils.background import Background2D
from scipy.spatial import KDTree
from astropy.utils import lazyproperty
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.table import vstack

from .image import read_image,source_extract,fwhm_estimate,_make_cutouts
from .plot import show_image
from .invariantfeatures import _generate_invariants

def parse_image(image_in):
    """
    Read and parse an astronomical image. Currently, image formats with .fits, generic image format(such as .bmp), .npy, and numpy array are supported.

    Usage:
        >>> from starextractor import parse_image
        >>> imagefile = 'obs/fits/img_00000.fits'
        >>> # imagefile = 'obs/bmp/img_00000.bmp'
        >>> # imagefile = 'obs/npy/img_00000.npy'
        >>> data = parse_image(imagefile)
        >>> # data = parse_image(image_array)
    Inputs:
        image_in -> [str or numpy array] image file or numpy array
    Outputs:
        data -> AstroImage object with attributes as follows:
            image_raw -> [2d array of float] Original grayscale image
            image -> [2d array of float] Image of background subtracted
            res -> [tuple of float] Resolution, such as [1024,1024]
            bkg -> [2d array of float] Grayscale background 
            bkg_rms -> [2d array of float] Background noise levels
            fwhm -> [float] Full width at half maximum (FWHM) of the gaussian kernel for the image.
            sigma -> [float] Sigma corresponding to FWHM. The relationship between the two is FWHM = 2sqrt{2ln2}sigma ≈ 2.355 sigma
            _offset -> [tuple of float] Pixel coordinates of the center of the image w.r.t. the center of the lower left pixel, such as [511.5,511.5]
            _bkg_sigma -> [float] Median of the background noise levels
            _if_rectangle -> [tuple] The lower left corner point, width and height of the rectangle constructed in the center part of the image, once the edge area of the image is masked.
            -if_region -> [2d array of bool] Masked array once the edge area of the image is masked.
    """ 
    # Read the astronomical image
    if type(image_in) is str:    
        image_raw = read_image(image_in)
    else:
        image_raw = image_in
                
    res = yc,xc = image_raw.shape

    # Calculate the offset of the image center from the origin
    offset = np.array([xc/2,yc/2]) - 0.5
            
    # Calculate the background with classic SExtractorBackground()
    bkg = Background2D(image_raw, (yc//16, xc//16))
    image = image_raw - bkg.background # image of background-subtracted
    bkg_rms = bkg.background_rms
    bkg_sigma = bkg.background_rms_median

    # A rectangle in the center part of the image is constructed, once the edge area of the image is masked.
    if_region = np.ones(res,dtype=bool)
    bb,lb = yc//4,xc//4
    ub,rb = bb*3,lb*3
    width,height = xc//2,yc//2
    if_rectangle = ([lb,bb],width,height) # The lower left corner point, width and height of the rectangle
    if_region[bb:ub, lb:rb] = False # Pixels in the edge areas of the image are masked
    info = {'image_raw':image_raw,'res':res,'_offset':offset,'image':image,'bkg':bkg.background,'bkg_rms':bkg_rms,'_bkg_sigma':bkg_sigma,'_if_rectangle':if_rectangle,'_if_region':if_region}

    return AstroImage(info) 

class AstroImage(object):
    """
    Class AstroImage
        Attributes:
            image_raw -> [2d array of float] Original grayscale image
            image -> [2d array of float] Image of background subtracted
            res -> [tuple of float] Resolution, such as [1024,1024]
            bkg -> [2d array of float] Grayscale background 
            bkg_rms -> [2d array of float] Background noise levels
            fwhm -> [float] Full width at half maximum (FWHM) of the gaussian kernel for the image
            sigma -> [float] Sigma corresponding to FWHM. The relationship between the two is FWHM = 2sqrt{2ln2}sigma ≈ 2.355 sigma
            _offset -> [tuple of float] Pixel coordinates of the center of the image w.r.t. the center of the lower left pixel, such as [511.5,511.5]
            _bkg_sigma -> [float] Median of the background noise level
            _if_rectangle -> [tuple] The lower left corner point, width and height of the rectangle constructed in the center part of the image, once the edge area of the image is masked.
            -if_region -> [2d array of bool] Masked array once the edge area of the image is masked.
        Methods:
            - show -> Show the original grayscale image.
            - find_source -> Search for star spots, and extract their centroids and do photometry.
    """    
    def __init__(self,info):  
        """
        Initialize an instance of class AstroImage.
        """
        for key in info.keys():
            setattr(self, key, info[key])    

    def __repr__(self):
        """
        Returns a more information-rich string representation of the AstroImage object.
        """
        return '<AstroImage object: RES = {} FWHM ≈ {:.2f}>'.format(self.res,self.fwhm)     

    @lazyproperty
    def fwhm(self):
        """
        Estimate the Full width at half maximum (FWHM) of the gaussian kernel for the image.
        """
        image,bkg_sigma,_if_region = self.image,self._bkg_sigma,self._if_region
        try:
            fwhm = fwhm_estimate(image,bkg_sigma,_if_region)
        except:
            raise Exception('The central area of the image requires at least an approximately circular Gaussian spot.')    
        return fwhm

    def find_source(self,threshold=5,max_control_points=100,fwhm=None,edgemask=False,phot='aperture'):
        """
        Search for star spots in an image, extracting their centroids and doing photometry.

        Usage: 
            >>> from starextractor import parse_image
            >>> imagefile = 'obs/fits/img_00000.fits' 
            >>> data = parse_image(imagefile)
            >>> data.show()
            >>> sources = data.find_source()
            >>> sources.show()
        Inputs:    
            threshold -> [float,optional,default=5.0] The detection threshold at the n-sigma noise level, above which sources are selected.
            max_control_points -> [int,optional,default=100] Maximum number of sources to extract.
            fwhm -> [float,optional,default=None] Full-width half-maximum (FWHM) of the Gaussian kernel in units of pixels. If None, FWHM is estimated by fitting a radial profile.
            edgemask -> [bool,optional,default=False] whether mask the edge area of the image. If True, the edge area is masked.
            phot -> [str,optional,default='aperture'] Method of photometry. Avaliable options include 'aperture', 'psf' and 'dao'. 
            The aperture photometry is widely used in general star field photometry, but less suitble for crowded star fields and when stars are located on the image boundaries.
            The Point Spread Function(PSF)/Pixel Response Function(PRF) photometry solved the problem encountered by the aperture photometry, and can estimate the brightness of stars more accurately.
            The DAOPHOT Photometry, essentially iterative PSF photometry, is useful for crowded fields where faint stars are very close to bright stars. 
            The faint stars may not be detected until after the bright stars are subtracted, so it can find more faint stars by a number of iterations.
        Outputs:
            sources -> Source object with attributes as follows:
                xy -> [2d array of float] Pixel coordinates of the star centroids w.r.t. the center of image
                brightness -> [array of float]  Brightness(sum of grayvalues) of star spots
                snr -> [array of float] Signal Noise Ratio(SNR) of star spots
                phot -> [str] Method of photometry
                edgemask -> [bool] whether the edge area of the image is masked.
                _offset -> [array of float] Pixel coordinates of the center of the image w.r.t. the center of the lower left pixel
                _image_raw -> [2d array of float] Original grayscale image
                _apertures -> CircularAperture objects for sources
                _if_rectangle -> [tuple] The rectangle constructed in the center part of the image, once the edge area of the image is masked.
        """
        if fwhm is None: 
            fwhm = self.fwhm
        else:
            self.fwhm = fwhm

        if edgemask:
            _if_region = self._if_region
        else:
            _if_region = None     

        # Convert FWHM to Sigma by FWHM = 2sqrt{2ln2}sigma ≈ 2.355 sigma
        self.sigma = sigma = fwhm * gaussian_fwhm_to_sigma

        image,bkg_sigma,bkg_rms,offset = self.image,self._bkg_sigma,self.bkg_rms,self._offset

        # Detect stars in an image using the DAOFIND algorithm.
        daofind = source_extract(bkg_sigma,fwhm,threshold,max_control_points)

        if phot == 'aperture':
            # Aperture photometry
            star_spots = daofind(image,mask=_if_region)
            xy_centroids = np.transpose((star_spots['xcentroid'], star_spots['ycentroid']))
            # Generate a circular aperture based on the pixel centroid coordinates and FWHM of the star spots
            apertures = CircularAperture(xy_centroids,r=fwhm) 
            # Calculate the starlight flux within the aperture
            phot_tbl = aperture_photometry(image,apertures) 
            brightness = phot_tbl['aperture_sum'].value
            # Calculate the background noise within the aperture
            noise_tbl = aperture_photometry(bkg_rms,apertures)
            noise = noise_tbl['aperture_sum'].value/np.sqrt(fwhm**2*np.pi)
        elif phot in ['psf','dao']:
            # Generate an integrated gaussian PRF model based on the sigma of the gaussian kernel
            psf_model = IntegratedGaussianPRF(sigma=sigma)
            fit_shape = daofind.kernel.shape # Size of of the gaussian kernel in pixels
            # Point Spread Function(PSF)/Pixel Response Function(PRF) photometry
            psfphot = PSFPhotometry(psf_model,fit_shape,finder=daofind,aperture_radius=fwhm)
            phot_tbl = psfphot(image,error=bkg_rms,mask=_if_region)
            # create the residual image, i.e. subtracting the fitted PSF models from the background-subtracted image
            image_resi = psfphot.make_residual_image(image,fit_shape) 
            
            if phot == 'dao':
                # Iterative PSF Photometry, which is an implementation of the DAOPHOT algorithm
                # Considering that the residual image may introduce fictitious sources, it is necessary to shrink the filtering threshold for source search.
                daofind = source_extract(bkg_sigma,fwhm,threshold,max_control_points,sharplo=0.5,sharphi=1,roundlo=-1,roundhi=1)
                psfphot = IterativePSFPhotometry(psf_model,fit_shape,daofind,aperture_radius=fwhm)
                phot_resi_tbl = psfphot(image_resi,error=bkg_rms,mask=_if_region)
                # combine tables of iterative photometry
                phot_tbl['iter_detected'] = 0
                phot_resi_tbl.meta = {} # prevent merge conflicts
                phot_tbl = vstack([phot_tbl, phot_resi_tbl])

            # Filter out invalid sources with flag values greater than 1:
            # 0 : normal source
            # 1 : one or more pixels in the ``fit_shape`` region were masked
            # 2 : the fit x and/or y position lies outside of the input data
            # 4 : the fit flux is less than or equal to zero
            # 8 : the fitter may not have converged
            # 16 : the fitter parameter covariance matrix was not returned
            valid_flags = phot_tbl['flags'] <= 1
            phot_tbl = phot_tbl[valid_flags]
            xypos = phot_tbl['x_fit','y_fit']
            xy_centroids = np.transpose((phot_tbl['x_fit'], phot_tbl['y_fit']))
            apertures = CircularAperture(xy_centroids,r=fwhm) # just for plot apertures in original image.   
            brightness = phot_tbl['flux_fit'].value # Calculate the fitted starlight flux within the kernel square
            noise = _make_cutouts(bkg_rms,xypos,fit_shape).sum(axis=(1,2))/fit_shape[0] # Calculate the background noise within the kernel square  
        else:
            raise Exception('Currently, feasible photometry methods include Circular Aperture photometry, Point Spread Function(PSF)/Pixel Response Function(PRF) photometry, and DAOPHOT(Iterative PSF) photometry.')      
        
        xy = xy_centroids - offset # move the origin to the center of image from the center of the lower left pixel.
        snr = brightness/noise # Signal Noise Ratio

        di = np.argsort(brightness)[::-1] # descending index
        dict_values = xy[di],apertures[di],brightness[di],snr[di],offset,edgemask,self.image_raw,self._if_rectangle,phot
        dict_keys = 'xy','_apertures','brightness','snr','_offset','edgemask','_image_raw','_if_rectangle','phot'

        info = dict(zip(dict_keys, dict_values))

        return Source(info)      
    
    def show(self,mode='image_raw',fig_out=None):
        """
        Show grayscale image of original/background/background-subtracted/background noise levels.

        Usage:
            >>> data.show('bkg')
        Inputs:
            mode -> [str,optional,default='image_raw'] Type of grayscale image to display. Available options include
                'image_raw' : Original grayscale image
                'image' : Image of background subtracted
                'bkg' : Grayscale background 
                'bkg_rms' : Background noise levels
            fig_out -> [str,optional,default=None] Path of the output figure.
        Outputs:
            figures    
        """
        if mode == 'image_raw':
            img = self.image_raw
        elif mode == 'bkg':
            img = self.bkg  
        elif mode == 'image':
            img = self.image
        elif mode == 'bkg_rms':
            img = self.bkg_rms   
        else:
            raise Exception("Type of grayscale image is not supported. Available options include 'image_raw', 'image', 'bkg', and 'bkg_rms'.")           

        if fig_out is None:
            show_image(img,origin='lower') 
        else:    
            plot_kwargs = {'figname':fig_out}
            show_image(img,origin='lower',**plot_kwargs)          

class Source(object):
    """
    Class Source
    """
    
    def __init__(self,info): 
        """
        Initialize an instance of class Source.
        """ 
        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
        """
        Returns a more information-rich string representation of the Source object.
        """
        return "<Source object: NUM = {:d} OFFSET = {} EDGEMASK = {} PHOT = '{:s}'>".format(len(self.xy),self._offset,self.edgemask,self.phot.upper())

    def show(self,fig_out=None):
        """
        Show original image with sources marked.

        Usage:
            >>> sources.show()
        Inputs:
            fig_out -> [str,optional,default=None] Path of the output figure.
        Outputs:
           figures    
        """
        if fig_out is None:
            plot_kwargs = {'mark':(self._apertures,'blue')}
        else:
            plot_kwargs = {'mark':(self._apertures,'blue'),'figname':fig_out}
        if self.edgemask:
            mask_rectangle = self._if_rectangle
        else:
            mask_rectangle = None    
                
        show_image(self._image_raw,origin='lower',mask_rectangle=mask_rectangle,**plot_kwargs)

    def invariantfeatures(self,max_control_points=None):
        """
        1. Calculate the unique invariants (L2/L1,L1/L0), where L2 >= L1 >= L0 are the three sides of the triangle composed of centroids.
        2. Construct the 2D Tree from the the unique invariants.
        3. Record an array of the indices of centroids that correspond to each invariant.

        Usage:
            >>> xy,asterisms,kdtree = sources.invariantfeatures()
        Inputs:
            max_control_points -> [int,optional,default=None] Maximum number of sources used to compute invariant features
        Outputs:
            xy -> [2D array(nx2) of float] Unique invariants (L2/L1,L1/L0)
            asterisms -> [2D array(nx23) of int] Indices of centroids that build the invariant triangle
            kdtree -> [scipy.spatial.KDTree] 2D Tree from the the unique invariants for quick nearest-neighbor lookup.
        """
        n = len(self.xy)
        if max_control_points is None or max_control_points > n: max_control_points = n
        xy = self.xy[:max_control_points]

        inv_uniq,asterisms = _generate_invariants(xy)
        if not inv_uniq.size:
            raise Exception('At least three sources are required to form an invariant triangle.')
        kdtree = KDTree(inv_uniq)

        return xy,asterisms,kdtree  