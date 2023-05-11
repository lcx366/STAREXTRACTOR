from photutils.aperture import CircularAperture,aperture_photometry

from .image import read_image,source_extract
from .plot import show_image

class AstroImage(object):
    """
    Class AstroImage
    """    

    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'instance of class AstroImage'

    def read_image(imagefile):
        """
        Read an astronomical image captured by cameras in fits format or in generic image format.

        Usage:
        >>> from sourceextractor import AstroImage
        >>> imagefile = 'obs/fits/img_00000.fits' # imagefile = 'obs/bmp/img_00000.bmp'
        >>> image = AstroImage.read_image(imagefile)

        Inputs:
        imagefile -> [str] image filename

        Outputs:
        image -> instance of class AstroImage with its attributes as follows:
            rawimage -> [2d array of float] raw grayscale image
            res -> [array of float] Camera resolution
        """    
        image_raw = read_image(imagefile)
        res = image_raw.shape
        info = {'image_raw':image_raw,'res':res}

        return AstroImage(info)  

    def find_source(self,max_control_points=60,fwhm=12,mask=False):
        """
        Search for star spots in an image, extracting their centroids and doing photometry.

        Usage: 
        >>> from sourceextractor import AstroImage
        >>> imagefile = 'obs/fits/img_00000.fits' # imagefile = 'obs/bmp/img_00000.bmp'
        >>> image = AstroImage.read_image(imagefile)
        >>> sources = image.find_source(imagefile,mask=True)

        Parameters:    
        max_control_points -> [int,optional,default=50] Maximum number of sources to extract
        fwhm -> [float,optional,default=15] Priori Full-width half-maximum (FWHM) of the Gaussian kernel in units of pixels used in DAOStarFinder
        mask -> [bool,optional,default=False] A True value indicates the edge area of an image is masked. Masked pixels are ignored when searching for stars.
        Outputs:
        sources -> Instance of class Source with its attributes as follows:
            xy_centroids -> [2d array of float] Pixel coordinates of the star centroids
            _offset -> [array of float] Pixel coordinates of the center of the image
            _image_raw -> [2d array of float] raw grayscale image
            _image -> [2d array of float] signal, i.e. subtracting the background gray value from the raw grayscale image
            _apertures -> A circular aperture defined in pixel coordinates.
            brightness -> [array of float]  Brightness(Grayvalues) of star spots
            snr -> [array of float] Signal Noise Ratio(SNR) of star spots
            mask_rectangle -> [None or tuple] If None, then no mask rectangle is generated; Else, a rectangle defined by the bottom left corner point, width and height is generated
        """
        image_raw = self.image_raw
        # Search for star spots in an image and extract them
        xy_centroids,_offset,_image,_bkg_rms,_mask_rectangle = source_extract(image_raw,max_control_points,fwhm,mask)

        # Do photometry
        _apertures = CircularAperture(xy_centroids, r=fwhm)
        phot_table = aperture_photometry(_image, _apertures) 
        noise_table = aperture_photometry(_bkg_rms, _apertures)

        brightness = phot_table['aperture_sum'].value
        noise = noise_table['aperture_sum'].value
        snr = brightness/noise # Signal Noise Ratio

        dict_values = xy_centroids,_offset,image_raw,_image,_apertures,brightness,snr,_mask_rectangle
        dict_keys = 'xy','_offset','image_raw','_image','_apertures','brightness','snr','_mask_rectangle'
        info = dict(zip(dict_keys, dict_values))

        return Source(info)      

    def show_image(self,fig_out=None):
        """
        Show raw image

        Parameters:
        fig_out -> [str,optional,default=None] In not None, save the output image to a file defined by fig_out

        """
        if fig_out is None:
            show_image(self.image_raw,origin='lower') 
        else:    
            plot_kwargs = {'figname':fig_out}
            show_image(self.image_raw,origin='lower',**plot_kwargs)         

class Source(object):
    """
    Class Source
    """
    
    def __init__(self,info):  

        self.info = info

        for key in info.keys():
            setattr(self, key, info[key])

    def __repr__(self):
    
        return 'Instance of class Source'

    def show_image(self,fig_out=None):
        """
        Show raw image with stars marked

        Inputs:
            fig_out -> [str,optional,default=None] In not None, save the image to a file defined by fig_out
        """
        if fig_out is None:
            plot_kwargs = {'mark':(self._apertures,'blue')}
        else:
            plot_kwargs = {'mark':(self._apertures,'blue'),'figname':fig_out}
        show_image(self.image_raw,origin='lower',mask_rectangle=self._mask_rectangle,**plot_kwargs)      
 
        