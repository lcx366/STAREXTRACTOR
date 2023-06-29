# Welcome to the STAREXTRACTOR package

[![PyPI version shields.io](https://img.shields.io/pypi/v/starextractor.svg)](https://pypi.python.org/pypi/starextractor/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/starextractor.svg)](https://pypi.python.org/pypi/starextractor/) [![PyPI status](https://img.shields.io/pypi/status/starextractor.svg)](https://pypi.python.org/pypi/starextractor/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/STAREXTRACTOR.svg)](https://GitHub.com/lcx366/STAREXTRACTOR/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/STAREXTRACTOR/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/STAREXTRACTOR.svg)](https://github.com/lcx366/STAREXTRACTOR/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/starextractor/badge/?version=latest)](http://starextractor.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/starextractor.svg?branch=master)](https://travis-ci.org/lcx366/starextractor)

This package is an archive of scientific routines for data processing related to the source extraction from an astronomical image.
Currently, operations on source extraction include:

1. Read an astronomical image in `.fits` format or in generic image format, such as `.bmp`;
2. Extract the centroid coordinates of the star spots and do photometry using `photutils`;
3. Show raw image with star spots marked;

## How to Install

On Linux, macOS and Windows architectures, the binary wheels can be installed using `pip` by executing one of the following commands:

```
pip install starextractor
pip install starextractor --upgrade # to upgrade a pre-existing installation
```

## How to use

### Read an astronomical image

Images can be in `.fits` format or in generic image format, such as `.bmp`.

```python
>>> from starextractor import AstroImage
>>> imagefile = 'obs/fits/img_00000.fits' 
>>> #imagefile = 'obs/bmp/img_00000.bmp'
>>> #imagefile = 'obs/npy/img_00000.npy'
>>> image = AstroImage.read_image(imagefile)
>>> image = AstroImage.read_image(image_array)
```

Print the raw grayscale image with the origin at the center of the bottom(the first row of array) left pixel.

```python
>>> print(image.image_raw,image.res) # original grayscale image and its resolution
```

### Show the raw image

```python
>>> image.show()
>>> #image.show('figs/raw_image.png') # save image to a file
```

<p align="middle">
  <img src="readme_figs/image_raw.png" width="500" />
</p>

### Extract the centroid coordinates of the star spots and do photometry

Estimate the centroids coordinates, brightness(sum of gray value within an aperture),  and SNR of star spots.

```python
>>> sources = image.find_source(fwhm=12,mask=True)
>>> print(sources.xy,sources.brightness,sources.snr,sources.offset)
```

### Calculate the triangle invariants and construct a 2D Tree; and record the asterism indices for each triangle.

```python
>>> sources.invariantfeatures()
>>> print(sources.invariants,sources.asterisms,sources.kdtree)
```

### Show the extracted sources in image

```python
>>> sources.show()
>>> #sources.show('figs/sources.png') # save image to a file
```

<p align="middle">
  <img src="readme_figs/centroids.png" width="500" />
</p>

## Change log

- **0.1.6 — Jun 29,  2023**
  - Added support for image file formats `.npy` and `numpy array` in function *AstroImage.read_image()*

- **0.1.5 — May 14,  2023**
  - The class `Centriod` is *deprecated*, and the class `Source` is used instead
  - Add method `.invariantfeatures()` to class `Source`, which calculates the triangle invariants and constructs a 2D Tree; and records the asterism indices for each triangle.

- **0.1.0 — Apr 5,  2023**
  
  - The ***starextractor*** package was released.

## Reference

- [photutils](https://photutils.readthedocs.io/en/stable/index.html)
- [Photometry Methods](http://srmastro.uvacreate.virginia.edu/astr313/lectures/photometry/photometry_methods.html)
- [Astroalign](https://astroalign.quatrope.org/en/latest/)