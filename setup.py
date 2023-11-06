from setuptools import setup,find_packages 

setup(
    name='starextractor',
    version='0.2.0',
    description='A package for extracting stars and doing photometry from an astronomical image',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/STAREXTRACTOR',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['source extraction','photometry','psf','daofind','fwhm'],
    python_requires = '>=3.10',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        'photutils',
        'astropy>=4.3.1',
        'Pillow',
        'scipy',
        'matplotlib'
        ],
)
