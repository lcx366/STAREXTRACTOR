"""
slrfield package

This package is an archive of scientific routines for data processing related to SLR(Satellite Laser Ranging).   
Currently, operations on SLR data include:

1. Download CPF(Consolidated Prediction Format) ephemeris files automatically from **CDDIS**(Crustal Dynamics Data Information System) or **EDC**(EUROLAS Data Center);
2. Parse the CPF ephemeris files;
3. Calculate the position of targets in GCRF;
4. Predict the azimuth, altitude, distance of targets, and the time of flight for laser pulse etc.;
"""

from .classes import AstroImage