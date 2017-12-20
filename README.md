# CHIMERA
Library of SolarSoft IDL routines to identify and extract coronal holes and their properties from EUV images of the Sun.

Note: This algorithm requires use of the IDL programming language and the SolarSoft libraries

Author: Tadhg M. Garton

This algorithm extracts CH boundaries and properties from 1 HMI lineof-sight magnetogram and 3 AIA EUV images centered at 171, 193, and 211A. Calls of this procedure are made in the form:

chimera, temp = temp, outpath = outpath

where temp describes the current location of the 3 AIA .fits files and the 1 HMI .fits file, and outpath describes the path location of any output files such as segmented images or property files. These files require the naming convention:

'/AIAsynoptic0171.f*'
'/AIAsynoptic0193.f*'
'/AIAsynoptic0211.f*'
'/HMI*mag.f*'
