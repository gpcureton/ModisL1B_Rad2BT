"""
   DESCRIPTION:
       Compute brightness temperature for a MODIS infrared band
       on Terra or Aqua.
   
       Spectral responses for each IR detector were obtained from MCST:
       ftp://ftp.mcst.ssai.biz/pub/permanent/MCST in directories
       PFM_L1B_LUT_4-30-99 (Terra) and FM1_RSR_LUT_07-10-01 (Aqua).
   
       An average spectral response for all detectors in each band was
       computed. The detector-averaged spectral response data were used
       to compute the effective central wavenumbers and temperature
       correction coefficients included in this module.
   
       NOTE: The plaform name ('Terra' or 'Aqua') is passed to this
       function via the common block defined in 'platform_name.inc'.
   
   INPUT PARAMETERS:
       RAD (REAL)      Planck radiance (units are determined by UNITS)
       PLATFORM (LONG) Platform Name: Terra=0 or Aqua=1
       BAND (LONG)     MODIS IR band number (20-25, 27-36)
       UNITS (LONG)    Flag defining radiance units
                       0 => milliWatts per square meter per
                            steradian per inverse centimeter
                       1 => Watts per square meter per
                            steradian per micron
   
   OUTPUT PARAMETERS:
       MODIS_BRIGHT  Brightness temperature (Kelvin)
                     Note that a value of -1.0 is returned if
                     RAD .LE. 0.0, or BAND is not in range 20-25, 27-36.
   
   !REVISION HISTORY:
       Liam.Gumley@ssec.wisc.edu (Fortran 77 implementation)
       geoff.cureton@ssec.wisc.edu (python wrapper)
   
   !TEAM-UNIQUE HEADER:
       Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
"""
#import mod_br
