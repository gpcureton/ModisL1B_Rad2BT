      REAL FUNCTION BRITE_M(V, R)

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c    Compute brightness temperature given monochromatic Planck radiance
c    (Radiance units: milliWatts per square meter per steradian per
c    inverse centimeter)
c
c!INPUT PARAMETERS:
c    V (REAL)          Wavenumber (inverse centimeters)
c    R (REAL)          Monochromatic Planck radiance (milliWatts per
c                      square meter per steradian per
c                      inverse centimeter)
c
c!OUTPUT PARAMETERS:
c    BRITE_M (REAL)    Brightness temperature (Kelvin)
c
c!REVISION HISTORY:
c
c!TEAM-UNIQUE HEADER:
c    Liam.Gumley@ssec.wisc.edu
c
c!END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Include files ...........................
      !include 'fundamental_constants.inc'
c ... Planck constant (Joule second)
      double precision h
      parameter (h = 6.62606876d-34)

c ... Speed of light in vacuum (meters per second)
      double precision c
      parameter (c = 2.99792458d+08)

c ... Boltzmann constant (Joules per Kelvin)      
      double precision k
      parameter (k = 1.3806503d-23)

c ... Derived constants      
      double precision c1, c2
      parameter (c1 = 2.0d+0 * h * c * c)
      parameter (c2 = (h * c) / k)
c ... Finish include files ......................

c ... Arguments
      real v, r

c ... Local variables
      double precision vs

c ... Set default return value
      brite_m = -1.0
      
c ... Check input parameters and return if they are bad
      if (v .le. 0.0 .or. r .le. 0.0) return
                  
c ... Convert wavenumber to inverse meters
      vs = 1.0d+2 * dble(v)
      
c ... Compute brightness temperature
      brite_m = sngl(c2 *
     &  vs / log(c1 * vs**3 / (1.0d-5 * dble(r)) + 1.0d+0))
      
      END
