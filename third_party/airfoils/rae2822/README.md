# RAE 2822 measured coordinates

The retained [geom.txt](geom.txt) is the unmodified coordinate table from the
[NASA NPARC RAE 2822 validation archive](https://www.grc.nasa.gov/WWW/wind/valid/raetaf/raetaf.html).
NASA identifies it as the measured ordinates from Table 6.1 of Cook, McDonald
and Firmin, *Aerofoil RAE 2822—Pressure Distributions, and Boundary Layer and
Wake Measurements*, AGARD AR 138 (1979). Coordinates use the airfoil chord as
the reference length.

NASA's accompanying [geom.f90](geom.f90) explicitly negates the tabulated
lower-surface ordinate. The browser preset follows that convention: reverse
the 65 upper points from TE to LE, then append the 64 remaining lower points
from LE to TE. This gives **129 points / 128 panels**, a repeated sharp TE,
and every original nonzero ordinate. Signed zero is canonicalized to zero.
There is no resampling, scaling, smoothing or altered trailing-edge closure.

Original download URLs:

- [Coordinate table](https://www.grc.nasa.gov/WWW/wind/valid/raetaf/geom.txt)
- [Fortran coordinate converter](https://www.grc.nasa.gov/WWW/wind/valid/raetaf/geom.f90)

SHA-256:

```text
geom.txt  0bef162c8778247a4746f02a1b487f33083a1ae80fdc82ea58afc259c87da736
geom.f90  631d2d0be0b2cb9fa87dc2bf49f4f65a658794115a5649b0f68749312f33700e
```

The original files contain no explicit license statement. Attribution is
retained; the project's code license is not assigned to these upstream data.
Adding this geometry preset does not validate a computed flow against the
experiment or set its Mach number, incidence, Reynolds number or transition.
