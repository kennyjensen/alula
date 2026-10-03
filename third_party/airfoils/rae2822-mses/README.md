# RAE 2822 — MSES Case 13a sample geometry

`blade.rae` is the exact public coordinate file linked from
[Mark Drela's MSES website](https://web.mit.edu/drela/Public/web/mses/),
downloaded from [its original URL](https://web.mit.edu/drela/Public/web/mses/blade.rae).
It accompanies the site's [RAE 2822 Case 13a example](https://web.mit.edu/drela/Public/web/mses/rae2822.pdf).
The retained file is copied byte-for-byte from this repository's
[source receipt](../../../docs/rae2822/mses-website-reference/SOURCE.json).

The preset keeps all 129 supplied points, including the repeated sharp
trailing-edge endpoint, in the original `c=1` coordinate frame. It performs
no resampling, smoothing, shifting or sign change. The source's domain
header is `x=[-2,3], y=[-3,3.5]`; it is recorded as provenance and does not
silently replace the GUI mesh settings or farfield construction.

This is the sole visible **RAE 2822** preset, with internal ID
`rae2822-mses`. The unchanged NASA measured-coordinate dataset remains
internally available as `rae2822` for historical validation and saved inputs;
its ID is not remapped and it is not a second menu selection. The MIT upper ordinates differ by up to `0.00031 c` and the lower
ordinates by up to `0.00027 c`. No reason for those differences is assumed.
The [reference audit](../../../docs/RAE2822_MSES_REFERENCE.md) identifies
the published MSES conditions, formulation and experimental comparison.
Choosing the preset changes geometry and reference chord only; it does
not apply MSES flow controls, trips, drag coefficients or a saved solution.

Source SHA-256:
`f12c6730d2be65b9f439cb0885e31f87fcf4a53a57636687fae1d7459fcc6edd`.
Attribution: Mark Drela, MIT, public MSES sample input. The downloaded
coordinate file contains no explicit data-license declaration. No GPL
license or ownership claim is assigned to the upstream numerical data.
Geometry import and sharp-edge compatibility are not convergence or
aerodynamic-accuracy validation.
