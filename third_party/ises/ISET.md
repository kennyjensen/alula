# ISET algorithm map

This is an independently worded map of the inspected program, not copied
Fortran. Page numbers refer to the thesis's printed pages. The main driver
on pp.171–173 establishes the order below.

| Stage | Original routine | Required correspondence in this project |
| --- | --- | --- |
| Read geometry and station distributions | READIN | Preserve the distinction between body, inlet/outlet and transverse distributions. Record units, topology and array order. |
| Normalize coordinates and parameterize the contour | NORMIT, SPLNIT | Check normalization, spline knots, LE/TE indexing and contour orientation before meshing. |
| Construct prescribed outlines | OUTLIN | Reconstruct body and stagnation-path station placement first. Derive the common fourth-root Xi labels from the same physical outline. |
| Fill the interior | NORLIN | Reproduce the transverse coordinate interpolation and the isolated-airfoil farfield row ordering. |
| Reconcile contour and grid stations | RESPLI | Account explicitly for contour resampling before comparing boundary coordinates or residuals. |
| Relax both physical coordinates | ELLIP | Original secant metrics, metric updates between lines, tridiagonal solves and source boundary-update order. |
| Initialize flow and optional viscous state | QCALC, BLINIT, OFFSET; ELLIP again when coupled | A prepared mesh is not a converged Euler/BL solution. Preserve mass labels and distinguish initialization from subsequent simultaneous coupling. |
| Write state | OUTPUT | Define a browser-safe state format with explicit geometry, labels, conditions and convergence status. |

The [machine-readable inventory](ROUTINES.json) records all twelve routines
from pp.171–189 and four spline helpers from pp.278–279, with their exact
page ranges and dependencies. OUTLIN is on pp.177–179, RESPLI pp.179–180,
NORLIN pp.180–181 and ELLIP pp.186–189. SPLINE occupies p.278, its scalar
tridiagonal helper SOLV pp.278–279, and SEVAL and DEVAL p.279.

The inspected driver references `STATE.INC` and `ISET.INC`. SPLNIT on p.176
instead names `ISES.INC`; this discrepancy is recorded without silently
correcting it. ELLIP's Thomas elimination is inline and needs no separately
named tridiagonal routine. Reconstructing
the driver alone without its arrays, include files and called routines
would not yield an executable reference. OUTPUT remains unindexed; this
inventory does not establish a complete initializer or full ISES implementation.
READIN's BLADE.DAT OPEN on unit 3 is commented out on p.173, while its
SPOS.DAT OPEN on unit 4 is active on p.174.

RESPLI replaces the blade samples with the boundary-grid samples and then
resplines them. It precedes ELLIP; it is not postprocessing of an already
converged mesh. This geometry change must be explicit in a faithful
initializer comparison.

## Verified spline reconstruction

SPLINE on p.278 constructs a natural cubic spline in the supplied parameter
`S`. Its output `XP` stores first derivatives, not second derivatives; the
contour arrays `XPB` and `YPB` therefore contain dx/ds and dy/ds. For knot
derivatives m and adjacent secant slopes delta, its endpoint equations are
`2*m[1] + m[2] = 3*delta[1]` and
`m[N-1] + 2*m[N] = 3*delta[N-1]`. These impose zero second derivatives at
the first and last knots. SEVAL evaluates the resulting cubic Hermite segment;
DEVAL evaluates its analytic first derivative (p.279).

SPLNIT (p.176) and RESPLI (pp.179–180) apply one spline to a smooth contour,
with natural conditions at its two TE endpoints. Coincident TE coordinates
do not make this spline periodic. A sharp LE has a duplicated coordinate
and parameter, and the two branches are splined separately; the duplicate
parameter must not enter one spline system. RESPLI recomputes `S` from the
cumulative straight-line lengths between resampled contour knots, rather
than from exact spline arc length. The rounded-contour reconstruction now
checks these equations against an independent dense first-derivative solve;
the sharp-LE split-spline branch remains unfinished.

## What is already implemented

The [isolated-airfoil reconstruction](../../docs/ISET_SOURCE_INITIALIZATION.md)
now implements NORMIT/SPLNIT/OUTLIN/NORLIN/RESPLI for a rounded LE and closed
TE, then enters source ELLIP through the actual solver chart. Its tiny case
has 108 positive initial/final quads; 18 ISET and 38 SLOR checks pass. No
original Fortran execution, multielement or physical acceptance is claimed.

The [elliptic kernel](../../src/geometry/elliptic-streamtube-grid.js) contains
the separately named `giles-1985` secant stencil, shared Thomas factorization
for x/y, linewise metric updates and fourth-root station labels. Its literal
horizontal exterior-row copying and the more general `giles-indexed-y`
condition are separately tested. The
[multielement adapter](../../src/euler/streamtube-elliptic-initializer.js)
shares a primary-body Xi array across passages; that choice is a declared
multielement extension, not recovered modern MSET source.

The browser currently uses potential-block boundary preparation and fixed
farfield nodes. It is not a complete reconstruction of OUTLIN/NORLIN/RESPLI.
Its default mesh now passes the public geometry regression; physical and
settings-envelope acceptance remain open.

## Physical station spacing before ELLIP: first implementation

Giles's OUTLIN uses supplied normalized physical inlet locations (`SINL`,
printed pp.174–177). The current upstream map instead uses potential with
a positive endpoint derivative. Near stagnation, distance varies as the
square root of potential difference. This produces incompatible spacing
when joined to a surface parameterized by arc length.

Saved main-LE boundary measurements are:

| Segment | Physical length |
| --- | ---: |
| Incoming cut, station 31 to 32 | 0.04726142254 |
| First upper surface interval | 0.00606211976 |
| First lower surface interval | 0.01027364967 |

The upper length ratio's fourth root is 0.5984521487, exactly the assigned
Xi interval ratio. This discontinuity exists before SLOR. Separately,
independent block-count normalization jumps the underlying rank interval
by a factor of 2.57556.

The first optional physical-x implementation matched discrete x intervals
across independently distributed cut blocks. Its default case had 16 invalid
unsmoothed cells despite 2,709 positive smoothed quads; that report is retained.
The [current full-outline reconciliation](../../docs/ISET_PHYSICAL_STATIONS.md)
joins cut-x and wall-arc maps and shares slope variables between cut copies
before collecting resolution demands. All 2,583 initial/final quads are now
positive. The [unchecked-mesh correction](../../docs/UNSMOOTHED_GRID_CORRECTION.md)
now connects this preparation to the GUI. Maximum smoothing movement still
reaches 0.4789 chord, so the reliability gate remains open.
These are declared multielement extensions, not
literal OUTLIN reconstruction or exact streamline-arclength parameterization.

The isolated-airfoil OUTLIN/NORLIN geometry and RESPLI order now have a
reconstructed baseline. The latest saved-node audit finds no invalid traced
cells, but linear interpolation of its same boundaries still has 36 invalid
cells. Resolve source outer-boundary construction and station transitions
before another solve. Keep the original
fourth-root Xi rule unchanged. Match distance increments, not tangents:
the incoming stagnation streamline and wall meet at a real turning corner.

The source line equations, copy order and one complete small initialization
are now checked. Next extend the boundary correspondence to multielement
passages with explicit shared-boundary and continuity checks. Further
alternative smoother development is paused.

## Farfield distinction to preserve

NORLIN (pp.180–181) initially places the exterior rows on straight lines whose
y values can vary with station index. ELLIP (pp.186–187) retains each indexed
exterior y value and copies only x from the adjacent interior row. This keeps
crossline end segments vertical at convergence, but generally changes the
physical exterior curve when those y values are nonconstant. Restricting this
operation to horizontal rows was a restriction of our earlier implementation,
not a requirement in the printed ELLIP listing. Its vertical derivative
condition is not normal Neumann data on a prescribed sloped curve.

A normal projection onto a curved farfield is a different condition. Our tested
normal-curve experiment created a corner incompatibility with the fixed
vertical inlet and failed the spacing screen, so it is not the GUI default.
See [the retained result](../../docs/CURVED_FARFIELD_SLOR.md). A vertical
intersection with a fixed curve is closer to the original direction, but
reevaluating y would still be an explicit extension, requiring an oblique
Xi boundary condition in the independent reference.
