# ISES / ISET source-reference workspace

Requested 7 September 2026. Reference: Michael B. Giles, *Newton solution
of steady two-dimensional transonic flow*, MIT PhD thesis, 1985, GTL 186.

- [MIT repository record](https://hdl.handle.net/1721.1/15250)
- [Original scanned PDF](https://dspace.mit.edu/server/api/core/bitstreams/24779c2e-66a8-48f9-8207-53ab301c004f/content)
- [Source identity and reuse status](SOURCE.json)
- [ISET routine map and implementation plan](ISET.md)
- [Verified routine pages, calls and include files](ROUTINES.json)
- [ISES Newton update, SMOVE and MSES source-to-code audit](../../docs/ISES_NEWTON_COMPARISON.md)
- [First independent physical-station implementation and its limits](../../docs/ISET_PHYSICAL_STATIONS.md)
- [Reconstructed isolated ISET initialization and solver-chart checks](../../docs/ISET_SOURCE_INITIALIZATION.md)

**No Fortran transcription is present.** MIT's item record restricts
reproduction/distribution without written permission; no permissive source
license was found for this listing. The repository's GPL does not relicense
the thesis. These files contain our source index and algorithm notes.

The listing starts after the appendix introduction on printed p.170. ISET's
main program occupies pp.171–173 and its initializer routines continue
through p.189. The appendix also contains the flow solver and boundary-layer
routines; an ISET-only reconstruction is not the whole ISES program.

Our implementation target is the ISET initialization sequence, including
its boundary construction and parameterization before ELLIP. Existing
JavaScript implementations remain under `src/`; they must not be labeled
as compiled or executed original Fortran. A supplied listing can be handled
as a separate transcription with page-level uncertainty records. Any
redistributable upstream source must retain its own license and provenance.

The independent JavaScript is executable now:

- [ISET geometry stages](../../src/geometry/iset-airfoil-initializer.js):
  NORMIT → SPLNIT → OUTLIN → NORLIN → RESPLI for one rounded-LE, closed-TE foil.
- [Quadrilateral solver-chart adapter](../../src/euler/iset-initializer.js):
  initial mesh snapshots and the source ELLIP boundary/stencil path.

From the repository root:

```sh
npm run test:streamtube:iset
npm run diagnose:streamtube:iset-airfoil
```

The recorded 108-cell case checks geometry and shared cuts. It does not
establish converged Euler/BL physics or multielement initialization. The
[saved multielement correspondence audit](../../docs/ISET_PHYSICAL_STATIONS.md#next-bounded-comparison)
identifies the remaining boundary-distribution conflict without a flow solve.
