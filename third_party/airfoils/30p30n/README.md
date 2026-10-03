# McDonnell Douglas 30P-30N source geometry

`30p-30n.dat` preserves the complete, unchanged bytes of the [UIUC Airfoil Data Site coordinate file](https://m-selig.ae.illinois.edu/ads/coord_seligFmt/30p-30n.dat). Its [updates-directory copy](https://m-selig.ae.illinois.edu/ads/coord_updates/30p-30n.dat) was byte-identical when retrieved on 2026-09-13. [UIUC's update history](https://m-selig.ae.illinois.edu/ads_history.html) credits Brent Pomeroy with contributing this geometry and rigging information on 2016-04-08. SHA-256: `1b973fea5899d4558cad9643c823795049ad98450a1c472599a9f01cf11f00ad`.

This is an exact copy of that public coordinate sample. We have not established that every sampled point equals the original wind-tunnel manufacturing surface. It must not be described as a reconstructed NASA master geometry or as the modified BANC acoustic geometry.

The three comment-delimited blocks are `# Slat`, `# Main Element`, and `# Flap`, in that order. Their coordinates already describe the deployed assembly: use **x/c, y/c with stowed reference chord c = 1**, without element rotations, translations, normalization by element chord, or coordinate resampling. The header specifies:

| Rigging | Slat | Flap |
|---|---:|---:|
| Deflection magnitude | 30° | 30° |
| Gap, % stowed chord | 2.95 | 1.27 |
| Overhang, % stowed chord | −2.50 | +0.25 |

The NASA-authored [AIAA 2002-0845 paper, §2, printed p. 2](https://fun3d.larc.nasa.gov/papers/aiaa2002-0845.pdf) documents this rigging for the MDA LB546 configuration: LS12 slat, W10BB main element, and F22 flap, with a 22-inch (0.5588 m) stowed chord in Langley's LTPT. Its rounded deployed leading-edge x/c locations (−0.0854, 0.0438, 0.8715) agree with the sampled minima (−0.085382, 0.043807, 0.871510). The reported study used nominal M = 0.2 and Re = 9 million. Those are experimental context, not automatic settings or an aerodynamic validation of this preset. The paper also documents three-dimensional tunnel effects near maximum lift. [NASA catalog record](https://ntrs.nasa.gov/citations/20020012701).

## Trailing-edge limitation

The raw sample contains **201 slat points, 221 main points, and 242 flap points**. All coordinates are finite and all counts fit the current per-contour limit. Slat and main have repeated trailing-edge points and pass `prepareContour` unchanged. The flap is open:

- First point: (1.128307, −0.145799); next: (1.127975, −0.145521).
- Last point: (1.130859, −0.140497); previous: (1.130625, −0.140133).
- Endpoint separation: **0.005884208357969709 c**, approximately 13.6 times either terminal surface interval.

Appending the first point would add a finite base segment. The resulting base/surface angles at its endpoints are approximately 75.76° and 121.56°, rather than the sharp cusp expected by the solver. Thus the incompatibility is not fixed by supplying a missing duplicate point. We have not added a base panel, moved the endpoints, filled the slat cove, or sharpened the flap. The preset preserves the supplied geometry for preview/export and blocks solve/mesh construction until this geometric restriction is addressed. A future sharp-edge approximation needs a distinct name and a documented change, separate from this source.

`geometry-audit.json` records exact bounds/endpoints and current contour-validation results; `audit.mjs` reproduces it without a mesh or flow calculation. An audit `passed` means the source inspection matched its expectations, not that a 30P30N aerodynamic solution passed.

## Independent workshop files and rights

`workshop/3-element-airfoil.IGS` and `.iges` are unchanged [original high-order CFD workshop CAD downloads](https://cfd.ku.edu/hiocfd/). [Workshop case C3.1](https://cfd.ku.edu/hiocfd/case_c3.1.html) specifies c = 0.5588 m, M = 0.2, α = 16°, Re = 9 million and fully turbulent flow. These CAD files contain spline surfaces and explicit finite base lines; they are a separate representation and are not used to replace any UIUC coordinate. Their exported unit labels differ and have not been used for an automatic conversion. No pointwise equivalence is claimed.

No explicit coordinate/CAD license was found in the downloaded files or reviewed publisher pages. UIUC's pages carry the Applied Aerodynamics Group copyright notice. Attribution and original bytes are retained; the data are **not relabeled GPL or public domain**. NASA's catalog marks the cited paper “Public Use Permitted”; that paper notice does not establish the licensing of separately hosted coordinate files. Full URLs, source hashes, interpretation and remaining limitations are recorded in `provenance.json`.
