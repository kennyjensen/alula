# NLR main/flap geometry

Import `blade.nlr7301-wind` through the workbench's coordinate-file control.
It contains the main element and flap from NASA's archived WIND example,
with their original coordinates and spacing. Use reference chord 1 for those
coordinates. Source, hash, extraction ranges and units are recorded in
`nlr7301-wind-provenance.json`.

This is a geometry example. The archive's case page, input file and force-file
labels disagree about incidence and turbulence treatment, and the page has no
experimental data or refinement study. Its flow outputs are therefore not used
as validation targets here.

An incompressible geometry test at incidence 0°, Re=1 million and trips
0.05/0.05 converges with four BL surfaces and two curved wakes. Those settings
are an alula test case, not a reconstruction of the archived WIND run.

Regenerate the extraction from the pinned source grid with:

```sh
node scripts/reference/extract-nlr-geometry.js
```

An optional argument supplies an already downloaded `nlrflap.x` file instead.
The source hash must match before any generated file is written.
