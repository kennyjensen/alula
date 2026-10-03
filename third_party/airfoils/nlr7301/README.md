# NLR 7301 with flap — SAAB/FLOWNET 2.6% gap geometry

This directory retains the public SAAB structured-grid archive linked by the
[ERCOFTAC UFR 3-01 description](https://www.kbwiki.ercoftac.org/w/index.php?title=UFR_3-01_Description).
The preset contains the exact two solid wall contours from `t/nlrgap26.dat`
inside that archive. It is available for viewing and coordinate download.
**Both elements have finite trailing-edge bases, which the current sharp-edge
solver cannot represent. Build and solve are disabled for this exact preset.**

The supplied reference frame has chord `c = 1`; no geometric transform is
applied. The [ERCOFTAC test-case description](https://www.kbwiki.ercoftac.org/w/index.php?title=UFR_3-01_Test_Case)
specifies a 20° flap deflection, 5.3% overlap and two gap settings, 1.3% and
2.6% of chord. This preset selects the 2.6% archive member. Its reported test
conditions are incidence 13.1°, Mach 0.185 and Reynolds number 2.51 million;
the dimensional reference chord is 0.57 m. Those are experimental context,
not preset flow controls or an assertion that this solver reproduces the data.
The shaped main-element shroud does not permit retraction of the flap.

Drela's [MSES page](https://web.mit.edu/drela/Public/web/mses/) links an
[NLR 7301 flapped example](https://web.mit.edu/drela/Public/web/mses/nlf26.pdf).
Its matching coordinate input is not linked. The public filenames
`blade.nlr7301`, `blade.nlf26`, `blade.nlr` and `blade.nlf` returned HTTP 404
on 13 September 2026. The preset is explicitly the SAAB/FLOWNET grid geometry;
pointwise identity with Drela's unpublished example input is not established.

## Exact extraction and retained bases

The archive is downloaded unchanged from
[ERCOFTAC's linked object](https://kbwiki-images.s3.amazonaws.com/5/5f/UFR3-01_nlr7301grid.zip).
`NLRdescr.pdf` is copied unchanged from the same archive. It documents the
nine-block file format: block dimensions followed by `x,y` pairs with the
`i` index varying fastest. `extract.py` reads exactly that format.

The extractor enumerates block-boundary segments and cancels pairs whose
endpoint coordinates are exactly identical. This leaves 1,920 unpaired edges
with degree two at every vertex and exactly three closed loops: the main
element (456 distinct points), flap (248) and exterior domain boundary
(1,216). Only the two solid loops enter the preset. Their orientation is made
counterclockwise and each starts at the supplied upper trailing-edge corner.
A duplicate endpoint closes each coordinate list, giving 457 and 249 entries.
No existing point is removed, resampled, shifted or fitted.

The main base endpoints are `(0.9436, 0.01499)` and `(0.9436, 0.01410)`;
the base endpoint distance is `0.00089 c`. The flap endpoints are
`(1.201771, -0.1037033)` and `(1.201377, -0.104784)`, separated by
`0.001150281917618465 c`. All 32 source mesh segments on each base are
retained. These are actual finite bases, not missing repeated closing points.
The existing contour preparer rejects both upper-corner starts because they
do not meet its sharp trailing-edge opening-angle assumption. The preset's
explicit unsupported reason additionally blocks browser build/solve.

`geometry.json` records every source block-side range and extracted point;
`nlr7301-gap26.dat` contains the same points in the repository's MSES-style
coordinate format. Its four-number domain header is a file-format display
default, not an extracted SAAB farfield or solver setting. Reproduce with:

```sh
python3 third_party/airfoils/nlr7301/extract.py
```

## Attribution and hashes

Source grid: SAAB Military Aircraft, made available for FLOWNET through
ERCOFTAC; UFR description contributor Jan Vos, CFS Engineering SA.
The current ERCOFTAC pages display © ERCOFTAC 2004 and a
`CC BY 4.0 (AI/ML-training & TDM reserved)` footer. The archive itself does
not contain a separate data-license declaration. No GPL license or ownership
claim is assigned to the upstream numerical data; the original extractor
code is GPL-2.0-or-later.

* Original archive SHA-256: `28aa28b47d981681df3b6300a41a6a9e1c029514ca28cb91bda699b32b7bfc62`
* Original `t/nlrgap26.dat` SHA-256: `dd930f6761a24a00d9a97925cc7967045572f43d0162381099ff030795d442f5`
* Extracted coordinate file SHA-256: `96376e8ee665e5c63df8f07467ad7f4d3fd8557f63f36e74712e3bd493c833ff`

The original archive also retains the 1.3% case, which is not selected or
transformed into the 2.6% preset. No experimental or CFD accuracy claim is
made from coordinate extraction alone.
