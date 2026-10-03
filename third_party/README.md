# Third-party material

alula is licensed GPL-2.0-or-later; this does not replace upstream licenses.
See [../NOTICE.md](../NOTICE.md) for the component and data attributions.

Place newly vendored upstream source and reference material in a named directory
here. Preserve upstream copyright and license files unchanged. Add a README
recording author, source URL, version or revision, license, files actually used,
and any local changes; record original-file hashes when practical.

Vendored runtime components and the remaining port are located here:

- [XFOIL / Vibefoil JavaScript port](../src/viscous/xfoil/README.md): GPL-2.0-or-later.
- [cdt2d mesher bundle](cdt2d/README.md): MIT components.
- [SuiteSparse KLU WebAssembly](klu/README.md): LGPL and BSD components, with complete selected sources and linked runtime notices.

Keep notices beside these components and include them in distributions. Adapted
files retain upstream notices and identify modifications. Reference documents
and geometry data retain their own rights; do not assign them the application's
license by default. Only selected references and benchmark coordinate documents
are copied by `scripts/build.js`, not this entire reference collection.
