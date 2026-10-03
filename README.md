# alula

This repository contains a browser-based airfoil analysis laboratory.  It has a
panel solver, integral boundary layers, and an experimental conservative
streamtube Euler plus boundary-layer solver.  The application is intended for
development and numerical investigation; convergence and a small residual do
not by themselves establish physical accuracy.

## Development

```sh
npm install
npm run dev
```

Build the browser application with `npm run build`. Run the tests with
`npm test` and the project validation set with `npm run validate`.

`docs/` contains local documentation, historical results and large archived
fixtures. It is ignored by Git and is not included in the static build. A fresh
checkout can build and run tests that use fixtures under `tests/fixtures/`;
tests that read archived files under `docs/` require a separate copy of that
local archive. Generated reports in `docs/`, `results/`, `artifacts/`,
`ci-artifacts/`, `test-results/` and `playwright-report/` are also ignored.

The RAE status table is generated from browser-form verification receipts. To
verify the default inviscid 64×9 case and update its local status:

```sh
RAE_APP_RECEIPT=/tmp/rae-app.json npm run verify:rae:app
npm run status:rae -- --record=/tmp/rae-app.json
```

After solver changes, `npm run status:rae` marks stale local receipts unresolved
until reverified. CI runs the 64×9 browser regression and retains its receipts
as a workflow artifact. The generated status table stays local under `docs/`.

See the [roadmap](ROADMAP.md) and [third-party notices](NOTICE.md) for project
scope and attribution.
