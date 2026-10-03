# Roadmap

The immediate work is to make the numerical methods dependable and measurable.

1. Make the RAE 2822 transonic inviscid and coupled cases converge reliably on
   supported grids, with reproducible iteration checkpoints.
2. Improve streamtube-grid initialization, smoothing, and state transfer
   without weakening conservative residuals.
3. Validate panel and Euler-plus-boundary-layer results against published and
   independently reproducible reference cases.
4. Keep solver status, diagnostics, and plots clear enough to distinguish a
   converged result from a provisional iterate.

Historical experiment logs were intentionally removed during repository
cleanup.  Tests and the fixtures they require remain the source of regression
evidence.
