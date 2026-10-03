// SPDX-License-Identifier: GPL-2.0-or-later
// Panel allocation adapts Vibefoil solver_worker.js by Kenny Jensen (GPL-2.0+).
// See xfoil/README.md and the project NOTICE.md for attribution.
import { createMatrix } from './xfoil/arrays.js';
import { scalc, segspl, seval } from './xfoil/spline.js';
import { lefind } from './xfoil/xgeom.js';
import { ncalc, apcalc, ggcalc } from './xfoil/xpanel.js';
import { tecalc } from './xfoil/xfoil.js';
import { buildBlContext, specal, viscal } from './xfoil/xoper.js';

// Geometry-only version of native TECALC. The prescribed base source and
// BL dead-air gap must use these SAME, unnormalized spline derivatives.
// No panel matrix, wake tracing, or flow solve is performed here.
export function panelTrailingEdgeGeometry(points) {
  if (!Array.isArray(points) || points.length < 3 || points.some(p => !Number.isFinite(p.x) || !Number.isFinite(p.y)))
    throw new Error('Supply finite, ordered wetted contour points.');
  const n = points.length, X = Float64Array.from(points, p => p.x), Y = Float64Array.from(points, p => p.y);
  const S = new Float64Array(n), XP = new Float64Array(n), YP = new Float64Array(n);
  scalc(X, Y, S, n);
  for (let i = 1; i < n; i++) if (!(S[i] > S[i - 1])) throw new Error('Wetted contour stations must be distinct.');
  segspl(X, XP, S, n); segspl(Y, YP, S, n);
  const ctx = { N: n, X, Y, XP, YP }; tecalc(ctx);
  if (![ctx.ANTE, ctx.ASTE, ctx.DSTE, ...XP, ...YP].every(Number.isFinite) || ctx.ANTE < 0)
    throw new Error('Invalid native trailing-edge projection; require upper-TE to lower-TE contour orientation.');
  return { width: ctx.ANTE, tangentialProjection: ctx.ASTE, magnitude: ctx.DSTE, sharp: ctx.SHARP,
    center: { x: .5 * (X[0] + X[n - 1]), y: .5 * (Y[0] + Y[n - 1]) },
    upperDerivative: { x: XP[0], y: YP[0] }, lowerDerivative: { x: XP[n - 1], y: YP[n - 1] },
    tangentDerivative: { x: .5 * (-XP[0] + XP[n - 1]), y: .5 * (-YP[0] + YP[n - 1]) },
    convention: 'Native SCALC/SEGSPL/TECALC on the ordered wetted panel contour; derivatives are not normalized.' };
}

export function createPanelContext(points, alpha, { geometryOnly = false } = {}) {
  if (typeof geometryOnly !== 'boolean') throw new Error('Invalid panel geometry-only control.');
  const n = points.length; const nw = Math.floor(n / 12) + 10; const total = n + nw;
  const ctx = { N: n, NW: nw, WAKLEN: 1, PI: Math.PI, QOPI: 1 / (4 * Math.PI),
    HOPI: 1 / (2 * Math.PI), ALFA: alpha * Math.PI / 180, QINF: 1,
    LIMAGE: false, YIMAGE: 0, SHARP: true, ANTE: 0, ASTE: 0, DSTE: 0,
    LWAKE: false, LWDIJ: false, LADIJ: false,
    GAMU: createMatrix(geometryOnly ? 0 : n + 1, 2), QINVU: createMatrix(geometryOnly ? 0 : total, 2),
    AIJ: createMatrix(geometryOnly ? 0 : n + 1, n + 1), BIJ: createMatrix(geometryOnly ? 0 : n + 1, total), AIJPIV: new Int32Array(geometryOnly ? 0 : n + 1) };
  for (const key of ['X', 'Y', 'XP', 'YP', 'S', 'NX', 'NY', 'APANEL', 'SIG', 'DZDM', 'DQDM', 'SNEW']) ctx[key] = new Float64Array(total);
  for (const key of ['QF0', 'QF1', 'QF2', 'QF3', 'DZDG', 'DZDN', 'DQDG']) ctx[key] = new Float64Array(n);
  for (const key of ['GAM', 'GAM_A']) ctx[key] = new Float64Array(n + 1);
  for (const key of ['QINV', 'QINV_A', 'QVIS']) ctx[key] = new Float64Array(total + 1);
  points.forEach((p, i) => { ctx.X[i] = p.x; ctx.Y[i] = p.y; });
  scalc(ctx.X, ctx.Y, ctx.S, n);
  segspl(ctx.X, ctx.XP, ctx.S, n); segspl(ctx.Y, ctx.YP, ctx.S, n);
  ncalc(ctx.X, ctx.Y, ctx.S, n, ctx.NX, ctx.NY);
  const sle = lefind(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, n);
  ctx.XLE = seval(sle, ctx.X, ctx.XP, ctx.S, n); ctx.YLE = seval(sle, ctx.Y, ctx.YP, ctx.S, n);
  ctx.XTE = (ctx.X[0] + ctx.X[n - 1]) / 2; ctx.YTE = (ctx.Y[0] + ctx.Y[n - 1]) / 2;
  ctx.CHORD = Math.hypot(ctx.XTE - ctx.XLE, ctx.YTE - ctx.YLE);
  tecalc(ctx);
  if (geometryOnly) ctx.SLE = sle;
  else { apcalc(ctx); ggcalc(ctx); }
  return ctx;
}

// Internal raw contexts are exposed only for native parity and derivative tests.
export function runCoupled(points, { alpha = 0, reynolds = 1e6, ncrit = 9,
  trips = [1, 1], maxIterations = 100, onIteration } = {}) {
  const panel = createPanelContext(points, alpha);
  specal(panel, alpha * Math.PI / 180);
  const bl = buildBlContext(points.length, panel, ncrit);
  bl.XSTRIP[1] = trips[0]; bl.XSTRIP[2] = trips[1]; bl.onIteration = onIteration;
  const solution = viscal(bl, panel, alpha * Math.PI / 180, reynolds, { maxIter: maxIterations });
  return { panel, bl, ...solution };
}
