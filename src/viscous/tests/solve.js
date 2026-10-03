// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour } from '../../geometry/airfoil.js';
import { runCoupled } from '../context.js';
import { ueset } from '../xfoil/xbl.js';
import { clcalc } from '../xfoil/xfoil.js';
import { seval } from '../xfoil/spline.js';

export const coupledLimitations = 'Single-element, incompressible XFOIL integral boundary-layer model with global Newton coupling. Native code parity is verified; multielement viscous flow, Euler shocks, and post-stall predictions are not validated.';

export function solveCoupled({ elements, alpha = 0, mach = 0, reynolds = 1e6, ncrit = 9,
  trips = [1, 1], maxIterations = 100, referenceChord = 1,
  momentReference = { x: referenceChord / 4, y: 0 } }, { onIteration } = {}) {
  if (!Array.isArray(elements) || elements.length !== 1) throw new Error('Coupled boundary layers currently require one element. Select a single airfoil or use potential flow for an assembly.');
  if (mach !== 0) throw new Error('The validated coupled baseline requires Mach 0.');
  if (!Number.isFinite(alpha) || Math.abs(alpha) > 20) throw new Error('The coupled baseline accepts incidence between −20° and 20°.');
  if (!Number.isFinite(reynolds) || reynolds < 5e4 || reynolds > 5e7) throw new Error('Enter a Reynolds number between 50,000 and 50,000,000, based on the reference chord.');
  if (!Number.isFinite(ncrit) || ncrit < 0 || ncrit > 14) throw new Error('Ncrit must be between 0 and 14.');
  if (!Array.isArray(trips) || trips.length !== 2 || trips.some(v => !Number.isFinite(v) || v <= 0 || v > 1)) throw new Error('Transition trips must be in (0, 1] of the element chord.');
  if (!Number.isInteger(maxIterations) || maxIterations < 1 || maxIterations > 200) throw new Error('Use 1–200 Newton iterations.');
  if (!Number.isFinite(referenceChord) || referenceChord <= 0 || !Number.isFinite(momentReference?.x) || !Number.isFinite(momentReference?.y)) throw new Error('Invalid reference chord or moment reference.');
  const points = prepareContour(elements[0].points);
  if (points.length < 81 || points.length > 401) throw new Error('The coupled baseline requires 80–400 panels. Use at least 160 for transition resolution.');
  // XFOIL coefficients and Re use one coordinate unit as reference length.
  const normalized = points.map(p => ({ x: p.x / referenceChord, y: p.y / referenceChord }));
  const raw = runCoupled(normalized, { alpha, reynolds, ncrit, trips, maxIterations, onIteration });
  const { panel: p, bl: b, qvis, qinv, history } = raw;
  const surfaces = []; const wake = [];
  for (let side = 1; side <= 2; side++) {
    const stations = [];
    for (let k = 2; k <= b.NBL[side]; k++) {
      const i = b.IPAN[k][side] - 1; const isWake = k > b.IBLTE[side];
      const station = { index: i, x: p.X[i] * referenceChord, y: p.Y[i] * referenceChord,
        s: b.XSSI[k][side] * referenceChord, ue: b.UEDG[k][side],
        theta: b.THET[k][side] * referenceChord, deltaStar: b.DSTR[k][side] * referenceChord,
        h: b.DSTR[k][side] / b.THET[k][side], cf: 2 * b.TAU[k][side],
        amplification: k < b.ITRAN[side] ? b.CTAU[k][side] : null,
        ctau: k >= b.ITRAN[side] ? b.CTAU[k][side] : null,
        regime: isWake ? 'wake' : k >= b.ITRAN[side] ? 'turbulent' : 'laminar' };
      if (!isWake) {
        station.displacement = { x: station.x + station.deltaStar * p.NX[i], y: station.y + station.deltaStar * p.NY[i] };
        stations.push(station);
      } else if (side === 2) wake.push(station);
    }
    const str = b.SST + (side === 1 ? -1 : 1) * b.XSSITR[side];
    surfaces.push({ side: side === 1 ? 'upper' : 'lower', transition: b.XOCTR[side],
      transitionPoint: { x: referenceChord * seval(str, p.X, p.XP, p.S, p.N),
        y: referenceChord * seval(str, p.Y, p.YP, p.S, p.N) }, stations });
  }
  const finite = [...surfaces.flatMap(s => s.stations), ...wake].every(s =>
    [s.x, s.y, s.ue, s.theta, s.deltaStar, s.h, s.cf].every(Number.isFinite) && s.theta > 0 && s.h > 1 && s.ue > 0);
  if (!finite || ![b.CL, b.CM, b.CD, b.CDF, b.RMSBL, ...qvis].every(Number.isFinite)) throw new Error('The coupled solution became nonfinite or inadmissible. Refine the panels or reduce incidence.');
  // Independently evaluate the displacement/edge-velocity interaction equation.
  const saved = b.UEDG.map(row => row.slice());
  ueset(b); let couplingResidual = 0;
  for (let side = 1; side <= 2; side++) for (let k = 2; k <= b.NBL[side]; k++) couplingResidual = Math.max(couplingResidual, Math.abs(b.UEDG[k][side] - saved[k][side]));
  b.UEDG = saved;
  const coeffs = clcalc(p.N, p.X, p.Y, p.GAM, p.GAM_A, b.ALFA, 0, 1,
    momentReference.x / referenceChord, momentReference.y / referenceChord);
  if (![coeffs.cl, coeffs.cm, coeffs.cdp].every(Number.isFinite)) throw new Error('The force or moment calculation became nonfinite. Check reference coordinates.');
  const cp = points.map((point, i) => ({ ...point, cp: 1 - qvis[i + 1] ** 2,
    cpInviscid: 1 - qinv[i + 1] ** 2, qt: qvis[i + 1] }));
  const warnings = [];
  if (!raw.converged) warnings.push(`Newton iteration limit reached (${maxIterations}); coefficients are withheld.`);
  if (points.length < 161) warnings.push('Coarse surface spacing can underresolve transition and separation bubbles. Refine to at least 160 panels.');
  if (Math.abs(alpha) > 4) warnings.push('Outside the initial incidence validation range. Strong separation and stall require further validation.');
  if (reynolds < 5e5 || reynolds > 3e6) warnings.push('Outside the initial Reynolds-number validation range (0.5–3 million).');
  if (surfaces.some(s => s.stations.some(v => v.cf < 0))) warnings.push('Local separation is predicted. Check panel refinement and boundary-layer thickness.');
  const converged = raw.converged && couplingResidual < 1e-5;
  if (raw.converged && !converged) warnings.push('The edge-velocity interaction residual remains above tolerance; coefficients are withheld.');
  return { model: 'xfoil-coupled-incompressible', status: converged ? 'solved' : 'unconverged',
    alpha, mach, reynolds, ncrit, trips, referenceChord, momentReference, panelCount: points.length - 1,
    cl: converged ? coeffs.cl : null, cm: converged ? coeffs.cm : null, cd: converged ? b.CD : null,
    cdf: converged ? b.CDF : null, boundaryLayer: { surfaces, wake }, history,
    elements: [{ name: elements[0].name ?? 'Airfoil', points, cp }],
    diagnostics: { rmsUpdate: b.RMSBL, updateTolerance: 1e-4, couplingResidual, couplingTolerance: 1e-5,
      iterations: history.length, transition: surfaces.map(s => s.transition), wakeStations: wake.length,
      pressureIntegralDrag: coeffs.cdp }, warnings, limitations: coupledLimitations };
}
