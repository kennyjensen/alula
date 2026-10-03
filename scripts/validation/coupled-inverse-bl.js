// SPDX-License-Identifier: GPL-2.0-or-later
// Research conditional inverse BL solution. All displaced-wall thicknesses
// and Euler variables are retained. No production solver imports this file.
import { solveNewton } from '../../src/numerics/newton.js';
import { evaluateTransitionInterval, checkAutomaticTransition } from '../../src/viscous/transition-interval.js';
import { evaluateLeadingTransitionInterval } from '../../src/viscous/leading-transition-interval.js';

export function prepareCoupledInverseBL(system, state, { tolerance = 2e-12, maxIterations = 30 } = {}) {
  const { bl, ne } = system, kernel = bl.kernel;
  if (bl.transitionMode !== 'automatic' || state.length !== system.n || !state.every(Number.isFinite))
    throw new Error('Invalid conditional inverse BL input.');
  const geometry = bl.geometry(state.subarray(0, ne)), x = state.slice();
  const states = bl.stations.map(({ id }) => ({ s: geometry.coordinates[id].s, aux: x[ne + 4 * id],
    theta: bl.scale * x[ne + 4 * id + 1], deltaStar: bl.scale * x[ne + 4 * id + 2], ue: x[ne + 4 * id + 3] }));
  const intervals = [], phase = [], summary = { updates: 0, maximumLocalResidual: 0, surfaces: bl.surfaces.length, wakes: bl.wakes.length };
  const put = (id, p) => {
    states[id] = p; x[ne + 4 * id] = p.aux; x[ne + 4 * id + 1] = p.theta / bl.scale;
    x[ne + 4 * id + 3] = p.ue;
  };
  const solve = (id, previous, type, tripS) => {
    const original = states[id], upstream = states[previous], similarity = type === 'similarity';
    const guess = { ...original };
    if (similarity) guess.aux = 0;
    else if (type === 'laminar-or-transition' && ['laminar', 'similarity'].includes(bl.stations[id].regime))
      guess.aux = kernel.station({ ...guess, aux: .03 }, 'turbulent').transitionShear;
    if (!similarity && !(guess.aux > 0)) guess.aux = .03;
    const keys = similarity ? ['theta', 'ue'] : ['aux', 'theta', 'ue'];
    const columns = similarity ? [1, 3] : [0, 1, 3], rows = similarity ? [1, 2] : [0, 1, 2];
    const at = z => ({ ...guess, ...Object.fromEntries(keys.map((key, i) => [key, Math.exp(z[i])])) });
    const evaluate = (p, jacobian = false) => {
      if (type === 'leading-transition') return evaluateLeadingTransitionInterval(kernel, { downstream: p, tripS }, { jacobian });
      if (type === 'laminar-or-transition') {
        // When natural/forced transition is beyond this station, force only
        // this local block at its downstream endpoint. Its turbulent portion
        // has zero length, so momentum/shape remain the laminar equations.
        // This keeps Ctau's meaning continuous inside the local Newton solve;
        // the artificial endpoint trip is never committed to the surface.
        return evaluateTransitionInterval(kernel, { upstream, downstream: p, tripS: Math.min(tripS, p.s) }, { jacobian });
      }
      return kernel.interval({ upstream, downstream: p, regime: type });
    };
    const scaled = (value, row) => value * (row === 0 ? 20 : 1);
    const result = solveNewton({ initial: keys.map(k => Math.log(guess[k])), tolerance, maxIterations,
      admissible: z => {
        const p = at(z);
        if (!z.every(Number.isFinite) || !(p.theta < p.deltaStar)) return false;
        try { return kernel.station(p, similarity ? 'laminar' : type === 'wake' ? 'wake' : 'turbulent').machSquared < 1; }
        catch { return false; }
      },
      residual: z => { const v = evaluate(at(z)); return rows.map(r => scaled(v.residual[r], r)); },
      jacobian: z => {
        const p = at(z), v = evaluate(p, true), d = v.partials?.downstream ?? v.downstream;
        return Float64Array.from(rows.flatMap(r => columns.map((col, j) => scaled(d[r][col] * p[keys[j]], r))));
      } });
    const p = at(result.x);
    if (!result.converged) throw Object.assign(new Error(`Conditional inverse BL failed at ${id} (${type}): ${result.reason}`),
      { diagnostic: { id, previous, type, original, upstream, final: p, history: result.history, completed: intervals.length } });
    summary.updates += result.history.length - 1;
    summary.maximumLocalResidual = Math.max(summary.maximumLocalResidual, result.history.at(-1).residual);
    put(id, p);
    const record = { id, previous, type, history: result.history, input: { upstream, downstream: p, regime: type, tripS } };
    intervals.push(record);
    return record;
  };
  for (const [k, surface] of bl.surfaces.entries()) {
    const tripS = geometry.surfaceData[k].tripS; let turbulent = false;
    for (const [j, id] of surface.ids.entries()) {
      const previous = j ? surface.ids[j - 1] : id;
      const type = turbulent ? 'turbulent' : j === 0 ? (tripS <= states[id].s ? 'leading-transition' : 'similarity') : 'laminar-or-transition';
      const record = solve(id, previous, type, tripS);
      if (type === 'similarity' || type === 'turbulent') continue;
      const check = type === 'leading-transition' ? { transition: true } : checkAutomaticTransition(kernel,
        { upstream: states[previous], downstream: states[id], tripS });
      if (check.transition) { phase[k] = j; turbulent = true; record.input.regime = type === 'leading-transition' ? type : 'transition'; }
      else {
        put(id, { ...states[id], aux: check.amplification });
        record.input = { ...record.input, downstream: states[id], regime: 'laminar' };
      }
    }
    if (!turbulent) throw new Error('Missing terminal transition in inverse BL preparation.');
  }
  for (const w of bl.wakes) {
    const [upper, lower] = bl.surfaces.filter(s => s.body === w.body).map(s => states[s.ids.at(-1)]), id = w.ids[0];
    const theta = upper.theta + lower.theta, deltaStar = upper.deltaStar + lower.deltaStar;
    const matched = { theta, deltaStar, aux: (upper.theta * upper.aux + lower.theta * lower.aux) / theta };
    // This TE matching station has no displacement-map entry. Its thickness
    // is a sum, while its Ue remains a retained outer-coupling variable.
    if (bl.thicknessMap.some(row => row.has(ne + 4 * id + 2))) throw new Error('TE matching would move the displaced grid.');
    put(id, { ...states[id], ...matched }); x[ne + 4 * id + 2] = deltaStar / bl.scale;
    intervals.push({ id, input: { regime: 'te', upper, lower, downstream: states[id], matched } });
    for (const id of w.ids.slice(1)) solve(id, id - 1, 'wake');
  }
  return { x, phase, intervals, summary };
}
