// SPDX-License-Identifier: GPL-2.0-or-later
// A forced trip before the first resolved BL station. The virtual upstream
// state follows the local linear stagnation velocity and the same similarity
// equations. Its dependence on downstream Ue/s and trip distance is included
// in the composite local derivatives. Resolution must still be checked.
import { solveNewton } from '../numerics/newton.js';
import { boundedResidualPartial } from '../numerics/bounded-partial.js';
import { evaluateTransitionInterval } from './transition-interval.js';

export function evaluateLeadingTransitionInterval(kernel, { downstream, tripS }, { jacobian = false } = {}) {
  const scale = 1 / Math.sqrt(kernel.parameters.reynolds);
  const evaluate = (state, trip) => {
    if (!(trip > 0 && trip <= state.s)) throw new Error('Leading transition must precede the first BL station.');
    const s = .5 * trip, slope = state.ue / state.s, ue = slope * s;
    const at = z => ({ s, ue, aux: 0, theta: scale * z[0], deltaStar: scale * z[1] });
    const root = solveNewton({ initial: [.292 / Math.sqrt(slope), .647 / Math.sqrt(slope)],
      tolerance: 1e-12, maxIterations: 15, admissible: z => z[0] > 0 && z[1] > z[0],
      residual: z => kernel.interval({ downstream: at(z), regime: 'similarity' }).residual.slice(1),
      jacobian: z => Float64Array.from(kernel.interval({ downstream: at(z), regime: 'similarity' }).downstream.slice(1).flatMap(r => [scale * r[1], scale * r[2]])) });
    if (!root.converged) throw new Error('Leading transition similarity state did not converge.');
    const upstream = at(root.x);
    return { ...evaluateTransitionInterval(kernel, { upstream, downstream: state, tripS: trip }), upstreamState: upstream };
  };
  const value = evaluate(downstream, tripS);
  if (!jacobian) return value;
  const packed = r => [...r.residual, r.transition.s], base = packed(value), columns = [];
  for (const [k, key] of ['aux', 'theta', 'deltaStar', 'ue', 's', 'tripS'].entries()) {
    const v = key === 'tripS' ? tripS : downstream[key];
    const h = Math.cbrt(Number.EPSILON) * Math.max(Math.abs(v), k === 0 ? .01 : k < 3 ? 1e-7 : 1e-12);
    const f = w => {
      const r = key === 'tripS' ? evaluate(downstream, w) : evaluate({ ...downstream, [key]: w }, tripS);
      if (r.transition.forced !== value.transition.forced) throw new Error('Leading transition derivative crossed its active branch.');
      return packed(r);
    };
    columns.push(boundedResidualPartial(f, v, { step: h, base,
      lower: key === 's' ? tripS : key === 'tripS' ? 0 : -Infinity,
      upper: key === 'tripS' ? downstream.s : Infinity }));
  }
  return { ...value, partials: {
    upstream: Array.from({ length: 3 }, () => [0, 0, 0, 0, 0]),
    downstream: Array.from({ length: 3 }, (_, r) => columns.slice(0, 5).map(c => c[r])), trip: columns[5].slice(0, 3),
    location: { upstream: [0, 0, 0, 0, 0], downstream: columns.slice(0, 5).map(c => c[3]), trip: columns[5][3] },
  } };
}
