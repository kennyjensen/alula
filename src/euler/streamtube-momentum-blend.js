// SPDX-License-Identifier: GPL-2.0-or-later
// Research ISMOM4-style local momentum/entropy blend. MSES 3.05 manual
// §§1.2.4–1.2.5 specifies deltaq*Delta(rho*q/p), a compression preference,
// and continuous blending, but does NOT publish the complete switch.
// Our choices are a C2 quintic, an OR of the two endpoint loss indicators,
// and a log-speed compression scale sqrt(epsilonP). epsilonP is mandatory;
// there is no inferred MSES tolerance or hidden default.
//
// Pass physical section states and their already shared transport speeds.
// R1/R2 must have identical units: normally raw cell pressure residuals,
// with R2 = -mean(p)*Delta(S_tilde), as evaluated by streamtube-cell.js.
// This helper changes only the streamwise equation, not the physical gas.
// A hybrid root is not a conservative momentum root. The local indicators
// can miss weak shocks and affine portions of a second-order shock profile.

function inputs(parameters) {
  const { states, transportSpeeds, streamwiseResidual, isentropicResidual, epsilonP } = parameters ?? {};
  if (!Array.isArray(states) || states.length !== 2 || !states.every(s => s &&
    ['rho', 'q', 'p'].every(k => Number.isFinite(s[k]) && s[k] > 0))
    || !Array.isArray(transportSpeeds) || transportSpeeds.length !== 2
    || !transportSpeeds.every(q => Number.isFinite(q) && q >= 0)
    || ![streamwiseResidual, isentropicResidual, epsilonP].every(Number.isFinite) || epsilonP <= 0)
    throw new Error('Invalid streamtube momentum-blend inputs; epsilonP must be explicit and positive.');
  return { states: states.map(({ rho, q, p }) => ({ rho, q, p })), transportSpeeds: [...transportSpeeds],
    streamwiseResidual, isentropicResidual, epsilonP };
}

function step(t) {
  if (t <= 0) return { value: 0, derivative: 0 };
  if (t >= 1) return { value: 1, derivative: 0 };
  // Evaluate symmetrically to avoid polynomial overshoot near one.
  const x = t <= .5 ? t : 1 - t;
  const y = x * x * x * (10 + x * (-15 + 6 * x));
  return { value: t <= .5 ? y : 1 - y, derivative: 30 * t * t * (1 - t) * (1 - t) };
}

function evaluate(p) {
  const phi = p.states.map(s => s.rho * (s.q / s.p));
  const deltaPhi = phi[1] - phi[0];
  const bias = p.states.map((s, i) => p.transportSpeeds[i] - s.q);
  const lossIndicators = bias.map(q => -q * deltaPhi);
  // Relative differences retain near-equal speeds; logs handle wide ratios.
  const relative = (p.states[0].q - p.states[1].q) / p.states[1].q;
  const compression = Math.abs(relative) < .5 ? Math.log1p(relative)
    : Math.log(p.states[0].q) - Math.log(p.states[1].q);
  if (![...phi, deltaPhi, ...bias, ...lossIndicators, compression].every(Number.isFinite))
    throw new Error('Nonfinite streamtube momentum-blend indicator.');
  const losses = lossIndicators.map(l => step(l / p.epsilonP));
  const compressionScale = Math.sqrt(p.epsilonP), gate = step(compression / compressionScale);
  const either = losses[0].value + (1 - losses[0].value) * losses[1].value;
  const fraction = gate.value * either;
  const residual = fraction === 0 ? p.isentropicResidual : fraction === 1 ? p.streamwiseResidual
    : fraction * p.streamwiseResidual + (1 - fraction) * p.isentropicResidual;
  if (!Number.isFinite(residual)) throw new Error('Nonfinite streamtube momentum-blend residual.');
  return { value: { fraction, residual, lossIndicators, compression },
    phi, deltaPhi, bias, losses, compressionScale, gate, either };
}

export function evaluateStreamtubeMomentumBlend(parameters) {
  return evaluate(inputs(parameters)).value;
}

export function linearizeStreamtubeMomentumBlend(parameters) {
  const p = inputs(parameters), prepared = evaluate(p);
  const { value, phi, deltaPhi, bias, losses, compressionScale, gate, either } = prepared;
  // The private numeric snapshot is independent of callers changing params
  // or the public value returned alongside apply().
  const fraction = value.fraction, compression = value.compression;
  const lossIndicators = [...value.lossIndicators];
  const apply = ({ states = [{}, {}], transportSpeeds = [0, 0],
    streamwiseResidual = 0, isentropicResidual = 0, epsilonP = 0 } = {}) => {
    if (!Array.isArray(states) || states.length !== 2 || !states.every(s => s &&
      ['rho', 'q', 'p'].every(k => s[k] === undefined || Number.isFinite(s[k])))
      || !Array.isArray(transportSpeeds) || transportSpeeds.length !== 2 || !transportSpeeds.every(Number.isFinite)
      || ![streamwiseResidual, isentropicResidual, epsilonP].every(Number.isFinite))
      throw new Error('Invalid streamtube momentum-blend tangent.');
    const state = states.map(s => ({ rho: s.rho ?? 0, q: s.q ?? 0, p: s.p ?? 0 }));
    const dPhi = p.states.map((s, i) => phi[i] *
      (state[i].rho / s.rho + state[i].q / s.q - state[i].p / s.p));
    const dDeltaPhi = dPhi[1] - dPhi[0];
    const dLoss = bias.map((q, i) => -(transportSpeeds[i] - state[i].q) * deltaPhi - q * dDeltaPhi);
    const dCompression = state[0].q / p.states[0].q - state[1].q / p.states[1].q;
    // Do not multiply zero plateau derivatives by an overflowing ratio.
    const dSteps = losses.map((l, i) => l.derivative === 0 ? 0 : l.derivative *
      (dLoss[i] - lossIndicators[i] / p.epsilonP * epsilonP) / p.epsilonP);
    const dGate = gate.derivative === 0 ? 0 : gate.derivative *
      (dCompression - .5 * compression * (epsilonP / p.epsilonP)) / compressionScale;
    let dFraction = dGate * either + gate.value *
      ((1 - losses[1].value) * dSteps[0] + (1 - losses[0].value) * dSteps[1]);
    // Preserve the selected equation exactly, including its Jacobian.
    if (fraction === 0 || fraction === 1) dFraction = 0;
    const residual = fraction === 0 ? isentropicResidual : fraction === 1 ? streamwiseResidual
      : fraction * streamwiseResidual + (1 - fraction) * isentropicResidual
        + dFraction * (p.streamwiseResidual - p.isentropicResidual);
    if (![dFraction, residual, ...dLoss, dCompression].every(Number.isFinite))
      throw new Error('Nonfinite streamtube momentum-blend derivative.');
    return { fraction: dFraction, residual, lossIndicators: dLoss, compression: dCompression };
  };
  return { value, apply };
}
