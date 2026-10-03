import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';
import { selectSurfaceTransition, prepareSurfaceTransition, chooseMixedTransitionShear, minimizeMixedTransitionShear } from '../src/viscous/transition-selection.js';

const native = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/automatic-transition.json', import.meta.url)));
const close = (a, b, limit = 3e-10) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= limit * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b} (limit ${limit})`);
const kernel = c => createIntegralKernel({ ...c.parameters, exactJacobian: true });

test('automatic-transition fixtures bind executed original Fortran and its complete local harness', () => {
  assert.match(native.provenance.scope, /Executed unmodified original Fortran/);
  for (const [path, hash] of Object.entries(native.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${path}`, import.meta.url))).digest('hex'), hash, path);
});

test('free roots, earlier/later trips and mixed transition residuals match original Fortran', () => {
  for (const { input, expected } of native.cases) {
    const result = evaluateTransitionInterval(kernel(input), input);
    assert.equal(result.transition.forced, expected.selectedForced, input.name);
    close(result.transition.s, expected.s); close(result.transition.amplification, expected.amplification);
    result.residual.forEach((v, i) => close(v, expected.residual[i]));
    if (input.name.startsWith('natural') || input.name.startsWith('later')) assert.equal(result.transition.forced, false);
    else assert.equal(result.transition.forced, true);
  }
});

test('all local transition/root derivatives match independent fourth-order Fortran differences', t => {
  let worst = { error: 0 };
  const values = native.cases.map(({ input }) => evaluateTransitionInterval(kernel(input), input, { jacobian: true }));
  for (const c of native.derivatives) {
    const p = values[c.case].partials, col = ['aux', 'theta', 'deltaStar', 'ue', 's'].indexOf(c.key);
    const actual = c.side ? [...p[c.side].map(row => row[col]), p.location[c.side][col]] : [...p.trip, p.location.trip];
    actual.forEach((v, r) => {
      const error = Math.abs(v - c.expected[r]) / Math.max(1, Math.abs(v), Math.abs(c.expected[r]));
      if (error > worst.error) worst = { error, case: native.cases[c.case].input.name, side: c.side, key: c.key, row: r };
      close(v, c.expected[r], 2e-6);
    });
  }
  t.diagnostic(JSON.stringify({ columns: native.derivatives.length, worst }));
});

test('coordinate-chain derivatives include motion of natural transition, with translation invariance', () => {
  const c = native.cases.find(c => c.input.name === 'natural-M0.2').input, k = kernel(c);
  const b = evaluateTransitionInterval(k, c, { jacobian: true }), p = b.partials;
  // A model grid coordinate changes both station distances and BL primitives.
  const up = [.17, 2e-5, -3e-5, .08, .03], down = [.01, -1e-5, 4e-5, -.04, -.02], fields = ['aux', 'theta', 'deltaStar', 'ue', 's'];
  const change = h => Object.fromEntries(['upstream', 'downstream'].map((side, i) => [side,
    Object.fromEntries(fields.map((f, j) => [f, c[side][f] + h * (i ? down : up)[j]]))]));
  const h = 1e-5, a = evaluateTransitionInterval(k, { ...c, ...change(h) }), z = evaluateTransitionInterval(k, { ...c, ...change(-h) });
  for (let row = 0; row < 3; row++) {
    const exact = up.reduce((s, v, j) => s + p.upstream[row][j] * v + p.downstream[row][j] * down[j], 0);
    close(exact, (a.residual[row] - z.residual[row]) / (2 * h), 2e-6);
  }
  const motion = up.reduce((s, v, j) => s + p.location.upstream[j] * v + p.location.downstream[j] * down[j], 0);
  assert.ok(Math.abs(motion) > .01); close(motion, (a.transition.s - z.transition.s) / (2 * h), 2e-6);
  close(p.location.upstream[4] + p.location.downstream[4], 1, 2e-7);
  assert.deepEqual(p.trip, [0, 0, 0]); assert.equal(p.location.trip, 0);
});

test('derivatives stay on the selected branch immediately on either side of a natural/trip switch', () => {
  const c = native.cases.find(c => c.input.name === 'natural-M0.2').input, k = kernel(c);
  const free = evaluateTransitionInterval(k, c, { jacobian: true });
  const later = evaluateTransitionInterval(k, { ...c, tripS: free.transition.s + 1e-10 }, { jacobian: true });
  const earlier = evaluateTransitionInterval(k, { ...c, tripS: free.transition.s - 1e-10 }, { jacobian: true });
  assert.equal(later.transition.forced, false); assert.equal(earlier.transition.forced, true);
  for (const side of ['upstream', 'downstream']) {
    later.partials[side].flat().forEach((v, i) => close(v, free.partials[side].flat()[i], 2e-6));
    later.partials.location[side].forEach((v, i) => close(v, free.partials.location[side][i], 2e-6));
    earlier.partials.location[side].forEach(v => close(v, 0, 2e-7));
  }
  close(earlier.partials.location.trip, 1, 2e-7);
});

test('surface selection matches independently propagated Fortran for Ncrit and earlier-trip changes', t => {
  let locationError = 0, amplificationError = 0;
  for (const c of native.profiles) {
    const before = structuredClone(c.states), result = selectSurfaceTransition(kernel(c), c.states, c);
    assert.deepEqual(c.states, before);
    assert.equal(result.index, c.expected.index); assert.equal(result.kind, c.expected.kind);
    // Original TRCHEK stops its N root at a 5e-5 correction; the JS wrapper
    // additionally closes the same laminar N equation. Compare to that native
    // stopping accuracy rather than claiming bitwise identity of the wrappers.
    close(result.s, c.expected.s, 5e-5);
    locationError = Math.max(locationError, Math.abs(result.s - c.expected.s));
    result.amplification.forEach((v, i) => {
      close(v, c.expected.amplification[i], 5e-5);
      amplificationError = Math.max(amplificationError, Math.abs(v - c.expected.amplification[i]));
    });
  }
  t.diagnostic(JSON.stringify({ profiles: native.profiles.length, locationError, amplificationError }));
});

test('laminar surfaces and trips before the first station are distinct transition states', () => {
  const states = native.profiles[0].states.map(s => ({ ...s, deltaStar: 2.3 * s.theta }));
  const k = createIntegralKernel({ reynolds: 1e4, mach: .2, ncrit: 9 });
  const laminar = selectSurfaceTransition(k, states);
  assert.equal(laminar.index, null); assert.equal(laminar.kind, 'laminar'); assert.equal(laminar.s, null);
  assert.ok(laminar.amplification.every(n => n >= 0 && n < k.parameters.ncrit));
  const early = selectSurfaceTransition(k, states, { tripS: .002 });
  assert.equal(early.index, 0); assert.equal(early.kind, 'forced'); assert.equal(early.s, .002);
  assert.throws(() => selectSurfaceTransition(k, states, { tripS: -1 }), /controls/);
  assert.throws(() => selectSurfaceTransition(k, [states[1], states[0]]), /increase/);
});

test('entering terminal transition initializes the exact original-Fortran shear equation without changing thickness rows', () => {
  const terminal = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/terminal-transition.json', import.meta.url)));
  for (const c of terminal.cases) {
    const k = kernel(c), states = [c.upstream, { ...c.downstream, aux: 12 * c.expected.shear }], before = structuredClone(states);
    const p = prepareSurfaceTransition(k, states, { previousIndex: 0, tripS: c.tripS });
    assert.equal(p.kind, 'forced'); assert.equal(p.index, 1); assert.equal(p.s, c.downstream.s);
    close(p.auxiliary[1], c.expected.shear, 3e-12);
    assert.equal(p.converted.at(-1).terminalShearInitialization, true);
    const value = evaluateTransitionInterval(k, { upstream: { ...states[0], aux: p.auxiliary[0] },
      downstream: { ...states[1], aux: p.auxiliary[1] }, tripS: c.tripS });
    close(value.residual[0], 0, 1e-12);
    close(value.residual[1], c.expected.transition[1]); close(value.residual[2], c.expected.transition[2]);
    assert.deepEqual(states, before);
  }
});

test('unchanged terminal shear is not reset and interior trips use their mixed equation', () => {
  const c = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/terminal-transition.json', import.meta.url))).cases[0];
  const k = kernel(c), states = [c.upstream, { ...c.downstream, aux: .3 }];
  const unchanged = prepareSurfaceTransition(k, states, { previousIndex: 1, tripS: c.tripS });
  assert.equal(unchanged.changed, false); assert.equal(unchanged.auxiliary[1], .3); assert.deepEqual(unchanged.converted, []);
  const interior = prepareSurfaceTransition(k, states, { previousIndex: 0, tripS: .9 });
  assert.equal(interior.kind, 'forced'); assert.equal(interior.s, .9);
  const value = evaluateTransitionInterval(k, { upstream: { ...states[0], aux: interior.auxiliary[0] },
    downstream: { ...states[1], aux: interior.auxiliary[1] }, tripS: .9 });
  close(value.residual[0], 0, 1e-12);
  assert.ok(interior.converted.every(p => !p.terminalShearInitialization));
});

test('crossings in both directions prepare N/shear conversions without changing thickness or velocity', () => {
  const c = native.profiles.find(c => c.parameters.ncrit === 9 && !c.tripS), k = kernel(c);
  const base = selectSurfaceTransition(k, c.states), states = structuredClone(c.states);
  for (let j = 0; j < base.index; j++) states[j].aux = base.amplification[j];
  const before = structuredClone(states);
  const upstream = prepareSurfaceTransition(k, states, { previousIndex: base.index, tripS: .18 });
  assert.equal(upstream.changed, true); assert.ok(upstream.index < base.index);
  upstream.converted.forEach(({ index, aux, to }) => {
    assert.equal(to, 'turbulent');
    if (index === upstream.index) close(evaluateTransitionInterval(k, {
      upstream: { ...states[index - 1], aux: upstream.auxiliary[index - 1] },
      downstream: { ...states[index], aux }, tripS: .18 }).residual[0], 0, 1e-12);
    else close(aux, c.checks[index - 1].shear);
  });
  const moved = states.map((s, j) => ({ ...s, aux: upstream.auxiliary[j] }));
  const downstream = prepareSurfaceTransition(k, moved, { previousIndex: upstream.index });
  assert.equal(downstream.index, base.index); assert.equal(downstream.changed, true);
  downstream.converted.forEach(({ index, aux, to }) => {
    if (index === downstream.index) { assert.equal(to, 'transition'); return; }
    assert.equal(to, 'laminar'); close(aux, base.amplification[index]);
  });
  assert.deepEqual(states, before);
  const noChange = prepareSurfaceTransition(k, states, { previousIndex: base.index });
  assert.equal(noChange.changed, false); assert.deepEqual(noChange.auxiliary, states.map(s => s.aux));
});

test('a new station distribution can initialize N even when its selected transition interval is unchanged', () => {
  const c = native.profiles.find(c => c.parameters.ncrit === 9 && !c.tripS), k = kernel(c);
  const selected = selectSurfaceTransition(k, c.states);
  const states = c.states.map((s, j) => ({ ...s, aux: j < selected.index ? 0 : .04 }));
  const before = structuredClone(states), controls = { previousIndex: selected.index };
  const unchanged = prepareSurfaceTransition(k, states, controls);
  assert.deepEqual(unchanged.auxiliary, states.map(s => s.aux));
  const j = selected.index;
  assert.throws(() => evaluateTransitionInterval(k, { upstream: states[j - 1], downstream: states[j] }), /outside this active interval/);
  const mapped = prepareSurfaceTransition(k, states, { ...controls, reinitializeAmplification: true });
  assert.equal(mapped.index, selected.index); assert.equal(mapped.changed, false);
  assert.deepEqual(mapped.converted, []); assert.deepEqual(states, before);
  for (let i = 0; i < states.length; i++) assert.equal(mapped.auxiliary[i], i < j ? selected.amplification[i] : states[i].aux);
  const value = evaluateTransitionInterval(k, { upstream: { ...states[j - 1], aux: mapped.auxiliary[j - 1] }, downstream: states[j] });
  assert.equal(value.transition.forced, false); close(value.transition.s, selected.s);
});

test('natural transition moves across stations with Ncrit and restores laminar amplification on return', () => {
  const high = native.profiles.find(c => c.parameters.ncrit === 9 && !c.tripS);
  const low = native.profiles.find(c => c.parameters.ncrit === 6 && !c.tripS);
  const states = high.states.map((s, j) => ({ ...s, aux: j < high.expected.index ? high.expected.amplification[j] : .03 }));
  const toLow = prepareSurfaceTransition(kernel(low), states, { previousIndex: high.expected.index });
  assert.equal(toLow.kind, 'natural'); assert.equal(toLow.index, low.expected.index);
  assert.ok(toLow.converted.length >= 2);
  toLow.converted.forEach(({ index, aux, to }) => {
    assert.equal(to, 'turbulent'); if (index !== toLow.index) close(aux, high.checks[index - 1].shear);
  });
  const toHigh = prepareSurfaceTransition(kernel(high), states.map((s, j) => ({ ...s, aux: toLow.auxiliary[j] })), { previousIndex: toLow.index });
  assert.equal(toHigh.kind, 'natural'); assert.equal(toHigh.index, high.expected.index);
  toHigh.converted.forEach(({ index, aux, to }) => {
    if (index === toHigh.index) { assert.equal(to, 'transition'); return; }
    assert.equal(to, 'laminar'); close(aux, high.expected.amplification[index], 5e-5);
  });
});

test('RAE downstream transition selects the shear preserving turbulent transport', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/rae64x11-downstream-transition-profile.json', import.meta.url)));
  const before = structuredClone(f), k = createIntegralKernel(f.parameters);
  const p = prepareSurfaceTransition(k, f.states, { ...f.controls, previousIndex: f.previousIndex });
  assert.equal(p.kind, 'natural'); assert.equal(p.index, 31); assert.equal(p.changed, true);
  const i = p.index;
  const selected = chooseMixedTransitionShear(k, {
    upstream: { ...f.states[i - 1], aux: p.auxiliary[i - 1] }, downstream: f.states[i], following: f.states[i + 1],
    tripS: f.controls.tripS, initializedAux: p.auxiliary[i], shearWeight: 20,
  });
  assert.equal(selected.diagnostics.selected, 'existing');
  assert.ok(selected.diagnostics.existing.squaredResidual < selected.diagnostics.initialized.squaredResidual);
  const tailResidual = states => k.interval({ upstream: states[i], downstream: states[i + 1], regime: 'turbulent' }).residual;
  const oldResidual = tailResidual(f.states);
  const mapped = f.states.map((s, j) => ({ ...s, aux: j === i ? selected.aux : p.auxiliary[j] }));
  // The old scalar shear initializer changed this real downstream interval's
  // shape residual from 0.000264 to 0.00956 without changing thickness or Ue.
  assert.deepEqual(tailResidual(mapped), oldResidual);
  assert.equal(selected.aux, f.states[i].aux);
  // The opposite crossing still converts a laminar N into a shear guess.
  const upstream = prepareSurfaceTransition(k, f.states, { ...f.controls, previousIndex: 32 });
  assert.equal(upstream.index, 31);
  assert.equal(upstream.converted.find(c => c.index === 31).transitionShearInitialization, true);
  assert.deepEqual(f, before);
});

test('selection on four surfaces is independent and failed preparation cannot partially mutate callers', () => {
  const c = native.profiles[0], states = structuredClone(c.states), before = structuredClone(states);
  const k = kernel(c), expected = selectSurfaceTransition(k, states);
  const choices = [undefined, .015, .2, .4];
  for (const tripS of choices) selectSurfaceTransition(k, states, { tripS });
  assert.deepEqual(selectSurfaceTransition(k, states), expected);
  const station = k.station; let calls = 0;
  k.station = (s, regime) => {
    if (regime === 'turbulent' && ++calls === 2) throw new Error('deliberate invalid second conversion');
    return station(s, regime);
  };
  assert.throws(() => prepareSurfaceTransition(k, states, { previousIndex: states.length - 1, tripS: .015 }), /second conversion/);
  assert.equal(calls, 2); assert.deepEqual(states, before);
});

test('RAE upstream transition initializes shear using both affected native intervals', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/rae32x9-upstream-transition-intervals.json', import.meta.url)));
  const before = structuredClone(f), k = createIntegralKernel(f.parameters);
  const score = aux => {
    const endpoint = { ...f.downstream, aux };
    const a = evaluateTransitionInterval(k, { upstream: f.upstream, downstream: endpoint, tripS: f.tripS });
    const b = k.interval({ upstream: endpoint, downstream: f.following, regime: 'turbulent' });
    return { transition: a.transition, residual: a.residual,
      squared: [a.residual, b.residual].reduce((sum, row) => sum + row.reduce((s, v, i) => s + (v * (i === 0 ? f.shearWeight : 1)) ** 2, 0), 0) };
  };
  const initial = score(f.initializedAux), selected = minimizeMixedTransitionShear(k, f), result = score(selected.aux);
  assert.ok(Math.abs(initial.residual[0]) < 1e-12, 'The old initializer closes its shear row.');
  assert.ok(result.squared < .5 * initial.squared, 'Transport-aware initialization must remove the captured residual jump.');
  assert.deepEqual(result.transition, initial.transition, 'Only the shear guess changes, not the natural transition root.');
  assert.ok(selected.aux > 0 && selected.aux <= .25);
  assert.equal(selected.diagnostics.selected, 'minimum');
  assert.equal(selected.diagnostics.equationsChanged, false);
  for (const factor of [.99, 1.01]) assert.ok(result.squared <= score(selected.aux * factor).squared,
    'The selected shear is a local minimum of all six equations.');
  assert.deepEqual(f, before);
  const capped = minimizeMixedTransitionShear(k, { ...f, initializedAux: 1 });
  assert.ok(capped.aux > 0 && capped.aux <= .25, 'An oversized scalar root must enter the bounded shear domain.');
  assert.equal(capped.diagnostics.selected, 'minimum');
  assert.ok(score(capped.aux).squared <= score(.25).squared);
  assert.throws(() => minimizeMixedTransitionShear(k, { ...f, shearWeight: 0 }), /weight/);
  assert.throws(() => minimizeMixedTransitionShear(k, { ...f, initializedAux: NaN }), /initial transition shear/);
});
