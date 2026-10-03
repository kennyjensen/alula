import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { relaxXfoilProfiles } from '../src/viscous/xfoil-profile-relaxation.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-body.json', import.meta.url)));
const original = fixture.cases.find(c => c.name.startsWith('default-two-element-root'));

test('current-profile MRCHDU retains original Fortran parity including transition changes', t => {
  for (const [path, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(path)).digest('hex'), hash);
  let stations = 0, changedIntervals = 0, maxError = 0;
  for (const c of fixture.cases) {
    const input = structuredClone(c.input), result = relaxXfoilProfiles(input, c.parameters);
    assert.deepEqual(input, c.input); assert.equal(result.flowSolved, false);
    assert.deepEqual(result.localConvergenceWarnings, []);
    result.surfaces.forEach((a, side) => {
      assert.equal(a.transition, c.native.surfaces[side].transition);
      assert.equal(a.forced, c.native.surfaces[side].forced);
      assert.ok(Math.abs(a.s - c.native.surfaces[side].s) < 1e-12);
      changedIntervals += a.transition !== input.phases[side] ? 1 : 0;
    });
    const before = [...input.surfaces, input.wake], after = [...result.surfaces.map(s => s.states), result.wake];
    const native = [...c.native.surfaces.map(s => s.states), c.native.wake], changes = [...result.changes.surfaces, result.changes.wake];
    after.forEach((p, side) => p.forEach((a, i) => {
      assert.equal(a.s, before[side][i].s);
      for (const k of ['ue', 'theta', 'deltaStar', 'aux']) {
        const expected = native[side][i][k], error = Math.abs(a[k] - expected) / (k === 'aux' ? Math.max(.01, Math.abs(expected)) : Math.abs(expected));
        maxError = Math.max(maxError, error); assert.ok(error < 1e-10, `${c.name}/${side}/${i}/${k}: ${error}`);
        assert.equal(changes[side][i][k], a[k] - before[side][i][k]);
      }
      stations++;
    }));
  }
  assert.equal(stations, 2590); assert.equal(changedIntervals, 2);
  t.diagnostic(JSON.stringify({ stations, changedIntervals, maxError }));
});

test('finite gap is preserved in distance and merged total displacement, subtracted once by BL closure', t => {
  const input = structuredClone(original.input), normalGap = .1 * Math.min(...input.surfaces.map(p => p.at(-1).theta));
  const gap = createXfoilDeadAirGap({ normalGap, upperDerivative: { x: -1, y: .1 }, lowerDerivative: { x: 1, y: .1 } });
  input.normalGap = normalGap; input.tripS = input.surfaces.map(p => p.at(-1).s);
  input.wake.forEach(p => { p.wakeGap = gap.at(p.s - input.wake[0].s).gap; p.deltaStar += p.wakeGap; });
  const saved = structuredClone(input), result = relaxXfoilProfiles(input, original.parameters);
  assert.deepEqual(input, saved); assert.deepEqual(result.localConvergenceWarnings, []);
  const [upper, lower] = result.surfaces.map(p => p.states.at(-1)), first = result.wake[0];
  const kernel = createIntegralKernel(original.parameters), matching = kernel.trailingEdge(upper, lower, first, normalGap);
  const residual = matching.residual ?? matching;
  assert.ok(Math.max(...residual.map(Math.abs)) < 1e-12);
  assert.ok(Math.abs(first.deltaStar - upper.deltaStar - lower.deltaStar - normalGap) < 1e-12);
  result.wake.forEach((p, i) => {
    assert.equal(p.s, saved.wake[i].s); assert.equal(p.wakeGap, saved.wake[i].wakeGap);
    const properties = kernel.station(p, 'wake');
    assert.equal(properties.h, (p.deltaStar - p.wakeGap) / p.theta);
  });
  assert.equal(result.normalGap, normalGap); assert.deepEqual(result.tripS, saved.tripS);
  t.diagnostic(JSON.stringify({ normalGap, matchingResidual: Array.from(residual), wakeStations: result.wake.length }));
});

test('MRCHDU rejects unsupported trips and inconsistent phase/gap inputs without changing caller state', () => {
  const cases = [
    [a => { a.tripS = a.surfaces.map(p => .5 * p.at(-1).s); }, /terminal trips/],
    [a => { a.phases[0] = 0; }, /previous transition/],
    [a => { a.surfaces[0][a.phases[0]].aux = 0; }, /auxiliary state/],
    [a => { a.surfaces[0][a.phases[0]].aux = -1e-6; }, /auxiliary state/],
    [a => { a.surfaces[0][0].aux = NaN; }, /complete ordered/],
    [a => { a.surfaces[0][0].aux = -Infinity; }, /complete ordered/],
    [a => { a.wake[0].aux = 0; }, /Invalid MRCHDU wake state/],
    [a => { a.normalGap = 1e-4; }, /first wake gap/],
    [a => { a.wake[0].s += 1e-7; }, /first wake station/],
    [a => { a.surfaces[0][0].wakeGap = 1e-12; }, /solid surface/],
  ];
  for (const [change, pattern] of cases) {
    const input = structuredClone(original.input); change(input); const saved = structuredClone(input);
    assert.throws(() => relaxXfoilProfiles(input, original.parameters), pattern); assert.deepEqual(input, saved);
  }
  assert.throws(() => relaxXfoilProfiles(original.input, { ...original.parameters, velocityConvention: 'xfoil' }), /physical edge/);
});

test('finite negative old-laminar N remarches like original Fortran with strict returned states and unchanged caller input', t => {
  const reference = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-negative-laminar-input.json', import.meta.url)));
  const hash = p => createHash('sha256').update(fs.readFileSync(p)).digest('hex');
  for (const [p, expected] of Object.entries(reference.provenance.sha256)) assert.equal(hash(p), expected, p);
  for (const name of ['report', 'nativeOutput', 'nativeInput'])
    assert.equal(hash(reference.provenance[name]), reference.provenance[`${name}Hash`]);
  let comparisons = 0, maximumRelativeError = 0, negativeInputs = 0;
  for (const c of reference.cases) {
    assert.equal(hash(c.source.path), c.source.hash);
    if (c.source.checkpoint) assert.equal(hash(c.source.checkpoint), c.source.checkpointHash);
    const input = structuredClone(c.input), before = structuredClone(input), result = relaxXfoilProfiles(input, c.parameters);
    assert.deepEqual(input, before); assert.equal(result.flowSolved, false);
    assert.deepEqual(result.localConvergenceWarnings, []); assert.deepEqual(c.native.warnings, []);
    const kernel = createIntegralKernel(c.parameters);
    input.surfaces.forEach((profile, side) => profile.forEach((p, i) => {
      if (i < input.phases[side] && p.aux < 0) negativeInputs++;
    }));
    result.surfaces.forEach((s, side) => {
      assert.equal(s.transition, c.native.surfaces[side].transition);
      assert.equal(s.forced, c.native.surfaces[side].forced);
      assert.ok(Math.abs(s.s - c.native.surfaces[side].s) < 1e-12);
    });
    [...result.surfaces.map(s => s.states), result.wake].forEach((profile, part) => {
      const expected = [...c.native.surfaces.map(s => s.states), c.native.wake][part];
      assert.equal(profile.length, expected.length);
      profile.forEach((p, i) => {
        assert.equal(p.s, expected[i].s);
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const error = Math.abs(p[key] - expected[i][key])
            / (key === 'aux' ? Math.max(.01, Math.abs(expected[i][key])) : Math.abs(expected[i][key]));
          assert.ok(error < 1e-10, `${c.name}/${part}/${i}/${key}: ${error}`);
          comparisons++; maximumRelativeError = Math.max(maximumRelativeError, error);
        }
        const laminar = part < 2 && i < result.surfaces[part].transition;
        assert.ok(laminar ? p.aux >= 0 && p.aux < c.parameters.ncrit : p.aux > 0);
        kernel.station(p, part === 2 ? 'wake' : laminar ? 'laminar' : 'turbulent');
      });
    });
  }
  assert.equal(negativeInputs, 26); assert.equal(comparisons, 1048);
  t.diagnostic(JSON.stringify({ negativeInputs, comparisons, maximumRelativeError }));
});

test('old-laminar N above Ncrit is remarched with the previous phase, matching original Fortran', t => {
  const reference = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-old-phase-input.json', import.meta.url)));
  const hash = p => createHash('sha256').update(fs.readFileSync(p)).digest('hex');
  for (const [p, expected] of Object.entries(reference.provenance.sha256)) assert.equal(hash(p), expected, p);
  assert.equal(hash(reference.parent.path), reference.parent.hash);
  assert.equal(hash(reference.provenance.nativeOutput), reference.provenance.nativeOutputHash);
  assert.equal(hash(reference.provenance.report), reference.provenance.reportHash);
  const input = structuredClone(fixture.cases.find(c => c.name === reference.parent.name).input);
  const { side, index, oldPhase, before, after, ncrit } = reference.mutation;
  assert.equal(input.phases[side], oldPhase); assert.equal(index, oldPhase - 1);
  assert.equal(input.surfaces[side][index].aux, before); assert.ok(after > ncrit);
  input.surfaces[side][index].aux = after;
  const saved = structuredClone(input), actual = relaxXfoilProfiles(input, reference.parameters);
  assert.deepEqual(input, saved); assert.deepEqual(input.phases, [22, 48]);
  assert.deepEqual(actual.localConvergenceWarnings, []);
  const kernel = createIntegralKernel(reference.parameters);
  let comparisons = 0, maxError = 0;
  actual.surfaces.forEach((s, part) => {
    assert.equal(s.transition, reference.native.surfaces[part].transition);
    assert.equal(s.forced, reference.native.surfaces[part].forced);
    assert.ok(Math.abs(s.s - reference.native.surfaces[part].s) < 1e-12);
  });
  [...actual.surfaces.map(s => s.states), actual.wake].forEach((profile, part) => {
    const expected = [...reference.native.surfaces.map(s => s.states), reference.native.wake][part];
    assert.equal(profile.length, expected.length);
    profile.forEach((p, i) => {
      assert.equal(p.s, expected[i].s);
      for (const key of ['ue', 'theta', 'deltaStar', 'aux']) {
        const error = Math.abs(p[key] - expected[i][key]) / (key === 'aux' ? Math.max(.01, Math.abs(expected[i][key])) : Math.abs(expected[i][key]));
        comparisons++; maxError = Math.max(maxError, error);
        assert.ok(error < 1e-10, `${part}/${i}/${key}: ${error}`);
      }
      const laminar = part < 2 && i < actual.surfaces[part].transition;
      assert.ok(laminar ? p.aux >= 0 && p.aux < ncrit : p.aux > 0);
      kernel.station(p, part === 2 ? 'wake' : laminar ? 'laminar' : 'turbulent');
    });
  });
  const unchanged = relaxXfoilProfiles(original.input, reference.parameters);
  assert.deepEqual(actual.surfaces, unchanged.surfaces); assert.deepEqual(actual.wake, unchanged.wake);
  assert.equal(comparisons, 524);
  t.diagnostic(JSON.stringify({ comparisons, maxError, phases: actual.surfaces.map(s => s.transition), inputUnchanged: true }));
});

test('strict returned auxiliary/phase checks reject invalid native output', async () => {
  const url = new URL('../src/viscous/xfoil-profile-relaxation.js', import.meta.url);
  const source = fs.readFileSync(url, 'utf8');
  const call = 'ensureCtx(ctx); blpini(ctx); mrchdu(ctx);';
  assert.equal(source.split(call).length, 2);
  for (const statement of ['ctx.CTAU[2][1] = -1e-6;', 'ctx.CTAU[2][1] = ctx.AMCRIT;', 'ctx.CTAU[ctx.ITRAN[1]][1] = 0;',
    'ctx.CTAU[ctx.IBLTE[2] + 1][2] = 0;']) {
    // No march: inject an invalid returned value and exercise the actual
    // production output checks, leaving every input and gas check intact.
    const text = source.replace(call, `ensureCtx(ctx); blpini(ctx); ${statement}`)
      .replace(/from '(\.\/[^']+)'/g, (_, p) => `from '${new URL(p, url).href}'`);
    const module = await import(`data:text/javascript;base64,${Buffer.from(text).toString('base64')}`);
    const input = structuredClone(original.input), before = structuredClone(input);
    assert.throws(() => module.relaxXfoilProfiles(input, original.parameters), /returned (auxiliary state inconsistent|nonpositive wake shear)/);
    assert.deepEqual(input, before);
  }
});
