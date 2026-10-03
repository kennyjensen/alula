import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { relaxXfoilProfiles } from '../src/viscous/xfoil-profile-relaxation.js';
import { prepareCoupledMrchduProfiles } from '../src/euler/streamtube-coupled-mrchdu-predictor.js';
import { createIntegralKernel } from '../src/viscous/integral.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-input-shape.json', import.meta.url)));
const sourceURL = new URL('../src/viscous/xfoil-profile-relaxation.js', import.meta.url);
const importText = text => import(`data:text/javascript;base64,${Buffer.from(text.replace(/from '(\.\/[^']+)'/g,
  (_, p) => `from '${new URL(p, sourceURL).href}'`)).toString('base64')}`);

test('raw-shape starting profile matches original native MRCHDU and reports input-only recovery', () => {
  for (const [p, expected] of Object.entries(fixture.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(p)).digest('hex'), expected);
  const input = structuredClone(fixture.input), before = structuredClone(input);
  const result = relaxXfoilProfiles(input, fixture.parameters), kernel = createIntegralKernel(fixture.parameters);
  assert.deepEqual(input, before); assert.deepEqual(result.localConvergenceWarnings, []);
  assert.equal(result.inputShapeRecovery.count, fixture.inputRecoveries.length);
  assert.equal(result.inputShapeRecovery.minimumRawHk, Math.min(...fixture.inputRecoveries.map(s => s.rawHk)));
  assert.equal(result.inputShapeRecovery.method, 'native-mrchdu-local-shape');
  const expectedLocations = fixture.inputRecoveries.map(s => [s.part === 2 ? 'wake' : 'surface',
    s.part === 2 ? undefined : ['upper', 'lower'][s.part], s.index]);
  assert.deepEqual(result.inputShapeRecovery.stations.map(s => [s.part, s.side, s.index]), expectedLocations);
  let comparisons = 0;
  [...result.surfaces.map(s => s.states), result.wake].forEach((profile, part) => {
    const expected = part === 2 ? fixture.native.wake : fixture.native.surfaces[part].states;
    assert.equal(profile.length, expected.length);
    profile.forEach((p, i) => {
      assert.equal(p.s, expected[i].s);
      for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
        const scale = key === 'aux' ? Math.max(.01, Math.abs(expected[i][key])) : Math.abs(expected[i][key]);
        assert.ok(Math.abs(p[key] - expected[i][key]) / scale < 1e-10); comparisons++;
      }
      assert.ok(kernel.station(p, part === 2 ? 'wake' : i < result.surfaces[part].transition ? 'laminar' : 'turbulent').rawHk > 1);
    });
  });
  assert.equal(comparisons, 780);
  result.surfaces.forEach((s, i) => {
    assert.equal(s.transition, fixture.native.surfaces[i].transition);
    assert.equal(s.forced, fixture.native.surfaces[i].forced);
    assert.ok(Math.abs(s.s - fixture.native.surfaces[i].s) < 1e-12);
  });
});

test('valid input retains the complete archived numerical result and omits new diagnostics', async () => {
  const oldText = fs.readFileSync(new URL('../docs/rae2822/mrchdu-input-domain487/runtime-promotion/xfoil-profile-relaxation.js.before.txt', import.meta.url), 'utf8');
  const old = await importText(oldText);
  const controls = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-body.json', import.meta.url)));
  const good = controls.cases.find(c => c.name.startsWith('default-two-element-root'));
  const expected = old.relaxXfoilProfiles(good.input, good.parameters), actual = relaxXfoilProfiles(good.input, good.parameters);
  assert.deepEqual(actual, expected); assert.equal(Object.hasOwn(actual, 'inputShapeRecovery'), false);
});

async function stubbed(body) {
  const text = fs.readFileSync(sourceURL, 'utf8'), line = "import { mrchdu, blpini } from './xfoil/xbl.js';";
  assert.equal(text.split(line).length, 2);
  return importText(text.replace(line, `import { blpini } from './xfoil/xbl.js';\nconst mrchdu = ctx => { ${body} };`));
}
const parameters = { reynolds: 2.7e6, mach: .74, ncrit: 4 };
function simpleInput() {
  const profile = [{ s: .1, ue: 1, aux: 0, theta: .001, deltaStar: .0025 },
    { s: .2, ue: 1, aux: .03, theta: .001, deltaStar: .0025 }];
  return { surfaces: [profile, structuredClone(profile)], phases: [1, 1] };
}

test('thermal and primitive input errors reject before the local march', async () => {
  const module = await stubbed("throw new Error('Unexpected local march');");
  for (const change of [p => { p.ue = 10; }, p => { p.ue = NaN; }, p => { p.theta = -1; },
    p => { p.deltaStar = p.theta; }, p => { p.wakeGap = -1; }]) {
    const input = simpleInput(); change(input.surfaces[0][0]); const before = structuredClone(input);
    assert.throws(() => module.relaxXfoilProfiles(input, parameters), error => !error.message.includes('Unexpected local march'));
    assert.deepEqual(input, before);
  }
});

test('strict output validation rejects an invalid native-returned physical shape', async () => {
  const module = await stubbed('ctx.DSTR[2][1] = 1.05 * ctx.THET[2][1];');
  const input = simpleInput(), before = structuredClone(input);
  assert.throws(() => module.relaxXfoilProfiles(input, parameters), error =>
    error.code === 'BL_EDGE_STATE_DOMAIN' && /raw kinematic shape factor/.test(error.message)
    && error.diagnostics.surfaces[0].states[0].deltaStar === 1.05 * input.surfaces[0][0].theta);
  assert.deepEqual(input, before);
});

test('wake-input exception uses fluid shape while preserving the supplied finite gap', async () => {
  const module = await stubbed(`for (let i = ctx.IBLTE[2] + 1; i <= ctx.NBL[2]; i++) {
    const gap = ctx.WGAP[i - ctx.IBLTE[2]];
    if (ctx.DSTR[i][2] !== 1.05 * ctx.THET[i][2] + gap) throw new Error('Native input changed');
    ctx.DSTR[i][2] = 2.5 * ctx.THET[i][2] + gap;
  }`);
  const input = simpleInput(); input.normalGap = .002;
  input.wake = [.2, .3].map(s => ({ s, ue: 1, aux: .03, theta: .001, deltaStar: 1.05 * .001 + .002, wakeGap: .002 }));
  const before = structuredClone(input), result = module.relaxXfoilProfiles(input, parameters);
  assert.deepEqual(input, before); assert.equal(result.inputShapeRecovery.count, 2);
  assert.ok(result.inputShapeRecovery.stations.every(s => s.part === 'wake'));
  result.wake.forEach(p => { assert.equal(p.wakeGap, .002); assert.equal(p.deltaStar, 2.5 * .001 + .002); });
});

test('coupled predictor forwards detached input-shape diagnostics without changing source profiles', () => {
  const { input, parameters } = fixture, states = [], initialBL = [];
  const put = profile => profile.map(p => { const id = states.length; states.push(structuredClone(p));
    initialBL.push(p.aux, p.theta, p.deltaStar, p.ue); return id; });
  const surfaces = input.surfaces.map((p, i) => ({ body: 0, side: ['upper', 'lower'][i], transition: input.phases[i], ids: put(p) }));
  const wakes = [{ body: 0, ids: put(input.wake) }];
  const request = { bl: { transitionMode: 'automatic', surfaces, wakes, scale: 1, trips: [[1, 1]], kernel: { parameters } },
    states, initialBL, targetMach: parameters.mach };
  const before = structuredClone(request), result = prepareCoupledMrchduProfiles(request);
  assert.deepEqual(request, before);
  assert.equal(result.diagnostics.bodies[0].inputShapeRecovery.count, fixture.inputRecoveries.length);
  assert.deepEqual(result.transitionState, fixture.native.surfaces.map(s => s.transition));
  result.diagnostics.bodies[0].inputShapeRecovery.stations[0].rawHk = 99;
  assert.deepEqual(request, before);
});
