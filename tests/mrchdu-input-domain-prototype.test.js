import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { validateMrchduMarchInput } from '../scripts/validation/mrchdu-input-domain.js';

const kernel = createIntegralKernel({ reynolds: 2.7e6, mach: .74, ncrit: 4 });
const state = { s: .5, ue: 1, theta: .001, deltaStar: .00105, aux: .03 };

test('only raw-Hk-invalid local march input is admitted, with physical variables unchanged', () => {
  for (const gap of [0, .002]) {
    const p = { ...state, deltaStar: state.deltaStar + gap, wakeGap: gap }, before = structuredClone(p);
    const value = validateMrchduMarchInput(kernel, p, 'wake');
    assert.equal(value.marchInputOnly, true); assert.ok(value.rawHk < 1);
    assert.deepEqual(value.failedChecks, ['raw-hk']); assert.deepEqual(p, before);
    assert.throws(() => kernel.station(p, 'wake'), error => error.code === 'BL_EDGE_STATE_DOMAIN');
  }
});

test('good inputs retain the exact strict-kernel value', () => {
  const p = { ...state, deltaStar: .0025 };
  assert.deepEqual(validateMrchduMarchInput(kernel, p, 'turbulent'), kernel.station(p, 'turbulent'));
});

test('thermal, primitive and unexpected errors remain hard failures', () => {
  for (const p of [{ ...state, ue: 10 }, { ...state, theta: -1 }, { ...state, ue: NaN },
    { ...state, deltaStar: state.theta }, { ...state, wakeGap: -1 }]) {
    const before = structuredClone(p);
    assert.throws(() => validateMrchduMarchInput(kernel, p, 'wake'));
    assert.deepEqual(p, before);
  }
  const unexpected = new Error('unrelated');
  assert.throws(() => validateMrchduMarchInput({ station: () => { throw unexpected; } }, state, 'wake'), error => error === unexpected);
});

test('a native-returned raw-Hk-invalid profile still fails the untouched output checks', async () => {
  const url = new URL('../src/viscous/xfoil-profile-relaxation.js', import.meta.url);
  let text = fs.readFileSync(url, 'utf8');
  const nativeImport = "import { mrchdu, blpini } from './xfoil/xbl.js';";
  assert.equal(text.split(nativeImport).length, 2);
  // A controlled stub changes one returned physical thickness. It does not
  // execute MRCHDU or replace the production output validation loop.
  text = text.replace(nativeImport, "import { blpini } from './xfoil/xbl.js';\nconst mrchdu = ctx => { ctx.DSTR[2][1] = 1.05 * ctx.THET[2][1]; };");
  text = text.replace("    kernel.station(p, i < phases[side] ? 'laminar' : 'turbulent');",
    "    validateMrchduMarchInput(kernel, p, i < phases[side] ? 'laminar' : 'turbulent');");
  text = text.replace("      kernel.station(p, 'wake');", "      validateMrchduMarchInput(kernel, p, 'wake');");
  text = `import { validateMrchduMarchInput } from '${new URL('../scripts/validation/mrchdu-input-domain.js', import.meta.url).href}';\n`
    + text.replace(/from '(\.\/[^']+)'/g, (_, p) => `from '${new URL(p, url).href}'`);
  const module = await import(`data:text/javascript;base64,${Buffer.from(text).toString('base64')}`);
  const profile = [{ ...state, deltaStar: .0025, s: .1, aux: 0 }, { ...state, deltaStar: .0025, s: .2 }];
  const input = { surfaces: [profile, structuredClone(profile)], phases: [1, 1] }, before = structuredClone(input);
  assert.throws(() => module.relaxXfoilProfiles(input, kernel.parameters), error =>
    error.code === 'BL_EDGE_STATE_DOMAIN' && /raw kinematic shape factor/.test(error.message)
    && error.diagnostics.surfaces[0].states[0].deltaStar === 1.05 * state.theta);
  assert.deepEqual(input, before);
});
