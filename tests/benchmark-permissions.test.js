import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { benchmarkRestriction } from '../src/ui/tests/benchmark-permissions.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { builtinAirfoils } from '../src/geometry/builtin-airfoils.js';

const presets = ['single','flap','three','rae2822-mses','nlr7301','30p30n'];
const visiblePresets = presets.filter(id => id !== '30p30n');
const modes = ['inviscid','coupled','streamtube-grid','streamtube-bl'];
test('all six available airfoil presets allow attempts in all four modes', () => {
  for (const id of presets) {
    const preset = getBenchmarkAirfoil(id) ?? { elements: builtinAirfoils[id] }, before = structuredClone(preset);
    for (const mode of modes) assert.equal(benchmarkRestriction(preset, mode, { quadBoundaryLayers: true }), '');
    assert.deepEqual(preset, before);
  }
});
test('catalog geometry carries no solver-specific availability metadata', () => {
  for (const id of presets) {
    const preset = getBenchmarkAirfoil(id);
    if (preset) for (const key of ['solverSupportedModels','solverSupportedQuadModes','solverUnsupportedReason'])
      assert.equal(preset[key], undefined);
  }
});
test('historical metadata cannot restore a removed lock', () => {
  const old = { solverSupportedModels: [], solverSupportedQuadModes: [], solverUnsupportedReason: 'Old restriction' };
  for (const mode of modes) assert.equal(benchmarkRestriction(old, mode, { quadBoundaryLayers: true }), '');
  assert.equal(benchmarkRestriction(null, 'coupled'), '');
});

test('the menu offers one RAE geometry with its unchanged dataset identity', () => {
  const html = fs.readFileSync(new URL('../index.html', import.meta.url), 'utf8');
  const selector = html.match(/<select id="preset"[^>]*>([\s\S]*?)<\/select>/)[1];
  const options = [...selector.matchAll(/<option value="([^"]+)"([^>]*)>([^<]+)<\/option>/g)]
    .filter(([, , attributes]) => !attributes.includes('disabled'))
    .map(([, id, , label]) => ({ id, label }));
  assert.deepEqual(options.map(o => o.id).sort(), [...visiblePresets].sort());
  assert.deepEqual(options.filter(o => o.label.includes('RAE')), [{ id: 'rae2822-mses', label: 'RAE2822' }]);
  const visible = getBenchmarkAirfoil('rae2822-mses'), retained = getBenchmarkAirfoil('rae2822');
  assert.equal(visible.name, 'RAE 2822'); assert.equal(visible.elements[0].name, 'RAE 2822');
  assert.equal(visible.id, 'rae2822-mses'); assert.equal(retained.id, 'rae2822');
  assert.equal(visible.provenance.sourcePath, 'third_party/airfoils/rae2822-mses/blade.rae');
  assert.equal(retained.provenance.sourcePath, 'third_party/airfoils/rae2822/geom.txt');
  assert.notDeepEqual(visible.elements[0].points, retained.elements[0].points);
});
