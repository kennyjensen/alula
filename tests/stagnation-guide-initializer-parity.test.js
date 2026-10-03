// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { pathToFileURL } from 'node:url';
import { resolve } from 'node:path';
import { isDeepStrictEqual } from 'node:util';
import { createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const archive = 'docs/rae2822/guide-reference-current/before/';
function moduleURL(text, original, replacements = {}) {
  const base = pathToFileURL(resolve(original));
  return 'data:text/javascript;base64,' + Buffer.from(text.replace(/\bfrom\s+(['"])(\.[^'"]+)\1/g,
    (_, quote, name) => `from ${quote}${replacements[name] ?? new URL(name, base).href}${quote}`)).toString('base64');
}
const matcher = moduleURL(fs.readFileSync(archive + 'streamtube-potential-guides.js.txt', 'utf8'), 'src/euler/streamtube-potential-guides.js');
const old = await import(moduleURL(fs.readFileSync(archive + 'streamtube-body-initializer.js.txt', 'utf8'), 'src/euler/streamtube-body-initializer.js',
  { './streamtube-potential-guides.js': matcher }));

test('successful existing incoming guides keep their complete initialized geometry and diagnostics bit-identical', () => {
  const input = intrinsicBodyFixture({ elements: 1, alpha: 2, bodySegments: 8, tubes: 3, contourPanels: 80 });
  const before = structuredClone(input);
  const controls = { crosslinePlacement: 'potential', surfaceCrosslineMetric: 'arc', cutStationSpacing: 'physical-x',
    outerCrosslineSpread: .85, resolveSurfaceTurning: true, surfaceStationPlacement: 'passage-density',
    passageCountPlanning: true, interiorInitialization: 'linear' };
  const a = old.createPanelStreamtubeGrid(input, controls), b = createPanelStreamtubeGrid(input, controls);
  assert.ok(!b.diagnostics.geometricStagnationConnectors, 'ordinary successful guide unexpectedly used fallback');
  assert.ok(isDeepStrictEqual(input, before), 'caller input changed');
  for (const key of ['input', 'nodes', 'initial', 'diagnostics'])
    assert.ok(isDeepStrictEqual(a[key], b[key]), `existing successful ${key} changed`);
});
