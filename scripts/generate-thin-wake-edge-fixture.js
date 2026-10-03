// SPDX-License-Identifier: GPL-2.0-or-later
// Retain the independently checked scalar stencil as a fast regression.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { changedSources, sha256 } from './validation/provenance.js';
const source = 'docs/current-natural-wall-edge-decimal.json';
const r = JSON.parse(fs.readFileSync(source));
assert.equal(r.passed, true); assert.deepEqual(changedSources(r.sourceHashes), []);
assert.ok(r.checks.every(c => c.error < 1e-8));
const output = 'tests/fixtures/thin-wake-edge-decimal.json';
const fixture = { scope: 'Two thin wake-bank cells and one physical directional perturbation. Independent 80-digit residual/derivative reference; no complete-flow or physical-accuracy claim.',
  input: r.input, expected: { residual: r.oracle.residual, derivative: r.oracle.checks.at(-1).derivative },
  provenance: { source, sourceHash: sha256(source), sourceHashes: r.sourceHashes,
    generator: 'scripts/generate-thin-wake-edge-fixture.js', generatorHash: sha256('scripts/generate-thin-wake-edge-fixture.js') } };
fs.writeFileSync(output, JSON.stringify(fixture, null, 2) + '\n');
console.log(JSON.stringify({ output, bytes: fs.statSync(output).size }));
