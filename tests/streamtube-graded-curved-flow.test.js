// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';

test('graded thin/skew curved compressible cells retain smooth-field accuracy in both flow directions', t => {
  const temp = fs.mkdtempSync(path.join(os.tmpdir(), 'mses-graded-curved-'));
  try {
    const output = path.join(temp, 'study');
    const run = spawnSync(process.execPath, ['scripts/check-graded-curved-flow.js', `--output=${output}`],
      { cwd: new URL('..', import.meta.url), encoding: 'utf8', timeout: 10000, maxBuffer: 1024 * 1024 });
    assert.equal(run.status, 0, run.stderr || run.stdout || String(run.error));
    const r = JSON.parse(fs.readFileSync(path.join(output, 'report.json')));
    assert.equal(r.passed, true);
    assert.deepEqual(r.failures, []);
    assert.deepEqual(r.sourceChanges, []);
    assert.equal(r.ringleb.filter(f => f.gated !== false).length, 12);
    assert.equal(r.vortex.length, 2);
    assert.equal(r.operations.globalFlowEvaluations + r.operations.BLClosureEvaluations + r.operations.jacobians
      + r.operations.linearSolves + r.operations.newtonUpdates, 0);
    // A finite spacing jump does not become a smooth grid merely because all
    // intervals shrink. Keep it visible instead of counting it as a passing
    // second-order family or relaxing the smooth-family order requirement.
    const contrasts = r.ringleb.filter(f => f.gated === false);
    assert.equal(contrasts.length, 2);
    for (const f of contrasts) assert(f.orders.R1Density.at(-1) < 1.5);
    t.diagnostic(JSON.stringify({ seconds: r.seconds, localCells: r.operations.explicitLocalCells,
      smoothFamilies: 14, spacingJumpOrders: contrasts.map(f => f.orders.R1Density.at(-1)) }));
  } finally {
    fs.rmSync(temp, { recursive: true, force: true });
  }
});
