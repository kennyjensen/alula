// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';

// Controller test: synthetic barrier at high Mach; this is not a PDE validation.
test('one bounded lower-Mach detour advances alpha and restores requested Mach', async () => {
  const families = { euler: 0, boundaryLayer: 0, edgeMatching: 0 };
  const checkpoint = { families, restart: { input: { mach: .74, alpha: 2.1 }, options: {} }, continuation: {} };
  const root = cp => ({ converged: true, mesh: { quality: { valid: true } }, families,
    checkpoint: structuredClone(cp), history: [] });
  const transfers = [];
  globalThis.__detour = {
    initializeCoupledStreamtubeFromFlow: (mach, cp, { targetAlpha } = {}) => {
      transfers.push({ source: structuredClone(cp.restart.input), mach, alpha: targetAlpha });
      const next = structuredClone(cp); next.restart.input.mach = mach;
      if (targetAlpha !== undefined) next.restart.input.alpha = targetAlpha;
      return { checkpoint: next, system: {} };
    },
    solveCoupledStreamtubeIses: (_, { resume, maxIterations }) => {
      const p = resume.restart.input;
      const r = root(resume);
      // Crossing alpha=2.1 at .74 fails, but approaching at .72 succeeds.
      if (maxIterations && p.mach > .73 && p.alpha > 2.1 && p.alpha < 2.29)
        return { ...r, converged: false, reason: 'synthetic barrier', families: { ...families, boundaryLayer: 1 } };
      return r;
    },
    transitionRecoveryPlan: () => null,
  };
  try {
    const code = fs.readFileSync('src/euler/tests/streamtube-coupled-alpha.js', 'utf8')
      .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis.__detour;`);
    const { solveCoupledStreamtubeAlpha: solve } = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
    const r = solve(2.3, { initialCheckpoint: checkpoint, targetMach: .74, maxSubdivisions: 0 });
    assert.equal(r.converged, true);
    assert.equal(r.mach, .74); assert.equal(r.alpha, 2.3);
    const a = r.operatingPointContinuation.attempts;
    assert.equal(a.filter(x => x.stepMethod === 'lower-mach-detour').length, 1);
    assert.ok(a.some(x => x.mach === .72 && x.alpha > 2.1 && x.accepted));
    assert.deepEqual(transfers[2].source, checkpoint.restart.input, 'detour starts at retained root');
    assert.deepEqual(checkpoint.restart.input, { mach: .74, alpha: 2.1 });
  } finally { delete globalThis.__detour; }
});

test('alpha progress preserves Mach-transfer rejection history', async () => {
  const families = { euler: 0, boundaryLayer: 0, edgeMatching: 0 };
  const checkpoint = { restart: { input: { mach: .55, alpha: 1 }, options: {} }, continuation: {} };
  globalThis.__machCap = {
    initializeCoupledStreamtubeFromFlow: (mach, cp, { targetAlpha } = {}) => {
      // Reproduce a transfer that is inadmissible before Newton can start.
      if (mach - cp.restart.input.mach > .02) throw new Error('inadmissible transferred Hk');
      const next = structuredClone(cp); next.restart.input.mach = mach;
      next.restart.input.alpha = targetAlpha ?? cp.restart.input.alpha;
      return { checkpoint: next, system: {} };
    },
    solveCoupledStreamtubeIses: (_, { resume }) => ({ converged: true, mesh: { quality: { valid: true } },
      checkpoint: resume, families, history: [] }),
    transitionRecoveryPlan: () => null,
  };
  try {
    const code = fs.readFileSync('src/euler/tests/streamtube-coupled-alpha.js', 'utf8')
      .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis.__machCap;`);
    const { solveCoupledStreamtubeAlpha: solve } = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
    const r = solve(1.3, { initialCheckpoint: checkpoint, targetMach: .60, preferAlphaApproach: true, maxStages: 20 });
    assert.equal(r.converged, true, r.reason);
    const a = r.operatingPointContinuation.attempts;
    assert.equal(a[0].accepted, false);
    assert.equal(a[1].accepted, false);
    assert.equal(a[2].stepMethod, 'alpha'); assert.equal(a[2].accepted, true);
    assert.equal(a[3].stepMethod, 'mach'); assert.equal(a[3].accepted, true);
    assert.ok(a[3].mach < a[1].mach, 'successful alpha step must not reset rejected Mach increment');
    assert.deepEqual(checkpoint.restart.input, { mach: .55, alpha: 1 });
  } finally { delete globalThis.__machCap; }
});
