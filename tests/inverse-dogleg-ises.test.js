import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { prepareCoupledInverseBL } from '../scripts/validation/coupled-inverse-bl.js';
import { createInverseBLDogleg } from '../scripts/validation/inverse-bl-dogleg.js';
import { inverseDoglegIsesSource } from '../scripts/validation/inverse-dogleg-ises-source.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));

test('inverse dogleg ISES converges all four natural BLs and two wakes after perturbation and resumes the same trust radius', async t => {
  const folder = fs.mkdtempSync(path.join(os.tmpdir(), 'mses-inverse-dogleg-test-'));
  try {
    const originalPath = path.resolve('src/euler/streamtube-coupled-ises.js');
    const { code } = inverseDoglegIsesSource(fs.readFileSync(originalPath, 'utf8'));
    const destination = path.join(folder, 'ises.mjs');
    fs.writeFileSync(destination, code.replace(/from '([^']+)'/g, (_, p) => `from '${pathToFileURL(path.resolve(path.dirname(originalPath), p)).href}'`));
    const research = (await import(pathToFileURL(destination))).solveCoupledStreamtubeIses;
    const controls = { transitionMode: 'automatic', edgeMatching: 'section-velocity', reynolds: 1e6, ncrit: 9,
      tolerance: 1e-10, maxIterations: 12, stepAcceptance: 'admissible' };
    const root = solveCoupledStreamtubeIses(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), controls);
    assert.equal(root.converged, true, root.reason);
    const checkpoint = serialize(root.checkpoint), f = checkpoint.restart;
    const s = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
    const x = s.initial.slice(); x[s.euler.layout.densityCount - 1] += 1e-4;
    const p = prepareCoupledInverseBL(s, x); s.bl.restoreActive(p.phase); const v = s.evaluate(p.x);
    checkpoint.families = v.families; checkpoint.inverseDogleg = { radius: 1 };
    f.initialEuler = { x: p.x.slice(0, s.ne), nodes: v.outer.nodes, undisplacedNodes: v.outer.undisplacedNodes };
    f.initialBL = p.x.slice(s.ne); f.options.transitionState = p.phase;
    const initial = serialize(checkpoint), researchControls = { ...controls,
      inversePrepare: prepareCoupledInverseBL,
      researchModel: (system, state, value, options) => {
        const matrix = system.jacobian(state);
        const linear = solveSparseDirect(matrix, value.residual.map(v => -v), options);
        return { ...createInverseBLDogleg(system, matrix, value.residual, linear.x), linear };
      } };
    const full = research(undefined, { ...researchControls, resume: initial });
    assert.equal(full.converged, true, full.reason);
    assert.equal(full.boundaryLayer.surfaces.length, 4); assert.equal(full.boundaryLayer.wakes.length, 2);
    assert.ok(full.boundaryLayer.transitions.every(q => q.kind === 'natural'));
    assert.ok(full.history.slice(1).every(h => h.trustRegion.accepted && h.trustRegion.ratio > 1e-4));
    const first = research(undefined, { ...researchControls, resume: initial, maxIterations: 1 });
    assert.equal(first.history.length, 2); assert.ok(first.checkpoint.inverseDogleg.radius > 0);
    const resumed = research(undefined, { ...researchControls, resume: serialize(first.checkpoint) });
    assert.equal(resumed.converged, true, resumed.reason);
    assert.deepEqual(serialize(resumed.x), serialize(full.x));
    assert.deepEqual(serialize(initial), serialize(checkpoint));
    const zero = research(undefined, { ...researchControls, resume: serialize(root.checkpoint) });
    assert.equal(zero.converged, true); assert.equal(zero.linearDiagnostics.solves, 0);
    t.diagnostic(JSON.stringify({ unknowns: s.n, iterations: full.history.length - 1,
      families: full.families, exactResume: true, kinds: full.history.slice(1).map(h => h.stepKind) }));
  } finally { fs.rmSync(folder, { recursive: true, force: true }); }
});
