// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded validation of MSES's temporary dissipation schedule on a saved state.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { isMainThread, parentPort, workerData } from 'node:worker_threads';
import { runBoundedWorker } from './validation/bounded-worker.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';

const serial = v => JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x);

// BEGIN PURE SAVED-INPUT HELPERS (also exercised without importing the flow driver).
function loadShockBroadeningInput(saved) {
  assert.ok(saved && typeof saved === 'object', 'Supply saved checkpoint or controller data.');
  const controllers = [];
  if (saved.kind === 'mses-shock-broadening') controllers.push(saved);
  if (Object.hasOwn(saved, 'controller')) controllers.push(saved.controller);
  if (saved.result && Object.hasOwn(saved.result, 'controller')) controllers.push(saved.result.controller);
  for (const controller of controllers) {
    assert.ok(controller?.kind === 'mses-shock-broadening' && controller.version === 1
      && controller.checkpoint?.version === 1 && controller.checkpoint.restart && controller.checkpoint.continuation,
    'A saved controller must be complete; do not discard its density history.');
    assert.deepEqual(controller, controllers[0], 'Conflicting saved shock-broadening controllers.');
  }
  if (controllers.length) return controllers[0];
  const checkpoint = saved.checkpoint ?? saved.result?.checkpoint ?? saved;
  assert.ok(checkpoint?.version === 1 && checkpoint.restart && checkpoint.continuation,
    'Supply a complete coupled checkpoint or preserved shock-broadening controller.');
  return checkpoint;
}

function requireAdjacentShockCheckpoints(previousSaved, currentSaved) {
  const read = saved => {
    const checkpoint = loadShockBroadeningInput(saved);
    assert.notEqual(checkpoint.kind, 'mses-shock-broadening', 'A saved controller already owns its density history.');
    const candidates = [['details.history', saved.details?.history], ['history', saved.history],
      ['result.history', saved.result?.history], ['h', Array.isArray(saved.h) ? saved.h : undefined]];
    const available = candidates.filter(([, value]) => value !== undefined);
    assert.ok(available.length, 'Adjacent density history requires retained complete iteration histories.');
    const [path, history] = available[0];
    assert.ok(Array.isArray(history) && history.length > 0, 'Retained iteration history must be nonempty.');
    for (const [, other] of available) assert.deepEqual(other, history, 'Conflicting retained iteration histories.');
    for (let i = 0; i < history.length; i++) {
      assert.ok(Number.isInteger(history[i]?.iteration) && history[i].iteration >= 0,
        'Retained history needs nonnegative integer iteration numbers.');
      if (i) assert.equal(history[i].iteration, history[i - 1].iteration + 1, 'Retained history must be contiguous.');
    }
    const last = history.at(-1);
    if (saved.h && !Array.isArray(saved.h))
      assert.equal(saved.h.iteration, last.iteration, 'Saved h iteration disagrees with retained history.');
    assert.ok(checkpoint.families && Object.keys(checkpoint.families).length > 0,
      'Bind retained history to the checkpoint residual families.');
    for (const [family, value] of Object.entries(checkpoint.families)) {
      assert.ok(Number.isFinite(value), 'Checkpoint residual families must be finite.');
      assert.equal(last[family], value, `Retained ${family} history does not match its checkpoint.`);
    }
    return { checkpoint, history, path, iteration: last.iteration };
  };
  const previous = read(previousSaved), current = read(currentSaved);
  assert.equal(current.iteration, previous.iteration + 1, 'Density history must come from adjacent accepted updates.');
  assert.equal(current.history.length, previous.history.length + 1, 'Current history must extend the previous history once.');
  assert.deepEqual(current.history.slice(0, -1), previous.history, 'Adjacent accepted histories must share an exact prefix.');
  return { previousIteration: previous.iteration, currentIteration: current.iteration,
    previousHistoryPath: previous.path, currentHistoryPath: current.path,
    previousHistoryEntries: previous.history.length, currentHistoryEntries: current.history.length,
    historyPrefixExact: true, checkpointFamiliesExact: true };
}
// END PURE SAVED-INPUT HELPERS

if (isMainThread) {
  const args = Object.fromEntries(process.argv.slice(2).map(arg => {
    assert.match(arg, /^--[^=]+=.*/); const i = arg.indexOf('=');
    return [arg.slice(2, i), arg.slice(i + 1)];
  }));
  assert.ok(args.input && args.output);
  assert.ok(Object.keys(args).every(k => ['input', 'previous', 'output', 'iterations', 'timeout'].includes(k)));
  assert.equal(fs.existsSync(args.output), false, 'Preserve existing evidence.');
  const controls = { maxIterations: Number(args.iterations ?? 4), tolerance: 1e-10 };
  const timeoutMs = Number(args.timeout ?? 45000);
  assert.ok(Number.isInteger(controls.maxIterations) && controls.maxIterations >= 0 && controls.maxIterations <= 40);
  assert.ok(Number.isFinite(timeoutMs) && timeoutMs > 0);
  const saved = JSON.parse(fs.readFileSync(args.input));
  let controller = loadShockBroadeningInput(saved);
  if (controller.kind === 'mses-shock-broadening' && controller.tolerance !== undefined)
    controls.tolerance = controller.tolerance;
  assert.ok(Number.isFinite(controls.tolerance) && controls.tolerance > 0, 'Preserve a valid saved controller tolerance.');
  const inputs = { [args.input]: sha256(args.input) };
  let densityHistory;
  if (args.previous) {
    assert.notEqual(controller.kind, 'mses-shock-broadening', 'A saved controller already owns its density history.');
    const previousSaved = JSON.parse(fs.readFileSync(args.previous)), prior = loadShockBroadeningInput(previousSaved), cp = controller;
    const adjacency = requireAdjacentShockCheckpoints(previousSaved, saved);
    assert.equal(serial(prior.restart.input), serial(cp.restart.input), 'Adjacent states must use the same Euler equations.');
    const withoutPhase = options => { const next = { ...options }; delete next.transitionState; return next; };
    assert.equal(serial(withoutPhase(prior.restart.options)), serial(withoutPhase(cp.restart.options)),
      'Adjacent states must use the same physical BL conditions.');
    const nodes = cp.restart.initialEuler.nodes;
    const densityCount = (nodes[0].length - 1) * nodes.reduce((n, g) => n + g[0].length - 1, 0);
    const { maximumAppliedDensityChange, temporaryShockMcrit } = await import('./validation/coupled-shock-broadening.js');
    const d = maximumAppliedDensityChange(prior.restart.initialEuler.x, cp.restart.initialEuler.x, densityCount);
    controller = { version: 1, kind: 'mses-shock-broadening', checkpoint: cp,
      targetMcrit: cp.restart.input.upwind.mcrit, previousDensityChange: d, acceptedUpdates: 0 };
    densityHistory = { previous: args.previous, current: args.input, ...adjacency, densityCount, appliedDensityChange: d,
      scheduledMcrit: temporaryShockMcrit({ targetMcrit: controller.targetMcrit, densityChange: d }),
      convention: 'Maximum actual accepted physical density ratio change; the manual does not specify damped versus undamped.' };
    inputs[args.previous] = sha256(args.previous);
  }
  const output = args.output; fs.mkdirSync(output, { recursive: true });
  const write = (name, value) => fs.writeFileSync(`${output}/${name}.json`, serial(value) + '\n');
  write('input', controller);
  const derivedInput = { path: `${output}/input.json`, sha256: sha256(`${output}/input.json`) };
  inputs[derivedInput.path] = derivedInput.sha256;
  const report = { inProgress: true, passed: false, physicalAcceptance: false, controls, timeoutMs, inputs, derivedInput, densityHistory,
    sourceHashes: numericalSourceHashes(['scripts/check-coupled-shock-broadening.js', 'scripts/validation/bounded-worker.js',
      'scripts/validation/coupled-shock-broadening.js', 'scripts/validation/coupled-shock-broadening-driver.js']),
    scope: 'Validation-only temporary Mcrit continuation using the unchanged ordinary simultaneous Newton driver. No new BL march, mesh initialization or residual tolerance. Success requires exact requested-threshold equations and final convexity.' };
  write('report', report);
  const terminal = await runBoundedWorker(new URL(import.meta.url), { output, controls }, { timeoutMs,
    onProgress: p => { fs.appendFileSync(`${output}/progress.jsonl`, serial(p) + '\n'); console.log(serial(p)); } });
  write('state', terminal);
  const result = terminal.result;
  Object.assign(report, { inProgress: false, status: terminal.type, seconds: terminal.seconds,
    sourceChanges: changedSources(report.sourceHashes), changedInputs: Object.keys(inputs).filter(p => sha256(p) !== inputs[p]),
    inputUnchanged: terminal.inputUnchanged, error: terminal.message,
    result: result && { converged: result.converged, reason: result.reason, targetRestored: result.targetRestored,
      families: result.checkpoint?.families, controller: result.controller && {
        targetMcrit: result.controller.targetMcrit, previousDensityChange: result.controller.previousDensityChange,
        effectiveMcrit: result.controller.checkpoint.restart.input.upwind.mcrit,
        acceptedUpdates: result.controller.acceptedUpdates }, linearDiagnostics: result.linearDiagnostics } });
  report.passed = terminal.type === 'result' && result?.converged === true && result?.targetRestored === true
    && report.inputUnchanged === true && !report.sourceChanges.length && !report.changedInputs.length;
  write('report', report);
  const { sourceHashes, inputs: inputHashes, ...summary } = report;
  write('summary', summary); console.log(serial(summary));
  process.exitCode = report.passed ? 0 : 1;
} else {
  const { output, controls } = workerData;
  const write = (name, value) => fs.writeFileSync(`${output}/${name}.json`, serial(value) + '\n');
  let recordCount = 0;
  try {
    const { solveCoupledWithShockBroadening } = await import('./validation/coupled-shock-broadening-driver.js');
    const input = JSON.parse(fs.readFileSync(`${output}/input.json`)), before = serial(input);
    const result = solveCoupledWithShockBroadening(input, { ...controls,
      onRecord: (kind, value) => {
        const name = `${String(recordCount++).padStart(3, '0')}-${kind.replace(/[^a-z0-9-]/gi, '-')}`;
        write(name, value); parentPort.postMessage({ type: 'flow-stage', stage: kind, record: name });
      },
      onIteration: h => parentPort.postMessage({ type: 'iteration', iteration: {
        iteration: h.iteration, newIteration: h.newIteration, residual: h.residual,
        euler: h.euler, boundaryLayer: h.boundaryLayer, edgeMatching: h.edgeMatching,
        step: h.step, backtracks: h.backtracks, effectiveMcrit: h.effectiveMcrit,
        appliedDensityChange: h.appliedDensityChange, activeChange: h.activeChange,
      } }),
    });
    write('result', result);
    if (result.controller) write('controller', result.controller);
    parentPort.postMessage({ type: 'result', result, inputUnchanged: serial(input) === before });
  } catch (error) {
    write('exception', { message: error.message, stack: error.stack, diagnostics: error.diagnostics,
      controller: error.controller ?? error.controllerCheckpoint });
    parentPort.postMessage({ type: 'error', message: error.message });
  }
}
