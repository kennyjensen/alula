// SPDX-License-Identifier: GPL-2.0-or-later
// Default: parameter/topology planning only. Costly stages require named cases
// or --all. Each stage has a separate outcome; a mesh pass is not convergence.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { isMainThread, parentPort, workerData } from 'node:worker_threads';
import { gridRobustnessCases, buildGridRobustnessCase, planGridRobustnessCase,
  assertGridRobustnessGuiCoverage } from './validation/grid-robustness-cases.js';
import { runBoundedWorker } from './validation/bounded-worker.js';
import { assessSolverResult } from './validation/solver-reliability-acceptance.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';

const serial = x => JSON.stringify(x, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v);
const gates = () => Object.fromEntries(['plan', 'topology', 'geometry', 'gasStart', 'linearSolve', 'nonlinearConvergence', 'physicalAccuracy']
  .map(k => [k, { status: 'not-run' }]));
const writeNew = (p, x) => fs.writeFileSync(p, serial(x) + '\n', { flag: 'wx' });

export function parseGridRobustnessArgs(argv) {
  const options = {};
  for (const arg of argv) {
    const match = /^--([a-z-]+)=(.+)$/.exec(arg);
    const key = arg === '--all' ? 'all' : match?.[1];
    if (!key || !['stage', 'case', 'out', 'seconds', 'mode', 'all', 'max-iterations'].includes(key)
      || key in options) throw new Error(`Invalid or duplicate argument: ${arg}`);
    options[key] = key === 'all' ? true : match[2];
  }
  const stage = options.stage ?? 'plan', mode = options.mode ?? 'coupled';
  if (!['plan', 'topology', 'geometry', 'startup', 'solve'].includes(stage)) throw new Error('Unknown robustness stage.');
  if (!['coupled', 'euler'].includes(mode)) throw new Error('Mode must be coupled or euler.');
  if (options.all && options.case) throw new Error('Choose named cases or --all, not both.');
  const ids = options.case?.split(',');
  if (ids && (new Set(ids).size !== ids.length || ids.some(id => !gridRobustnessCases.some(c => c.id === id))))
    throw new Error('Unknown or duplicate case ID.');
  if (stage !== 'plan' && !ids && !options.all)
    throw new Error('Name --case=... or explicitly use --all for a stage that constructs numerical geometry or flow.');
  const seconds = Number(options.seconds ?? (stage === 'topology' ? 15 : 180));
  if (!Number.isFinite(seconds) || seconds <= 0 || seconds > 600) throw new Error('Per-case seconds must be in (0,600].');
  const maxIterations = options['max-iterations'] === undefined ? undefined : Number(options['max-iterations']);
  if (maxIterations !== undefined && (stage !== 'solve' || !Number.isInteger(maxIterations) || maxIterations < 0 || maxIterations > 100))
    throw new Error('An explicit iteration cap of 0–100 applies only to solve.');
  return { stage, mode, ids, seconds, maxIterations,
    output: options.out ?? `docs/grid-robustness/${new Date().toISOString().replace(/[:.]/g, '-')}-${stage}` };
}

export function assessGridRobustnessSolve(input, result) {
  const assessment = assessSolverResult(input, result);
  if (input.quadBoundaryLayers) {
    const actualNcrit = result?.checkpoint?.restart?.options?.ncrit ?? result?.actualNcrit ?? result?.conditions?.ncrit;
    assessment.checks.requestedNcrit = actualNcrit === input.ncrit && result?.ncritContinuation?.reachedTarget !== false;
    assessment.actualNcrit = actualNcrit;
  }
  assessment.passed = Object.values(assessment.checks).every(Boolean);
  assessment.failures = Object.keys(assessment.checks).filter(k => !assessment.checks[k]);
  return assessment;
}

export async function runGridRobustness(options) {
  const { stage, mode, ids, seconds, maxIterations, output } = options;
  assertGridRobustnessGuiCoverage(fs.readFileSync('index.html', 'utf8'));
  if (fs.existsSync(output)) throw new Error(`Preserve existing evidence: ${output}`);
  fs.mkdirSync(output, { recursive: true });
  const selected = gridRobustnessCases.filter(c => !ids || ids.includes(c.id));
  const sourceHashes = numericalSourceHashes(['index.html', 'scripts/check-grid-robustness.js',
    'scripts/validation/grid-robustness-cases.js', 'scripts/validation/bounded-solver-worker.js',
    'tests/oracles/streamtube-control-volume-geometry.js']);
  const report = { startedAt: new Date().toISOString(), inProgress: true, stage, mode, output,
    limits: { timeoutSecondsPerCase: seconds, maxIterations: maxIterations ?? 'public default' },
    sourceHashes, selected: selected.map(c => c.id), notRun: gridRobustnessCases.filter(c => !selected.includes(c)).map(c => c.id),
    scope: 'Only the named stage is certified. Planning never builds a mesh. Geometry uses independent volume/embedding checks. Startup tests one public cold seed; solve uses the actual browser Worker route.',
    methodReference: 'docs/MSET_GRID_CONTROLS_AUDIT_2026-09-15.md',
    currentConvergenceEnvelopeCertified: false, physicalAccuracy: 'not-run', cases: [] };
  writeNew(`${output}/start.json`, report);
  const save = () => fs.writeFileSync(`${output}/report.json`, serial(report) + '\n'); save();
  for (const spec of selected) {
    const entry = { id: spec.id, preset: spec.preset, tier: spec.tier, stage, gates: gates(), passed: false };
    const started = performance.now(), directory = `${output}/${spec.id}`; fs.mkdirSync(directory);
    try {
      const input = buildGridRobustnessCase(spec, { mode });
      if (maxIterations !== undefined) input.maxIterations = maxIterations;
      writeNew(`${directory}/input.json`, input); entry.inputSha256 = sha256(`${directory}/input.json`);
      entry.plan = planGridRobustnessCase(input); entry.gates.plan = { status: 'pass' };
      if (stage === 'plan') entry.passed = true;
      else {
        const worker = stage === 'solve' ? new URL('./validation/bounded-solver-worker.js', import.meta.url) : new URL(import.meta.url);
        const outcome = await runBoundedWorker(worker, { caseData: input, stage, directory }, {
          timeoutMs: 1000 * seconds, onProgress: message => {
            fs.appendFileSync(`${directory}/progress.jsonl`, serial(message) + '\n');
            console.log(serial({ case: spec.id, ...message }));
          },
        });
        writeNew(`${directory}/outcome.json`, outcome);
        entry.status = outcome.type; entry.reason = outcome.message; entry.code = outcome.code;
        if (outcome.gates) for (const [name, gate] of Object.entries(outcome.gates))
          if (gate.status !== 'not-run') entry.gates[name] = gate;
        if (stage === 'solve') {
          if (outcome.retained?.mesh) entry.gates.geometry = { status: outcome.retained.mesh.quality?.valid ? 'observed-valid' : 'fail',
            scope: 'Solver mesh-quality observation only; independent saved-node audit was not run by this solve stage.' };
          if (outcome.type === 'result') {
            const assessment = assessGridRobustnessSolve(input, outcome.result);
            entry.gates.nonlinearConvergence = { status: assessment.passed ? 'pass' : 'fail', ...assessment };
            entry.passed = assessment.passed;
            entry.actual = { mach: outcome.result.mach, ncrit: assessment.actualNcrit,
              cells: outcome.result.mesh?.cells?.length, grid: outcome.result.mesh?.initialization?.gridSpacing };
          } else entry.gates.nonlinearConvergence = { status: 'fail', reason: outcome.message ?? outcome.type };
        } else entry.passed = outcome.type === 'result' && outcome.passed === true;
        if (!entry.passed && !entry.reason) entry.reason = outcome.reason ?? outcome.type;
      }
    } catch (error) { entry.reason = error.message; entry.code = error.code;
      if (entry.gates.plan.status !== 'pass') entry.gates.plan = { status: 'fail', reason: error.message }; }
    entry.seconds = (performance.now() - started) / 1000; report.cases.push(entry); save();
    console.log(serial({ case: entry.id, stage, passed: entry.passed, seconds: entry.seconds, gates: entry.gates, reason: entry.reason }));
  }
  report.sourceChanges = changedSources(sourceHashes); report.inProgress = false; report.finishedAt = new Date().toISOString();
  report.passed = report.cases.every(c => c.passed) && !report.sourceChanges.length;
  report.selectedStagePassed = report.passed;
  // Even every selected nonlinear root would not cover continuous controls,
  // other airfoils or physical accuracy. Keep the envelope claim false.
  save(); return report;
}

async function numericalWorker() {
  const { caseData, stage, directory } = workerData, resultGates = gates();
  let currentStage = 'topology';
  try {
    if (stage === 'topology') {
      const { createInitialStreamtubeTopology } = await import('../src/geometry/streamtube-topology.js');
      const t = createInitialStreamtubeTopology(caseData, {
        surfaceIntervals: caseData.gridIntervals, surfaceChordExponent: caseData.gridChordExponent, tubes: caseData.gridTubes,
        upperTubes: caseData.gridUpperTubes, lowerTubes: caseData.gridLowerTubes, gapTubes: caseData.gridGapTubes,
        inletIntervals: caseData.gridInletIntervals, outletIntervals: caseData.gridOutletIntervals,
      });
      writeNew(`${directory}/topology.json`, t);
      resultGates.topology = { status: 'pass', rawIntervals: t.outerLower.length - 1,
        requestedTubes: t.weights.map(g => g.length), gridSpacing: t.gridSpacing,
        scope: 'Raw topology only; panel skeleton, tracing, smoothing, final refinement and gas were not run.' };
      parentPort.postMessage({ type: 'result', passed: true, gates: resultGates }); return;
    }
    currentStage = 'geometry';
    const { prepareStreamtubeAssembly } = await import('../src/euler/streamtube-result.js');
    const { coupledMachPlan } = await import('../src/euler/tests/streamtube-coupled-mach-assembly.js');
    const { directStreamtubeVolumeGeometry } = await import('../tests/oracles/streamtube-control-volume-geometry.js');
    const actualMach = caseData.quadBoundaryLayers ? coupledMachPlan(caseData).sourceMach : caseData.mach;
    const coldCase = { ...caseData, mach: actualMach };
    let lastPublication = -Infinity, initialSaved = false;
    const { system, ...prepared } = prepareStreamtubeAssembly(coldCase, { onMesh: (mesh, meshStage) => {
      if (performance.now() - lastPublication < 10000 && meshStage !== 'initial') return;
      lastPublication = performance.now();
      if (!initialSaved) { writeNew(`${directory}/first-mesh.json`, mesh); initialSaved = true; }
      fs.writeFileSync(`${directory}/latest-mesh.json`, serial(mesh) + '\n');
      parentPort.postMessage({ type: 'flow-stage', stage: `grid-${meshStage}`, cells: mesh.cells.length,
        valid: mesh.quality.valid, minCornerSine: mesh.quality.minCornerSine });
    } });
    writeNew(`${directory}/prepared.json`, prepared);
    const independent = directStreamtubeVolumeGeometry(prepared.nodes);
    writeNew(`${directory}/geometry-audit.json`, independent);
    const smoothingPass = !coldCase.gridEllipticSmoothing || prepared.diagnostics.gridSmoothing?.converged === true;
    resultGates.geometry = { status: prepared.mesh.quality.valid && independent.valid && !independent.concavePrimal.length && smoothingPass ? 'pass' : 'fail',
      cells: prepared.mesh.cells.length, nx: system.layout.nx, tubes: system.layout.tubes,
      quality: prepared.mesh.quality, independentValid: independent.valid, concaveCells: independent.concavePrimal.length,
      smoothingRequested: coldCase.gridEllipticSmoothing, smoothingConverged: prepared.diagnostics.gridSmoothing?.converged ?? null,
      actualMach, requestedMach: caseData.mach };
    if (resultGates.geometry.status !== 'pass') throw new Error('Independent final geometry or requested smoothing gate failed.');
    if (stage === 'startup') {
      currentStage = 'gasStart';
      const { initializeStreamtubeStartup } = await import('../src/euler/streamtube-startup.js');
      const startup = initializeStreamtubeStartup(system, prepared.initial);
      writeNew(`${directory}/startup.json`, { initial: startup.initial, diagnostics: startup.diagnostics });
      resultGates.gasStart = { status: 'pass', actualMach, requestedMach: caseData.mach,
        diagnostics: startup.diagnostics, selectedIsmom: caseData.eulerIsmom,
        scope: 'Complete selected Euler residual at one cold seed; no Euler iterations, BL initialization or Mach continuation.' };
    }
    parentPort.postMessage({ type: 'result', passed: true, gates: resultGates });
  } catch (error) {
    resultGates[currentStage] = { ...resultGates[currentStage], status: 'fail', reason: error.message,
      code: error.code, diagnostics: error.diagnostics };
    parentPort.postMessage({ type: 'error', message: error.message, code: error.code, stage: currentStage, gates: resultGates });
  }
}

if (!isMainThread) await numericalWorker();
else if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const report = await runGridRobustness(parseGridRobustnessArgs(process.argv.slice(2)));
  console.log(serial({ report: `${report.output}/report.json`, selected: report.selected.length, stage: report.stage, passed: report.passed,
    sourceChanges: report.sourceChanges, currentConvergenceEnvelopeCertified: false }));
  if (!report.passed) process.exitCode = 1;
}
