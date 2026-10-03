// SPDX-License-Identifier: GPL-2.0-or-later
// Plain diagnostic data only. The caller renders summary/message with
// textContent and owns request identity, timestamps and source URL capture.
const finite = value => typeof value === 'number' && Number.isFinite(value);
const firstNumber = (...values) => values.find(finite);
const present = value => value !== undefined && value !== null && value !== '';
const firstValue = (...values) => values.find(present);
const show = value => present(value) ? String(value) : 'unknown';
const number = value => finite(value) ? String(value) : 'unknown';
const automatic = value => present(value) ? String(value) : 'Auto';

// Submitted inputs are JSON-shaped. Handle accidental non-JSON diagnostic
// fields without allowing an error-reporting failure or retaining aliases.
function plain(value, { diagnostic = false } = {}, stack = new Set()) {
  if (value === null || typeof value === 'string' || typeof value === 'boolean') return value;
  if (typeof value === 'number') return Number.isFinite(value) ? value : diagnostic ? String(value) : null;
  if (typeof value === 'bigint') return String(value);
  if (!value || typeof value !== 'object') return undefined;
  if (stack.has(value)) return '[Circular]';
  stack.add(value);
  let result;
  if (Array.isArray(value) || ArrayBuffer.isView(value)) {
    result = Array.from(value, entry => plain(entry, { diagnostic }, stack) ?? null);
  } else {
    // Never embed a full solver result, mesh, flow or restart from an error.
    // The separately captured submitted input retains all of its geometry.
    const excluded = diagnostic ? new Set(['flow', 'mesh', 'checkpoint', 'restart', 'result', 'system', 'initialEuler', 'initialBL',
      'nodes', 'undisplacedNodes', 'vertices', 'cells', 'sections', 'residual']) : null;
    result = Object.fromEntries(Object.entries(value).flatMap(([key, entry]) => {
      if (excluded?.has(key) && !(key === 'residual' && typeof entry === 'number')) return [];
      if (diagnostic && key === 'lastRejectedStep') return [[key, selected(entry, rejectionKeys)]];
      const copied = plain(entry, { diagnostic }, stack);
      return copied === undefined ? [] : [[key, copied]];
    }));
  }
  stack.delete(value);
  return result;
}

function selected(object, keys) {
  if (!object || typeof object !== 'object') return {};
  return Object.fromEntries(keys.flatMap(key => {
    const value = plain(object[key], { diagnostic: true });
    return value === undefined ? [] : [[key, value]];
  }));
}

function errorCause(cause, seen = new Set()) {
  if (cause === undefined) return undefined;
  if (cause === null || typeof cause !== 'object') return plain(cause);
  if (seen.has(cause)) return '[Circular]';
  seen.add(cause);
  const result = selected(cause, ['name', 'message', 'code', 'stage', 'diagnostics', 'reason']);
  if (cause.cause !== undefined) result.cause = errorCause(cause.cause, seen);
  return result;
}

const progressKeys = ['type', 'stage', 'startupAttempt', 'attempt', 'iteration', 'mach', 'actualMach', 'targetMach', 'actualAlpha', 'targetAlpha',
  'requestedMach', 'actualNcrit', 'targetNcrit', 'coarseStage', 'residual', 'euler', 'boundaryLayer', 'edgeMatching', 'step', 'backtracks', 'elapsed',
  'progress', 'prescribedResidual', 'residualContext', 'dissipation', 'limiter'];
const rejectionKeys = ['stage', 'message', 'reason', 'code', 'diagnostics', 'step', 'maximumStep', 'iteration', 'backtracks', 'limiter'];

// Inspect only the current failure/cause chain, never earlier retry histories.
// Both thrown errors and returned precursor failures use these containers.
function localFailureSummary(error) {
  const lines = [], pending = [error], seen = new Set();
  if (present(error.code)) lines.push(`Failure code: ${error.code}`);
  while (pending.length) {
    const value = pending.shift();
    if (!value || typeof value !== 'object' || seen.has(value)) continue;
    seen.add(value);
    if (value.interfacePressure !== undefined) {
      const cell = value.cell ?? value.section;
      if (cell) lines.push(`Failing cell: i=${show(cell.i)}, group=${show(cell.group)}, tube=${show(cell.tube)} (solver indices)`);
      const pressure = value.interfacePressure;
      lines.push(`Interface pressure (solver units): ${pressure && typeof pressure === 'object'
        ? `lower ${show(pressure.lower)}, upper ${show(pressure.upper)}` : show(pressure)}`);
      return lines;
    }
    for (const key of ['lastRejectedStep', 'diagnostics', 'cause', 'precursor', 'failure', 'failureDiagnostic', 'stagnationDensityFallback'])
      if (value[key]) pending.push(value[key]);
  }
  return lines;
}

export function buildSolverErrorContext({ message, input, event = {}, progress = {}, mesh, label } = {}) {
  input ??= {}; event ??= {}; progress ??= {};
  const diagnostic = event.failureDiagnostic ?? {};
  const actualMach = firstNumber(event.actualMach, diagnostic.actualMach, event.conditions?.mach, event.mach,
    event.failure?.actualMach, progress.actualMach, progress.mach);
  const requestedMach = firstNumber(input.mach, event.targetMach, diagnostic.targetMach, progress.targetMach);
  const explicitStage = firstValue(event.stage, event.failure?.stage, diagnostic.stage);
  const stage = firstValue(explicitStage, progress.stage);
  const alphaContinuation = event.alphaContinuation;
  const actualAlpha = firstNumber(event.checkpoint?.restart?.input?.alpha, event.actualAlpha, alphaContinuation?.actualAlpha);
  const targetAlpha = firstNumber(input.alpha, event.targetAlpha, alphaContinuation?.targetAlpha);
  const targetNcrit = firstNumber(input.ncrit, event.targetNcrit, event.ncritContinuation?.targetNcrit, progress.targetNcrit);
  const sameStage = !present(explicitStage) || !present(progress.stage) || explicitStage === progress.stage;
  const actualNcrit = firstNumber(event.checkpoint?.restart?.options?.ncrit, event.conditions?.ncrit,
    event.actualNcrit, diagnostic.actualNcrit, event.ncritContinuation?.actualNcrit, event.failure?.actualNcrit,
    ...(sameStage ? [progress.actualNcrit] : []));
  const ncritContext = actualNcrit !== undefined || event.ncritContinuation !== undefined
    || event.actualNcrit !== undefined || progress.actualNcrit !== undefined;
  const explicitIteration = firstNumber(event.iteration, event.iteration?.iteration, event.iterations, event.diagnostics?.iterations);
  // An iteration number from a different stage is not the failed stage's
  // iteration. Keep that old observation only in report.progress.
  const iteration = explicitIteration ?? (!present(explicitStage) || !present(progress.stage) || explicitStage === progress.stage
    ? firstNumber(progress.iteration, progress.iteration?.iteration) : undefined);
  const sourceId = firstValue(typeof input.geometrySource === 'string' ? input.geometrySource : input.geometrySource?.id,
    input.presetId, input.airfoilId);
  const elements = Array.isArray(input.elements) ? input.elements : [];
  const airfoil = firstValue(label, input.name, elements.map((element, i) => element?.name ?? `Element ${i + 1}`).join(' + '), sourceId);
  const mode = input.flowModel === 'streamtube-grid' ? input.quadBoundaryLayers ? 'Quad Euler + boundary layers' : 'Quad Euler · inviscid'
    : input.flowModel === 'streamtube-bl' ? 'Quad Euler + boundary layers'
      : input.flowModel === 'coupled' ? 'Panel flow + boundary layers' : input.flowModel === 'inviscid' ? 'Panel flow · inviscid' : show(input.flowModel);
  const transition = input.transitionMode ?? (input.quadBoundaryLayers ? 'fixed-trip (default)' : 'not submitted');
  const trips = Array.isArray(input.materialTrips) ? input.materialTrips.map((pair, i) =>
    `${elements[i]?.name ?? `Element ${i + 1}`}: upper ${number(pair?.[0])}, lower ${number(pair?.[1])}`).join('; ')
    : Array.isArray(input.trips) ? `upper ${number(input.trips[0])}, lower ${number(input.trips[1])}` : 'not submitted';
  const smoothing = input.gridEllipticSmoothing === true ? `on (${input.gridSmoothingMethod ?? 'default'})`
    : input.gridEllipticSmoothing === false ? 'off' : 'not submitted';
  const currentMesh = event.mesh ?? mesh;
  const actualMesh = {
    cellCount: Number.isInteger(currentMesh?.cells?.length) ? currentMesh.cells.length : currentMesh?.cellCount,
    streamwiseIntervals: currentMesh?.initialization?.streamwiseSegments ?? currentMesh?.streamwiseSegments,
    tubes: plain(currentMesh?.initialization?.tubes ?? currentMesh?.tubes),
    quality: plain(currentMesh?.quality, { diagnostic: true }),
  };
  const error = { message: typeof message === 'string' ? message : String(message ?? event.message ?? event.reason ?? 'Unknown solver error'),
    ...selected(event, ['code', 'reason', 'solverStopReason', 'stack']) };
  if (present(stage)) error.stage = plain(stage);
  if (event.diagnostics !== undefined) error.diagnostics = plain(event.diagnostics, { diagnostic: true });
  if (event.failure !== undefined) error.failure = selected(event.failure, ['message', 'reason', 'code', 'stage', 'actualMach', 'targetMach', 'actualNcrit', 'targetNcrit', 'diagnostics']);
  if (event.ncritContinuation !== undefined) error.ncritContinuation = plain(event.ncritContinuation, { diagnostic: true });
  if (event.failureDiagnostic !== undefined) error.failureDiagnostic = plain(event.failureDiagnostic, { diagnostic: true });
  if (event.cause !== undefined) error.cause = errorCause(event.cause);
  if (event.precursor !== undefined) {
    error.precursor = selected(event.precursor, ['reason', 'code', 'stage', 'solverStopReason', 'quality']);
    if (event.precursor?.lastRejectedStep !== undefined)
      error.precursor.lastRejectedStep = selected(event.precursor.lastRejectedStep, rejectionKeys);
  }
  const rejected = event.lastRejectedStep ?? event.diagnostics?.lastRejectedStep ?? event.flow?.lastRejectedStep;
  if (rejected !== undefined) error.lastRejectedStep = selected(rejected, rejectionKeys);
  const report = plain({ schemaVersion: 1, label, input, airfoil: { name: airfoil, id: sourceId }, mode,
    requestedMach, actualMach, ...(event.operatingPointContinuation ? { operatingPointContinuation: event.operatingPointContinuation } : {}), ...(alphaContinuation ? { actualAlpha, targetAlpha, alphaContinuation } : {}), ...(ncritContext ? { actualNcrit, targetNcrit } : {}),
    stage, iteration, error, progress: selected(progress, progressKeys), mesh: actualMesh });
  const invalid = actualMesh.quality?.invalidCells;
  const lines = [
    `Airfoil: ${show(airfoil)}${present(sourceId) ? ` [${sourceId}]` : ''} · Mode: ${mode}`,
    `α ${number(input.alpha)}°${alphaContinuation ? ` · Actual α ${number(actualAlpha)}°` : ''} · Requested Mach ${number(requestedMach)} · Actual Mach ${number(actualMach)} · Re ${number(input.reynolds)} · ${ncritContext
      ? `Requested Ncrit ${number(targetNcrit)} · Actual Ncrit ${number(actualNcrit)}` : `Ncrit ${number(input.ncrit)}`}`,
    `Transition: ${transition} · Trips: ${trips} · ISMOM: ${automatic(input.eulerIsmom)}${input.eulerStartup ? ` · Euler startup: ${input.eulerStartup}` : ''}`,
    `Requested grid: ${number(input.gridIntervals)} surface intervals · ${number(input.gridTubes)} tubes · Inlet ${automatic(input.gridInletIntervals)} · Wake ${automatic(input.gridOutletIntervals)}`,
    `Tube overrides: upper ${automatic(input.gridUpperTubes)}, lower ${automatic(input.gridLowerTubes)}, gap ${automatic(input.gridGapTubes)} (tube counts)`,
    `SLOR: ${smoothing} · Spacing: ${show(input.gridSurfaceSpacing)} · Chord exponent: ${number(input.gridChordExponent)} · Stagnation aspect ratio: ${present(input.gridStagnationAspectRatio) ? number(input.gridStagnationAspectRatio) : 'off'}`,
  ];
  if (input.gridCurvatureSpacing) lines.push(`Curvature spacing: exponent ${number(input.gridCurvatureSpacing.exponent)}, LE ratio ${number(input.gridCurvatureSpacing.leadingSpacingRatio)}, TE ratio ${number(input.gridCurvatureSpacing.trailingSpacingRatio)}`);
  lines.push(`Actual mesh: ${number(actualMesh.cellCount)} cells · ${number(actualMesh.streamwiseIntervals)} streamwise intervals · tubes ${Array.isArray(actualMesh.tubes) ? `[${actualMesh.tubes.join(', ')}]` : 'unknown'} · invalid cells ${Array.isArray(invalid) ? `[${invalid.join(', ')}]` : 'unknown'}`);
  lines.push(`Stage: ${show(stage)} · Iteration: ${number(iteration)}`);
  lines.push(...localFailureSummary(error));
  return { summary: lines.join('\n'), report };
}
