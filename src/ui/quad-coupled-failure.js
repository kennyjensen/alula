// SPDX-License-Identifier: GPL-2.0-or-later
// Describe an already rejected grid. No geometry/gas/residual evaluation,
// changed acceptance flag, or inference of an unrecorded source revision.
export function quadCoupledDisplayedStateLabel(result) {
  if (result?.stateConverged === true && result.mesh?.quality?.valid === true) return 'retained converged state';
  return result?.failureDiagnostic?.stage === 'cold-coupled' ? 'failed startup state' : 'displayed state';
}

export function quadCoupledStartingFlowLabel(result, sourceMach) {
  if (!Number.isFinite(sourceMach)) return 'Source condition unavailable';
  const state = quadCoupledDisplayedStateLabel(result);
  const prefix = state === 'failed startup state' ? 'Failed startup at'
    : state === 'retained converged state' ? 'Converged' : 'Source';
  return `${prefix} Mach ${sourceMach.toFixed(3)}`;
}

export function quadCoupledGridFailure(raw, input = {}) {
  const mesh = raw?.mesh;
  if (mesh?.quality?.valid !== false) return null;
  const continuation = raw.machContinuation
    ?? (raw.continuation?.method === 'freestream-mach' ? raw.continuation : null);
  const actualMach = raw.conditions?.mach ?? raw.actualMach ?? raw.mach;
  const targetMach = continuation?.targetMach ?? raw.targetMach ?? input.mach;
  const explicitStage = raw.failure?.stage ?? continuation?.failureStage;
  const startup = ['cold-initialization', 'cold-coupled'].includes(explicitStage)
    || continuation?.coldBaselineUsed === true && raw.stateConverged !== true
      && actualMach === continuation.sourceMach && raw.initialization?.attempts?.length > 0
    || !continuation && !raw.refinement && raw.initialization?.attempts?.length > 0;
  const stage = startup ? 'cold-coupled' : raw.refinement ? 'refinement'
    : explicitStage ?? (continuation?.attempts?.length > 0 ? 'mach-continuation' : 'coupled');
  const stageLabel = stage === 'cold-coupled' ? 'Cold coupled startup'
    : stage === 'refinement' ? 'Coupled refinement'
      : stage === 'mach-continuation' ? 'Coupled Mach continuation'
        : stage === 'boundary-layer-initialization' ? 'Boundary-layer initialization' : 'Coupled solve';
  const ids = Array.isArray(mesh.quality.invalidCells) ? mesh.quality.invalidCells.slice() : [];
  const { tubes, streamwiseSegments: nx } = mesh.initialization ?? {};
  const bodies = raw.solverInput?.bodies ?? raw.bodies;
  const mapped = mesh.topology === 'intrinsic-quadrilateral-streamtubes'
    && Number.isInteger(nx) && nx > 0 && Array.isArray(tubes) && tubes.length > 0
    && tubes.every(n => Number.isInteger(n) && n > 0)
    && Array.isArray(mesh.cells) && mesh.cells.length === nx * tubes.reduce((a, b) => a + b, 0);
  const locations = ids.slice(0, 8).map(id => {
    if (!mapped || !Number.isInteger(id) || id < 0 || id >= mesh.cells.length) return { cell: id, resolved: false };
    let group = 0, local = id;
    while (local >= nx * tubes[group]) local -= nx * tubes[group++];
    const interval = Math.floor(local / tubes[group]), tube = local % tubes[group], boundaries = [];
    const add = (body, side) => {
      const b = bodies?.[body];
      if (!Array.isArray(bodies) || bodies.length !== tubes.length - 1 || !b
        || !Number.isInteger(b.leadingIndex) || !Number.isInteger(b.trailingIndex)
        || !(0 <= b.leadingIndex && b.leadingIndex < b.trailingIndex && b.trailingIndex < nx)) return;
      const element = Number.isInteger(b.element) ? b.element : null;
      const region = interval < b.leadingIndex ? 'upstream cut' : interval < b.trailingIndex ? 'surface' : 'wake';
      boundaries.push({ body, element, ...(element !== null && input.elements?.[element]?.name
        ? { name: input.elements[element].name } : {}), side, region,
        ...(region === 'wake' ? { wakeInterval: interval - b.trailingIndex } : {}) });
    };
    if (tube === 0 && group > 0) add(group - 1, 'upper');
    if (tube === tubes[group] - 1 && Array.isArray(bodies) && group < bodies.length) add(group, 'lower');
    return { cell: id, resolved: true, group, interval, tube, boundaries };
  });
  const first = locations[0];
  const locationText = !first ? '' : first.resolved
    ? ` First rejected cell ${first.cell}: passage ${first.group}, interval ${first.interval}, tube ${first.tube}`
      + (first.boundaries.length ? ` (${first.boundaries.map(b => `${b.name ?? (b.element === null ? `body ${b.body}` : `element ${b.element + 1}`)} ${b.side} ${b.region === 'wake' && b.wakeInterval === 0 ? 'first wake cell' : b.region}`).join('; ')})` : '') + '.'
    : ` First reported rejected cell: ${first.cell}; logical location unavailable.`;
  const machText = Number.isFinite(actualMach) ? ` at Mach ${actualMach.toFixed(3)}` : '';
  const countText = ids.length ? ` (${ids.length} rejected cell${ids.length === 1 ? '' : 's'})` : '';
  return {
    code: 'coupled-final-displacement-grid', stage,
    ...(explicitStage ? { reportedStage: explicitStage } : {}),
    message: `${stageLabel}${machText} failed: final displacement grid rejected${countText}.${locationText}`,
    solverReason: raw.reason, actualMach, targetMach,
    ...(raw.ncritContinuation || raw.actualNcrit !== undefined ? {
      actualNcrit: raw.checkpoint?.restart?.options?.ncrit ?? raw.conditions?.ncrit ?? raw.actualNcrit,
      targetNcrit: input.ncrit ?? raw.targetNcrit } : {}),
    invalidCells: ids, invalidCellCount: ids.length, locations, indicesAreZeroBased: true,
    minArea: mesh.quality.minArea, minCornerSine: mesh.quality.minCornerSine,
    returnedGridValid: false, returnedStateConverged: false,
    case: { alpha: input.alpha, reynolds: input.reynolds, ncrit: input.ncrit,
      transitionMode: input.transitionMode, materialTrips: input.materialTrips?.map(p => p.slice()),
      gridIntervals: input.gridIntervals, gridTubes: input.gridTubes,
      gridEllipticSmoothing: input.gridEllipticSmoothing, geometryId: input.geometrySource?.id },
    scope: 'Description of the returned rejected mesh; no geometry or residual recomputation. Source revision is not inferred from the application version.',
  };
}
