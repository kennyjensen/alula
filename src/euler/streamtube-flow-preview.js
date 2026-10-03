// SPDX-License-Identifier: GPL-2.0-or-later
// Display the actual discrete Euler velocity field: a tube's direction joins
// successive bank midpoints, and its speed is the section state used by mass
// and momentum conservation. q is nondimensionalized by U_infinity = 1.
// These are piecewise-linear tube centerlines, not a separately traced field.
export function streamtubeFlowSnapshot({ nodes, sections, diagnostics }, iteration = 0) {
  let minimumSpeedRatio = Infinity, maximumSpeedRatio = -Infinity;
  const lines = nodes.flatMap((region, group) => Array.from({ length: region[0].length - 1 }, (_, tube) => {
    const points = region.map(row => ({ x: (row[tube].x + row[tube + 1].x) / 2,
      y: (row[tube].y + row[tube + 1].y) / 2 }));
    const speedRatios = sections.map(section => {
      const speed = section[group][tube].q;
      minimumSpeedRatio = Math.min(minimumSpeedRatio, speed);
      maximumSpeedRatio = Math.max(maximumSpeedRatio, speed);
      return speed;
    });
    return { group, tube, points, speedRatios,
      machNumbers: sections.map(section => { const m2 = section[group][tube].machSquared;
        return Number.isFinite(m2) && m2 >= 0 ? Math.sqrt(m2) : null; }) };
  }));
  return { kind: 'euler-tube-centerlines', speedUnit: 'U_infinity', iteration,
    residual: diagnostics.residual, minimumSpeedRatio, maximumSpeedRatio, lines };
}

// Compare corresponding tube sections on the evolving grid. This is a
// discrete-iterate comparison, not an Eulerian difference at fixed x,y.
export function compareStreamtubeFlow(current, initial, previous = initial) {
  if ([initial, previous].some(reference => reference.lines.length !== current.lines.length
    || reference.lines.some((line, k) => line.group !== current.lines[k].group || line.tube !== current.lines[k].tube
      || line.speedRatios.length !== current.lines[k].speedRatios.length)))
    throw new Error('Flow comparison requires corresponding streamtube sections.');
  let maximumSpeedChangeFromPrevious = 0, minimumSpeedChange = Infinity, maximumSpeedChange = -Infinity;
  const lines = current.lines.map((line, k) => ({ ...line, speedChanges: line.speedRatios.map((q, i) => {
    const delta = q - initial.lines[k].speedRatios[i];
    minimumSpeedChange = Math.min(minimumSpeedChange, delta); maximumSpeedChange = Math.max(maximumSpeedChange, delta);
    maximumSpeedChangeFromPrevious = Math.max(maximumSpeedChangeFromPrevious, Math.abs(q - previous.lines[k].speedRatios[i]));
    return delta;
  }) }));
  return { ...current, lines, changeReference: 'initial corresponding tube section',
    minimumSpeedChange, maximumSpeedChange, maximumSpeedChangeFromPrevious };
}

// Existing flow data only: wall pressure plus the final exit section needed
// for provisional loads. No residual evaluation or numerical state mutation.
export function observableFlow(flow) {
  const sections = flow.sections ? Array(flow.sections.length).fill(null) : undefined;
  if (sections?.length) sections[sections.length - 1] = flow.sections.at(-1)
    .map(group => group.map(({ rho, q, p, enthalpy }) => ({ rho, q, p, enthalpy })));
  return { nodes: flow.nodes, undisplacedNodes: flow.undisplacedNodes,
    cells: flow.cells.map(row => row.map(group => group.map(cell => ({ interfacePressure: cell.interfacePressure })))),
    ...(sections ? { sections } : {}),
    ...(flow.allocation?.groups ? { allocation: { groups: flow.allocation.groups.map(group => group.map(({ massFlow }) => ({ massFlow }))) } } : {}) };
}
export function observableBL(bl) {
  return { scale: bl.scale,
    stations: bl.stations.map(({ id, kind, body, side, i }) => ({ id, kind, body, side, i })),
    surfaces: bl.surfaces.map(({ body, side, ids }) => ({ body, side, ids })),
    wakes: bl.wakes.map(({ body, ids }) => ({ body, ids })) };
}
