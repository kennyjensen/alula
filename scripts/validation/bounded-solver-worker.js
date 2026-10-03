// SPDX-License-Identifier: GPL-2.0-or-later
import { parentPort, workerData } from 'node:worker_threads';
if (workerData.stage === 'geometry') {
  try {
    const { validateAssembly } = await import('../../src/geometry/airfoil.js');
    const { createContourTopology } = await import('../../src/geometry/contour-topology.js');
    const topologies = workerData.caseData.elements.map(e => createContourTopology(e.points, e));
    validateAssembly(topologies.map(t => t.points));
    const wettedPanels = topologies.reduce((n, t) => n + t.surface.panels.length, 0);
    if (wettedPanels > 700) throw new Error(`Assembly exceeds the public 700 wetted panel budget (${wettedPanels}).`);
    parentPort.postMessage({ type: 'geometry-ready', geometry: { elements: topologies.length, wettedPanels,
      solidPanels: topologies.reduce((n, t) => n + t.points.length - 1, 0),
      trailingEdges: topologies.map(t => t.kind) } });
  } catch (error) { parentPort.postMessage({ type: 'error', message: error.message, code: error.code, stage: 'geometry' }); }
} else {
  // Run the same message handler as the browser, including its normal startup
  // and Mach continuation. The parent watchdog can interrupt synchronous JS.
  globalThis.self = { postMessage: message => parentPort.postMessage(message) };
  await import('../../src/worker/solver.js');
  await self.onmessage({ data: { id: 1, caseData: workerData.caseData, meshOnly: workerData.stage === 'mesh',
    bounds: { xMin: -2, xMax: 4, yMin: -2, yMax: 2 } } });
}
