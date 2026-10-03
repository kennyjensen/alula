// SPDX-License-Identifier: GPL-2.0-or-later
// Preserve the intrinsic grid's shared cuts and LE/TE endpoints when
// converting its logical groups to a physical quadrilateral mesh.
function physicalNodeKey(layout, g, i, j) {
  const node = layout.nodes[g][i][j], body = node.kind === 'wall' ? layout.bodies[node.body] : null;
  const endpoint = body && (i === body.leadingIndex || (!layout.displacedBoundaries && body.trailingEdge?.kind !== 'finite-base' && i === body.trailingIndex));
  const wakeBank = layout.displacedBoundaries && node.kind === 'cut' && i > layout.bodies[node.body].trailingIndex;
  // Column ids are already unique integers. Numeric Map keys avoid making
  // a new string for every interior point; boundary keys remain disjoint.
  return wakeBank ? `wake:${node.body}:${g}:${i}` : node.column !== null ? node.column : endpoint ? `end:${node.body}:${i}` : `wall:${node.body}:${node.side}:${i}`;
}

// The acceptance path needs the same shared-node check as a published mesh,
// but does not need vertex copies, cell indices or fixed-node sets.
export function validateStreamtubeMeshConnectivity(layout, nodes) {
  const points = new Map();
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length; i++) for (let j = 0; j < nodes[g][i].length; j++) {
    const p = nodes[g][i][j], key = physicalNodeKey(layout, g, i, j);
    if (!points.has(key)) points.set(key, p);
    const previous = points.get(key);
    if (previous.x !== p.x || previous.y !== p.y) throw new Error('Disconnected physical nodes on a shared streamtube cut or body endpoint.');
  }
}

export function streamtubeMeshConnectivity(layout, nodes) {
  const vertices = [], cells = [], ids = new Map(), fixed = new Set();
  const indices = nodes.map((group, g) => group.map((row, i) => row.map((p, j) => {
    const node = layout.nodes[g][i][j];
    const wakeBank = layout.displacedBoundaries && node.kind === 'cut' && i > layout.bodies[node.body].trailingIndex;
    const key = physicalNodeKey(layout, g, i, j);
    if (!ids.has(key)) { ids.set(key, vertices.length); vertices.push({ x: p.x, y: p.y }); }
    const id = ids.get(key), previous = vertices[id];
    if (previous.x !== p.x || previous.y !== p.y) throw new Error('Disconnected physical nodes on a shared streamtube cut or body endpoint.');
    // Geometry-only repair cannot move two wake banks independently while
    // preserving their prescribed gap. Euler still moves their shared center.
    if (node.kind === 'wall' || wakeBank || node.kind === 'farfield' || i === 0 || i === layout.nx) fixed.add(id);
    return id;
  })));
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < layout.nx; i++) for (let j = 0; j < layout.tubes[g]; j++)
    cells.push([indices[g][i][j], indices[g][i + 1][j], indices[g][i + 1][j + 1], indices[g][i][j + 1]]);
  return { vertices, cells, indices, fixed };
}
