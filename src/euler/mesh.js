// SPDX-License-Identifier: GPL-2.0-or-later
const edgeKey = (a, b) => a < b ? `${a}:${b}` : `${b}:${a}`;

// Exact edge vectors, signed polygon areas and explicit owner/neighbor topology.
export function buildMesh(vertices, connectivity, boundaries = []) {
  if (!vertices.length || !connectivity.length || vertices.some(v => !Number.isFinite(v.x) || !Number.isFinite(v.y))) throw new Error('Invalid mesh vertices.');
  const faces = []; const byEdge = new Map(); const tags = new Map();
  for (const b of boundaries) {
    const key = edgeKey(b.a, b.b);
    if (tags.has(key)) throw new Error('Duplicate mesh boundary tag.');
    tags.set(key, b);
  }
  const cells = connectivity.map((indices, index) => {
    if (indices.length < 3 || new Set(indices).size !== indices.length || indices.some(i => !Number.isInteger(i) || i < 0 || i >= vertices.length)) throw new Error('Invalid cell connectivity.');
    const origin = vertices[indices[0]]; let area2 = 0; let cx = 0; let cy = 0; let perimeter = 0;
    const refs = [];
    for (let k = 0; k < indices.length; k++) {
      const a = indices[k]; const b = indices[(k + 1) % indices.length];
      const p = vertices[a]; const q = vertices[b];
      const cross = (p.x - origin.x) * (q.y - origin.y) - (q.x - origin.x) * (p.y - origin.y);
      area2 += cross; cx += (p.x + q.x - 2 * origin.x) * cross; cy += (p.y + q.y - 2 * origin.y) * cross;
      const dx = q.x - p.x; const dy = q.y - p.y; const length = Math.hypot(dx, dy); perimeter += length;
      if (!(length > 0)) throw new Error('Zero-length mesh face.');
      const next = vertices[indices[(k + 2) % indices.length]];
      if (dx * (next.y - q.y) - dy * (next.x - q.x) <= 1e-14 * length ** 2) throw new Error('Inverted or nonconvex mesh cell.');
      const key = edgeKey(a, b); const existing = byEdge.get(key);
      if (existing !== undefined) {
        const face = faces[existing];
        if (face.neighbor !== null || face.a !== b || face.b !== a) throw new Error('Nonmanifold or inconsistently oriented mesh face.');
        face.neighbor = index; refs.push({ face: existing, sign: -1 });
      } else {
        byEdge.set(key, faces.length); refs.push({ face: faces.length, sign: 1 });
        faces.push({ a, b, owner: index, neighbor: null, nx: dy / length, ny: -dx / length, length,
          x: (p.x + q.x) / 2, y: (p.y + q.y) / 2, boundary: null });
      }
    }
    if (!(area2 > 1e-14 * perimeter ** 2)) throw new Error('Nonpositive or unresolved mesh cell area.');
    return { vertices: [...indices], faces: refs, area: area2 / 2,
      x: origin.x + cx / (3 * area2), y: origin.y + cy / (3 * area2), perimeter };
  });
  for (const face of faces) {
    const key = edgeKey(face.a, face.b); const tag = tags.get(key);
    if (face.neighbor === null) {
      if (!tag?.type) throw new Error(`Untagged boundary face ${key}.`);
      face.boundary = { ...tag }; tags.delete(key);
    } else if (tag) throw new Error('Internal face was tagged as a boundary.');
  }
  if (tags.size) throw new Error('Boundary tags reference absent faces.');
  let metricClosure = 0;
  for (const cell of cells) {
    let x = 0; let y = 0;
    for (const ref of cell.faces) { const f = faces[ref.face]; x += ref.sign * f.nx * f.length; y += ref.sign * f.ny * f.length; }
    metricClosure = Math.max(metricClosure, Math.hypot(x, y) / cell.perimeter);
  }
  return { vertices, cells, faces, diagnostics: { minArea: cells.reduce((a,c)=>Math.min(a,c.area),Infinity), metricClosure,
    internalFaces: faces.filter(f => f.neighbor !== null).length, boundaryFaces: faces.filter(f => f.neighbor === null).length } };
}

export function channelMesh({ nx = 12, ny = 4, length = 2, lower = () => 0, upper = () => 1, map } = {}) {
  if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(ny) || ny < 1 || !Number.isFinite(length) || length <= 0) throw new Error('Invalid channel mesh dimensions.');
  const vertices = []; const cells = []; const boundaries = [];
  const id = (i, j) => i * (ny + 1) + j;
  for (let i = 0; i <= nx; i++) {
    const x = length * i / nx; const lo = lower(x); const hi = upper(x);
    if (!(hi > lo)) throw new Error('Closed or inverted channel.');
    for (let j = 0; j <= ny; j++) {
      const y = lo + j / ny * (hi - lo);
      vertices.push(map ? map({ x, y, i, j, nx, ny }) : { x, y });
    }
  }
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) cells.push([id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]);
  for (let i = 0; i < nx; i++) {
    boundaries.push({ a: id(i, 0), b: id(i + 1, 0), type: 'wall', side: 'lower', station: i });
    boundaries.push({ a: id(i, ny), b: id(i + 1, ny), type: 'wall', side: 'upper', station: i });
  }
  for (let j = 0; j < ny; j++) {
    boundaries.push({ a: id(0, j), b: id(0, j + 1), type: 'inlet' });
    boundaries.push({ a: id(nx, j), b: id(nx, j + 1), type: 'outlet' });
  }
  const streamlineNodes = [];
  if (!map) for (let i = 1; i <= nx; i++) for (let j = 1; j < ny; j++) {
    streamlineNodes.push({ vertex: id(i, j), previous: id(i - 1, j), scale: (upper(length * i / nx) - lower(length * i / nx)) / ny });
  }
  return { ...buildMesh(vertices, cells, boundaries), structured: { nx, ny, length }, streamlineNodes };
}
