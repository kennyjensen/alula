import { cylinderHarmonicFixture } from './curved-harmonic-fields.js';

const point = p => ({ x: p.x, y: p.y });
const linear = (a, b, t) => ({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
const chord = (a, b) => [point(a), linear(a, b, 1 / 3), linear(a, b, 2 / 3), point(b)];

export function cubicCylinderGridFixture(options = {}) {
  const data = cylinderHarmonicFixture(options), { nodes, nx, nt } = data;
  const hermite = (side, k, a, b) => {
    const da = data.boundary[side](k, { interval: k }).derivative;
    const db = data.boundary[side](k + 1, { interval: k }).derivative;
    return [point(a), { x: a.x + da.x / 3, y: a.y + da.y / 3 },
      { x: b.x - db.x / 3, y: b.y - db.y / 3 }, point(b)];
  };
  const horizontal = Array.from({ length: nx }, (_, i) => Array.from({ length: nt + 1 }, (_, j) =>
    !j ? hermite('bottom', i, nodes[i][j], nodes[i + 1][j])
      : j === nt ? hermite('top', i, nodes[i][j], nodes[i + 1][j]) : chord(nodes[i][j], nodes[i + 1][j])));
  const vertical = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt }, (_, j) =>
    !i ? hermite('left', j, nodes[i][j], nodes[i][j + 1])
      : i === nx ? hermite('right', j, nodes[i][j], nodes[i][j + 1]) : chord(nodes[i][j], nodes[i][j + 1])));
  const controlPoints = Array.from({ length: nx }, (_, i) => Array.from({ length: nt }, (_, j) => {
    const p = [nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]];
    const bottom = horizontal[i][j], top = horizontal[i][j + 1], left = vertical[i][j], right = vertical[i + 1][j];
    return Array.from({ length: 4 }, (_, a) => Array.from({ length: 4 }, (_, b) => {
      // Entire faces come from canonical arrays; no rounded reconstruction
      // can give adjacent cells different shared-face coefficients.
      if (!a) return point(left[b]); if (a === 3) return point(right[b]);
      if (!b) return point(bottom[a]); if (b === 3) return point(top[a]);
      const s = a / 3, t = b / 3, q = linear(linear(p[0], p[1], s), linear(p[3], p[2], s), t);
      for (const [edge, base, weight] of [[bottom[a], linear(p[0], p[1], s), 1 - t],
        [top[a], linear(p[3], p[2], s), t], [left[b], linear(p[0], p[3], t), 1 - s],
        [right[b], linear(p[1], p[2], t), s]]) {
        q.x += weight * (edge.x - base.x); q.y += weight * (edge.y - base.y);
      }
      return q;
    }));
  }));
  const tangent = (a, b) => { const x = b.x - a.x, y = b.y - a.y, d = Math.hypot(x, y);
    if (!(d > 0)) throw new Error('A cylinder streamwise edge is degenerate.'); return { x: x / d, y: y / d }; };
  const directions = nodes.map((row, i) => row.map((p, j) => {
    const a = i ? tangent(nodes[i - 1][j], p) : { x: 0, y: 0 };
    const b = i < nx ? tangent(p, nodes[i + 1][j]) : { x: 0, y: 0 };
    const x = a.x + b.x, y = a.y + b.y, d = Math.hypot(x, y);
    if (!(d > 1e-12)) throw new Error('A cylinder transverse guide is unresolved.');
    return { x: -y / d, y: x / d };
  }));
  return { ...data, controlPoints, directions,
    boundaryApproximation: 'Cubic Hermite endpoint interpolation of analytic circle/cut/potential curves. Circle arcs are polynomial approximations, not exact circles.' };
}
