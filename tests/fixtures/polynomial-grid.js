// Exact tensor-Q2 Bernstein representation of z=w+.2*w^2. Adjacent cells
// take their face coefficients from one canonical control lattice.
export function conformalPolynomialGridFixture({ nx = 4, nt = 3 } = {}) {
  if (![nx, nt].every(n => Number.isInteger(n) && n >= 1)) throw new Error('Invalid polynomial fixture dimensions.');
  const line = count => {
    const nodes = Array.from({ length: count + 1 }, (_, i) => i / count), linear = [], quadratic = [];
    for (let k = 0; k <= 2 * count; k++) {
      const i = Math.floor(k / 2);
      linear.push(k % 2 ? .5 * (nodes[i] + nodes[i + 1]) : nodes[i]);
      // For u(s)=u0+(u1-u0)s, Bernstein coefficients of u(s)^2
      // are [u0^2, u0*u1, u1^2], not samples at s=0,1/2,1.
      quadratic.push(k % 2 ? nodes[i] * nodes[i + 1] : nodes[i] * nodes[i]);
    }
    return { linear, quadratic };
  };
  const u = line(nx), v = line(nt);
  const lattice = u.linear.map((a, k) => v.linear.map((b, l) => ({
    x: a + .2 * (u.quadratic[k] - v.quadratic[l]), y: b + .4 * a * b,
  })));
  const controlPoints = Array.from({ length: nx }, (_, i) => Array.from({ length: nt }, (_, j) =>
    Array.from({ length: 3 }, (_, a) => Array.from({ length: 3 }, (_, b) => ({ ...lattice[2 * i + a][2 * j + b] })))));
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => ({ ...lattice[2 * i][2 * j] })));
  const analyticAt = (i, j, s, t) => {
    const u = (i + s) / nx, v = (j + t) / nt;
    return { point: { x: u + .2 * (u * u - v * v), y: v + .4 * u * v },
      ds: { x: (1 + .4 * u) / nx, y: .4 * v / nx },
      dt: { x: -.4 * v / nt, y: (1 + .4 * u) / nt } };
  };
  const psiAt = p => {
    const a = 1 + .8 * p.x, b = .8 * p.y;
    return p.y / Math.sqrt((Math.hypot(a, b) + a) / 2);
  };
  return { nx, nt, nodes, controlPoints, analyticAt, psiAt, massFlows: Array(nt).fill(1 / nt),
    directions: nodes.map(row => row.map(() => ({ x: 0, y: 1 }))) };
}
