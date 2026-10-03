// Exact analytic harmonic domains; no panel field or grid generator is used.
const copy = p => ({ x: p.x, y: p.y });
const scale = (p, a) => ({ x: a * p.x, y: a * p.y });
const validCount = n => Number.isInteger(n) && n >= 1;
const checked = (k, count) => {
  if (!Number.isFinite(k) || k < 0 || k > count) throw new Error('Boundary index is outside the fixture.');
};

export function conformalPolynomialFixture({ nx = 4, nt = 4 } = {}) {
  if (!validCount(nx) || !validCount(nt)) throw new Error('Invalid polynomial fixture dimensions.');
  const map = (u, v) => ({ x: u + .2 * (u * u - v * v), y: v + .4 * u * v });
  const du = (u, v) => ({ x: 1 + .4 * u, y: .4 * v });
  const dv = (u, v) => ({ x: -.4 * v, y: 1 + .4 * u });
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => map(i / nx, j / nt)));
  const geometry = { at(i, j, s, t) {
    const u = (i + s) / nx, v = (j + t) / nt;
    return { point: map(u, v), ds: scale(du(u, v), 1 / nx), dt: scale(dv(u, v), 1 / nt) };
  } };
  const horizontal = v => k => { checked(k, nx); const u = k / nx;
    return { point: map(u, v), derivative: scale(du(u, v), 1 / nx) }; };
  const vertical = u => k => { checked(k, nt); const v = k / nt;
    return { point: map(u, v), derivative: scale(dv(u, v), 1 / nt) }; };
  const boundary = { bottom: horizontal(0), top: horizontal(1), left: vertical(0), right: vertical(1) };
  // Invert z=w+.2w^2. Imag(sqrt(1+.8z))/.4 = Im(z)/Re(sqrt(1+.8z)).
  const psiAt = p => {
    const real = 1 + .8 * p.x, imag = .8 * p.y;
    return p.y / Math.sqrt((Math.hypot(real, imag) + real) / 2);
  };
  return { nx, nt, nodes, boundary, geometry, massFlows: Array(nt).fill(1 / nt), psiAt };
}

export function cylinderHarmonicFixture({ nx = 12, nt = 6, height = 1.4 } = {}) {
  if (!validCount(nx) || nx % 6 || !validCount(nt) || !(height > 0))
    throw new Error('Cylinder fixture needs nx divisible by six and positive height/counts.');
  const sqrtUpper = (x, y) => ({
    x: Math.sqrt(Math.max(0, (Math.hypot(x, y) + x) / 2)),
    y: Math.sqrt(Math.max(0, (Math.hypot(x, y) - x) / 2)),
  });
  const inverse = (phi, psi) => {
    const a = sqrtUpper(phi - 2, psi), b = sqrtUpper(phi + 2, psi);
    return { x: .5 * (phi + a.x * b.x - a.y * b.y), y: .5 * (psi + a.x * b.y + a.y * b.x) };
  };
  const inverseDerivative = p => {
    const r2 = p.x * p.x + p.y * p.y, r4 = r2 * r2;
    const real = 1 - (p.x * p.x - p.y * p.y) / r4, imag = 2 * p.x * p.y / r4;
    const d2 = real * real + imag * imag;
    if (!(d2 > 0)) throw new Error('The inverse-potential derivative is singular at stagnation.');
    return { x: real / d2, y: -imag / d2 };
  };
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => inverse(-3 + 6 * i / nx, height * j / nt)));
  const exactInteger = (k, point, pick) => Number.isInteger(k) ? copy(pick(k)) : point;
  const bottom = (k, { interval } = {}) => {
    checked(k, nx); const i = interval ?? Math.min(nx - 1, Math.floor(k)), s = k - i;
    if (!Number.isInteger(i) || i < 0 || i >= nx || s < 0 || s > 1) throw new Error('Invalid cylinder boundary interval.');
    const a = nodes[i][0], b = nodes[i + 1][0]; let point, derivative;
    if (i >= nx / 6 && i < 5 * nx / 6) {
      const theta0 = Math.atan2(a.y, a.x), delta = Math.atan2(b.y, b.x) - theta0, theta = theta0 + s * delta;
      point = { x: Math.cos(theta), y: Math.sin(theta) };
      derivative = { x: -Math.sin(theta) * delta, y: Math.cos(theta) * delta };
    } else {
      point = { x: a.x + s * (b.x - a.x), y: 0 }; derivative = { x: b.x - a.x, y: 0 };
    }
    return { point: exactInteger(k, point, i => nodes[i][0]), derivative };
  };
  const top = k => {
    checked(k, nx); const point = inverse(-3 + 6 * k / nx, height), d = inverseDerivative(point);
    return { point: exactInteger(k, point, i => nodes[i][nt]), derivative: scale(d, 6 / nx) };
  };
  const end = i => k => {
    checked(k, nt); const point = inverse(-3 + 6 * i / nx, height * k / nt), d = inverseDerivative(point);
    return { point: exactInteger(k, point, j => nodes[i][j]), derivative: { x: -d.y * height / nt, y: d.x * height / nt } };
  };
  const psiAt = p => p.y * (1 - 1 / (p.x * p.x + p.y * p.y));
  return { nx, nt, nodes, boundary: { bottom, top, left: end(0), right: end(nx) },
    massFlows: Array(nt).fill(height / nt), psiAt, height };
}
