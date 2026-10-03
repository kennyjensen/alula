// SPDX-License-Identifier: GPL-2.0-or-later
// Historical sharp-TE geometry remains available for frozen numerical tests.
// User-facing generation uses naca4Standard below.
export function naca4(code = '0012', panels = 120, { trailingEdge = 'sharp' } = {}) {
  if (!['sharp', 'finite'].includes(trailingEdge)) throw new Error('Invalid NACA trailing-edge convention.');
  if (!/^\d{4}$/.test(code)) throw new Error('Use a four-digit NACA designation, such as 2412.');
  if (!Number.isInteger(panels) || panels < 20 || panels > 600 || panels % 2) {
    throw new Error('Panel count must be even and between 20 and 600.');
  }
  const m = Number(code[0]) / 100;
  const p = Number(code[1]) / 10;
  const thickness = Number(code.slice(2)) / 100;
  if (!(thickness > 0) || (m > 0 && p === 0)) throw new Error('Invalid NACA camber or thickness.');
  const upper = []; const lower = [];
  for (let i = 0; i <= panels / 2; i++) {
    const x = (1 - Math.cos(2 * Math.PI * i / panels)) / 2;
    const yt = trailingEdge === 'sharp' && i === panels / 2 ? 0 : 5 * thickness * (0.2969 * Math.sqrt(x)
      - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - (trailingEdge === 'finite' ? 0.1015 : 0.1036) * x ** 4);
    let yc = 0; let slope = 0;
    if (m > 0) {
      const d = x < p ? p ** 2 : (1 - p) ** 2;
      yc = m / d * (x < p ? 2 * p * x - x * x : 1 - 2 * p + 2 * p * x - x * x);
      slope = 2 * m / d * (p - x);
    }
    // XFOIL NACA4 adds thickness vertically to the mean line (naca.f).
    const angle = trailingEdge === 'finite' ? 0 : Math.atan(slope);
    upper.push({ x: x - yt * Math.sin(angle), y: yc + yt * Math.cos(angle) });
    lower.push({ x: x + yt * Math.sin(angle), y: yc - yt * Math.cos(angle) });
  }
  return upper.reverse().concat(lower.slice(1));
}

// Standard XFOIL thickness polynomial; separate upper/lower TE endpoints.
export function naca4Standard(code = '0012', panels = 120) {
  return naca4(code, panels, { trailingEdge: 'finite' });
}

export function transform(points, { chord = 1, angle = 0, x = 0, y = 0 } = {}) {
  if (![chord, angle, x, y].every(Number.isFinite) || chord <= 0) throw new Error('Invalid element transform.');
  const a = angle * Math.PI / 180;
  return points.map(point => ({
    x: x + chord * (point.x * Math.cos(a) - point.y * Math.sin(a)),
    y: y + chord * (point.x * Math.sin(a) + point.y * Math.cos(a)),
  }));
}

export function signedArea(points) {
  // Relative coordinates avoid cancellation under translation.
  const origin = points[0];
  return points.reduce((sum, p, i) => {
    const q = points[(i + 1) % points.length];
    return sum + (p.x - origin.x) * (q.y - origin.y) - (q.x - origin.x) * (p.y - origin.y);
  }, 0) / 2;
}

export function pointInside(point, polygon) {
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const a = polygon[i]; const b = polygon[j];
    if ((a.y > point.y) !== (b.y > point.y)
      && point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x) inside = !inside;
  }
  return inside;
}

export function segmentsTouch(a, b, c, d, epsilon) {
  if (Math.max(a.x, b.x) < Math.min(c.x, d.x) - epsilon
    || Math.max(c.x, d.x) < Math.min(a.x, b.x) - epsilon
    || Math.max(a.y, b.y) < Math.min(c.y, d.y) - epsilon
    || Math.max(c.y, d.y) < Math.min(a.y, b.y) - epsilon) return false;
  const side = (p, q, r) => {
    const cross = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
    const threshold = epsilon * Math.hypot(q.x - p.x, q.y - p.y);
    return Math.abs(cross) <= threshold ? 0 : Math.sign(cross);
  };
  return side(a, b, c) * side(a, b, d) <= 0 && side(c, d, a) * side(c, d, b) <= 0;
}

export function prepareContour(input, { lifting = true } = {}) {
  if (!Array.isArray(input) || input.length < 9 || input.length > 601) {
    throw new Error('Each contour needs 8–600 panels and a repeated closing point.');
  }
  let points = input.map(p => ({ x: p.x, y: p.y }));
  if (!points.every(p => Number.isFinite(p.x) && Number.isFinite(p.y))) throw new Error('Nonfinite geometry coordinate.');
  const scale = Math.max(...points.map(p => Math.hypot(p.x - points[0].x, p.y - points[0].y)));
  if (!(scale > 0)) throw new Error('Degenerate contour.');
  const epsilon = scale * 1e-9;
  if (Math.hypot(points[0].x - points.at(-1).x, points[0].y - points.at(-1).y) > epsilon) {
    throw new Error('This baseline requires a closed, sharp trailing edge. Blunt/open contours are not supported yet.');
  }
  points.pop();
  if (Math.abs(signedArea(points)) < scale * scale * 1e-8) throw new Error('Contour has zero or unresolved area.');
  if (signedArea(points) < 0) points = [points[0], ...points.slice(1).reverse()];
  points.push({ ...points[0] });
  const n = points.length - 1;
  for (let i = 0; i < n; i++) {
    if (Math.hypot(points[i + 1].x - points[i].x, points[i + 1].y - points[i].y) < epsilon) {
      throw new Error('Contour contains a duplicate point or unresolved panel.');
    }
    for (let j = i + 2; j < n; j++) {
      if (i === 0 && j === n - 1) continue;
      if (segmentsTouch(points[i], points[i + 1], points[j], points[j + 1], epsilon)) {
        throw new Error('Contour intersects itself.');
      }
    }
  }
  if (lifting) {
    const a = points[1]; const b = points[n - 1]; const te = points[0];
    const cosine = ((a.x - te.x) * (b.x - te.x) + (a.y - te.y) * (b.y - te.y))
      / (Math.hypot(a.x - te.x, a.y - te.y) * Math.hypot(b.x - te.x, b.y - te.y));
    if (cosine < 0.5) throw new Error('Start and end the contour at its sharp trailing edge (opening angle < 60°).');
  }
  return points;
}

export function validateAssembly(contours) {
  for (let i = 0; i < contours.length; i++) {
    for (let j = i + 1; j < contours.length; j++) {
      const a = contours[i]; const b = contours[j];
      const scale = Math.max(...[a, b].flatMap(c => c.map(p => Math.hypot(p.x - c[0].x, p.y - c[0].y))));
      for (let k = 0; k < a.length - 1; k++) {
        for (let l = 0; l < b.length - 1; l++) {
          if (segmentsTouch(a[k], a[k + 1], b[l], b[l + 1], scale * 1e-9)) {
            throw new Error(`Elements ${i + 1} and ${j + 1} intersect or touch. Increase their gap.`);
          }
        }
      }
      if (pointInside(a[0], b) || pointInside(b[0], a)) throw new Error('One element is inside another.');
    }
  }
}
