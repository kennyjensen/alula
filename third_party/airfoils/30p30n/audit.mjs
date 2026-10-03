// Source-coordinate inspection only. Does not construct a mesh or evaluate flow.
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { prepareContour, signedArea } from '../../../src/geometry/airfoil.js';

const here = new URL('./', import.meta.url);
const source = readFileSync(new URL('30p-30n.dat', here));
const elements = [];
for (const line of source.toString('utf8').split(/\r?\n/)) {
  if (/^# (Slat|Main Element|Flap)$/.test(line)) elements.push({ name: line.slice(2), points: [] });
  else if (line.trim() && !line.startsWith('#')) {
    const parts = line.trim().split(/\s+/).map(Number);
    if (!elements.length || parts.length !== 2 || !parts.every(Number.isFinite)) throw new Error('Unexpected coordinate row.');
    elements.at(-1).points.push({ x: parts[0], y: parts[1] });
  }
}
if (elements.length !== 3) throw new Error('Expected three source blocks.');
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const angle = (a, b, c) => {
  const dot = (b.x - a.x) * (c.x - a.x) + (b.y - a.y) * (c.y - a.y);
  return Math.acos(Math.max(-1, Math.min(1, dot / (distance(a, b) * distance(a, c))))) * 180 / Math.PI;
};
const report = {
  passed: true,
  meaning: 'Successful source inspection; NOT a solver or aerodynamic acceptance result.',
  sourceSha256: createHash('sha256').update(source).digest('hex'),
  sourceBytes: source.length,
  geometryChanged: false,
  flowEvaluations: 0,
  referenceChord: 1,
  elements: elements.map(({ name, points }) => {
    const first = points[0], second = points[1], penultimate = points.at(-2), last = points.at(-1);
    let solverAccepted = true, solverRejection = null;
    try { prepareContour(points); } catch (error) { solverAccepted = false; solverRejection = error.message; }
    return {
      name, pointCount: points.length, finite: points.every(p => Number.isFinite(p.x) && Number.isFinite(p.y)),
      bounds: { x: [Math.min(...points.map(p => p.x)), Math.max(...points.map(p => p.x))], y: [Math.min(...points.map(p => p.y)), Math.max(...points.map(p => p.y))] },
      first, second, penultimate, last,
      signedPolygonAreaWithStraightClosure: signedArea(points),
      endpointGap: distance(first, last),
      firstSurfaceInterval: distance(first, second), lastSurfaceInterval: distance(penultimate, last),
      ...(distance(first, last) > 0 ? { straightBaseOpeningAngleAtFirstDegrees: angle(first, second, last), straightBaseOpeningAngleAtLastDegrees: angle(last, penultimate, first) } : {}),
      solverAccepted, solverRejection,
    };
  }),
};
if (report.elements.map(e => e.pointCount).join(',') !== '201,221,242'
  || report.elements.some(e => !e.finite)
  || report.elements[0].endpointGap !== 0 || report.elements[1].endpointGap !== 0
  || !(report.elements[2].endpointGap > .0058 && report.elements[2].endpointGap < .0059)
  || report.elements.map(e => e.solverAccepted).join(',') !== 'true,true,false') throw new Error('Source or solver assumptions changed.');
writeFileSync(new URL('geometry-audit.json', here), JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ passed: report.passed, sourceSha256: report.sourceSha256, flowEvaluations: 0, elements: report.elements.map(({ name, pointCount, endpointGap, solverAccepted }) => ({ name, pointCount, endpointGap, solverAccepted })) }, null, 2));
