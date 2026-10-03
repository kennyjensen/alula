// SPDX-License-Identifier: GPL-2.0-or-later
// Independent saved-state arithmetic only. No numerical-runtime imports.
import { directUpwindChannelConservation } from './upwind-streamtube.js';

const zero = () => [0, 0, 0, 0];
const plus = (a, b) => a.map((v, k) => v + b[k]);
const minus = (a, b) => a.map((v, k) => v - b[k]);
const midpoint = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
const normal = (a, b) => ({ x: b.y - a.y, y: a.x - b.x });
const pressureFlux = (n, p) => [0, p * n.x, p * n.y, 0];
const sum = rows => rows.reduce(plus, zero());
const maximum = rows => Math.max(0, ...rows.map(Math.abs));
const require = (condition, message) => { if (!condition) throw new Error(message); };

export function directDisplacedAssemblyLedger({ nodes, sections, cells, massFlows, bodies,
  stagnationEnthalpy, gamma = 1.4 }) {
  require(Array.isArray(nodes) && nodes.length === bodies.length + 1, 'One passage per body interval is required.');
  const nx = nodes[0].length - 1, tubes = nodes.map(group => group[0].length - 1);
  require(nx >= 3 && sections.length === nx && cells.length === nx - 1, 'Incomplete midpoint control volumes.');
  const groups = nodes.map((group, g) => {
    require(group.length === nx + 1 && group.every(row => row.length === tubes[g] + 1), 'Inconsistent passage nodes.');
    const value = directUpwindChannelConservation({ nodes: group, sections: sections.map(row => row[g]),
      cells: cells.map(row => row[g]), massFlows: massFlows[g] }, stagnationEnthalpy, gamma);
    const { local, ...summary } = value;
    return { group: g, tubes: tubes[g], volumeCount: local.length, ...summary };
  });
  const categories = {}, pairs = new Map(), internal = zero();
  let internalFaces = 0, halfFaces = 0;
  const bucket = name => categories[name] ??= { faces: 0, normal: { x: 0, y: 0 },
    pressure: zero(), physicalAdvection: zero(), transportAdvection: zero() };
  const tally = (name, n, p, physical = zero(), transport = physical) => {
    const b = bucket(name); b.faces++; b.normal.x += n.x; b.normal.y += n.y;
    b.pressure = plus(b.pressure, pressureFlux(n, p));
    b.physicalAdvection = plus(b.physicalAdvection, physical);
    b.transportAdvection = plus(b.transportAdvection, transport);
  };
  for (let g = 0; g < nodes.length; g++) for (let i = 1; i < nx; i++) {
    for (let j = 1; j < tubes[g]; j++) {
      const a = midpoint(nodes[g][i - 1][j], nodes[g][i][j]);
      const b = midpoint(nodes[g][i][j], nodes[g][i + 1][j]);
      const n = normal(b, a);
      const dp = cells[i - 1][g][j - 1].interfacePressure.upper - cells[i - 1][g][j].interfacePressure.lower;
      const flux = pressureFlux(n, dp);
      flux.forEach((v, k) => { internal[k] += v; }); internalFaces++;
    }
    for (const side of ['lower', 'upper']) {
      const lower = side === 'lower', j = lower ? 0 : tubes[g];
      const p = cells[i - 1][g][lower ? 0 : tubes[g] - 1].interfacePressure[side];
      const bodyIndex = lower ? g - 1 : g, body = bodies[bodyIndex];
      const halves = [[midpoint(nodes[g][i - 1][j], nodes[g][i][j]), nodes[g][i][j]],
        [nodes[g][i][j], midpoint(nodes[g][i][j], nodes[g][i + 1][j])]];
      halves.forEach((half, h) => {
        const [a, b] = lower ? half : half.slice().reverse(), n = normal(a, b), interval = i - 1 + h;
        const kind = !body ? (lower ? 'outerLower' : 'outerUpper')
          : interval < body.leadingIndex ? 'upstreamCut'
          : interval < body.trailingIndex ? 'displacedWall' : 'wakeBank';
        tally(kind, n, p); halfFaces++;
        if (kind === 'upstreamCut' || kind === 'wakeBank') {
          const key = `${bodyIndex}:${i}:${h}`, entry = pairs.get(key) ?? [];
          entry.push({ body: bodyIndex, element: body.element, kind, group: g, side, a, b, n, p }); pairs.set(key, entry);
        }
      });
    }
  }
  for (let g = 0; g < nodes.length; g++) for (let j = 0; j < tubes[g]; j++) for (const inlet of [true, false]) {
    const i = inlet ? 0 : nx - 1, sign = inlet ? -1 : 1, s = sections[i][g][j];
    const a = midpoint(nodes[g][i][j], nodes[g][i + 1][j]);
    const b = midpoint(nodes[g][i][j + 1], nodes[g][i + 1][j + 1]);
    const n = normal(a, b); n.x *= sign; n.y *= sign;
    const c = midpoint(nodes[g][i][j], nodes[g][i][j + 1]);
    const d = midpoint(nodes[g][i + 1][j], nodes[g][i + 1][j + 1]);
    const length = Math.hypot(d.x - c.x, d.y - c.y), direction = { x: (d.x - c.x) / length, y: (d.y - c.y) / length };
    const m = sign * massFlows[g][j], q = cells[inlet ? 0 : nx - 2][g][j].transportSpeeds[inlet ? 0 : 1];
    const flux = speed => [m, m * speed * direction.x, m * speed * direction.y, m * stagnationEnthalpy];
    tally(inlet ? 'inlet' : 'outlet', n, s.p, flux(s.q), flux(q));
  }
  const cutPairs = { paired: 0, coincident: 0, noncoincident: 0, unpaired: [],
    maximumEndpointSeparation: 0, maximumNormalMismatch: 0, maximumPressureMismatch: 0,
    retainedCoincidentFlux: zero(), retainedNoncoincidentFlux: zero(), perBody: [], largestNoncoincident: [] };
  const byBody = new Map();
  for (const [key, pair] of pairs) {
    if (pair.length !== 2) { cutPairs.unpaired.push({ key, faces: pair.length }); continue; }
    cutPairs.paired++;
    const [a, b] = pair, error = Math.max(Math.hypot(a.a.x - b.b.x, a.a.y - b.b.y), Math.hypot(a.b.x - b.a.x, a.b.y - b.a.y));
    const scale = Math.max(1, ...[a.a, a.b, b.a, b.b].flatMap(p => [Math.abs(p.x), Math.abs(p.y)]));
    const coincident = error <= 256 * Number.EPSILON * scale;
    const flux = plus(pressureFlux(a.n, a.p), pressureFlux(b.n, b.p));
    cutPairs[coincident ? 'coincident' : 'noncoincident']++;
    const name = coincident ? 'retainedCoincidentFlux' : 'retainedNoncoincidentFlux'; cutPairs[name] = plus(cutPairs[name], flux);
    cutPairs.maximumEndpointSeparation = Math.max(cutPairs.maximumEndpointSeparation, error);
    cutPairs.maximumNormalMismatch = Math.max(cutPairs.maximumNormalMismatch, Math.hypot(a.n.x + b.n.x, a.n.y + b.n.y));
    cutPairs.maximumPressureMismatch = Math.max(cutPairs.maximumPressureMismatch, Math.abs(a.p - b.p));
    const bodyKey = `${a.body}:${a.kind}`, body = byBody.get(bodyKey) ?? { body: a.body, element: a.element,
      kind: a.kind, pairs: 0, coincident: 0, noncoincident: 0, pressure: zero(), lowerBank: zero(), upperBank: zero() };
    body.pairs++; body[coincident ? 'coincident' : 'noncoincident']++; body.pressure = plus(body.pressure, flux);
    for (const face of pair) {
      const key = face.side === 'upper' ? 'lowerBank' : 'upperBank'; body[key] = plus(body[key], pressureFlux(face.n, face.p));
    }
    byBody.set(bodyKey, body);
    if (!coincident) cutPairs.largestNoncoincident.push({ key, kind: a.kind, body: a.body, error, pressureFlux: flux });
  }
  cutPairs.perBody = [...byBody.values()];
  cutPairs.largestNoncoincident.sort((a, b) => maximum(b.pressureFlux) - maximum(a.pressureFlux));
  cutPairs.largestNoncoincident = cutPairs.largestNoncoincident.slice(0, 8);
  for (const b of Object.values(categories)) {
    b.physical = plus(b.pressure, b.physicalAdvection); b.transport = plus(b.pressure, b.transportAdvection);
  }
  const whole = { localPhysical: sum(groups.map(g => g.physical.total)), localTransport: sum(groups.map(g => g.biased.total)),
    externalPhysical: sum(groups.map(g => g.physical.external)), externalTransport: sum(groups.map(g => g.biased.external)),
    internalPressureMismatch: internal, internalFaces };
  const category = names => sum(names.map(name => categories[name]?.physical ?? zero()));
  const remote = category(['inlet', 'outlet', 'outerLower', 'outerUpper']);
  const cuts = category(['upstreamCut', 'wakeBank']), wall = category(['displacedWall']);
  const recoveredWall = minus(minus(minus(whole.localPhysical, remote), cuts), internal);
  const errors = {
    classifiedPhysicalExterior: minus(sum(Object.values(categories).map(b => b.physical)), whole.externalPhysical),
    classifiedTransportExterior: minus(sum(Object.values(categories).map(b => b.transport)), whole.externalTransport),
    physicalInternalCancellation: minus(minus(whole.localPhysical, whole.externalPhysical), internal),
    transportInternalCancellation: minus(minus(whole.localTransport, whole.externalTransport), internal),
    cutTally: minus(plus(cutPairs.retainedCoincidentFlux, cutPairs.retainedNoncoincidentFlux), cuts),
    displacedWallLedger: minus(wall, recoveredWall),
  };
  return { order: ['mass', 'xMomentum', 'yMomentum', 'totalEnthalpy'], groups, categories, cutPairs,
    whole: { ...whole, remote, cuts, displacedWall: wall, recoveredDisplacedWall: recoveredWall },
    algebraErrors: errors, maximumAlgebraError: maximum(Object.values(errors).flat()),
    counts: { passages: groups.length, volumes: groups.reduce((s, g) => s + g.volumeCount, 0), boundaryHalfFaces: halfFaces },
    identity: 'F_displaced_wall = sum(local physical defects) - F_remote - F_upstream_cut - F_wake_banks - internal pressure mismatch.',
    physicalAcceptance: false };
}

// Signed pressure force on a CCW solid chain; exact linear nodal pressure
// quadrature. The caller identifies all base segments explicitly.
export function directLinearPressureChain(points, pressures, origin = { x: 0, y: 0 }) {
  require(points.length === pressures.length && points.length >= 2 && pressures.every(Number.isFinite), 'Incomplete pressure chain.');
  let x = 0, y = 0, moment = 0, length = 0;
  for (let i = 1; i < points.length; i++) {
    const a = points[i - 1], b = points[i], dx = b.x - a.x, dy = b.y - a.y;
    const mean = (pressures[i - 1] + pressures[i]) / 2, difference = pressures[i] - pressures[i - 1];
    x -= mean * dy; y += mean * dx;
    moment -= mean * (dx * ((a.x + b.x) / 2 - origin.x) + dy * ((a.y + b.y) / 2 - origin.y))
      + difference * (dx * dx + dy * dy) / 12;
    length += Math.hypot(dx, dy);
  }
  return { x, y, noseUpMoment: moment, segments: points.length - 1, length };
}
