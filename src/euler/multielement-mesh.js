// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour, validateAssembly } from '../geometry/airfoil.js';
import { buildMesh } from './mesh.js';

// Verification topology: bodies must be x-monotone and vertically orderable.
// Exact polygon vertices are retained. Zero-thickness inlet/wake cuts are
// shared INTERNAL faces, never walls. Cross-lines are fixed in this backend;
// their arbitrary placement does not impose zero cross-flow.
export function multielementMesh(contours, { rows = 2, padding = 1, subdivisions = 1,
  outerBoundary = 'farfield',outerGrowth,wallStretch=3,wakePaths,normalGrowth,wakeNormalRelaxation } = {}) {
  if (!Array.isArray(contours) || !contours.length || contours.length > 3) throw new Error('Euler mesh needs one to three elements.');
  if (!Number.isInteger(rows) || rows < 1 || rows > 80 || !Number.isInteger(subdivisions) || subdivisions < 1 || subdivisions > 8 || !(padding > 0) || !Number.isFinite(padding)) throw new Error('Invalid multielement mesh controls.');
  if (!['farfield', 'channel'].includes(outerBoundary)) throw new Error('Unknown outer boundary model.');
  if((outerGrowth!==undefined&&(!(outerGrowth>1)||outerGrowth>3))||!(wallStretch>=1)||wallStretch>8)throw new Error('Invalid Euler mesh grading.');
  if(normalGrowth!==undefined&&(!(normalGrowth>1)||normalGrowth>3))throw new Error('Invalid normal mesh growth.');
  if(wakeNormalRelaxation!==undefined&&(!(wakeNormalRelaxation>0)||!Number.isFinite(wakeNormalRelaxation)))throw new Error('Invalid downstream normal-spacing relaxation.');
  const prepared = contours.map(c => prepareContour(c));
  validateAssembly(prepared);
  if(wakePaths&&(!Array.isArray(wakePaths)||wakePaths.length!==prepared.length))throw new Error('Supply one mesh wake path per body.');
  const bodies = prepared.map((points, element) => {
    const xmin = Math.min(...points.map(p => p.x)); const xmax = Math.max(...points.map(p => p.x));
    const lo = points.find(p => p.x === xmin); const hi = points.find(p => p.x === xmax);
    const width = xmax - xmin;
    if (!(width > 0)) throw new Error('Unresolved streamwise body extent.');
    const wake=wakePaths?.[element];
    if(wake){
      if(!Array.isArray(wake)||wake.length<2||wake.some(p=>!Number.isFinite(p.x)||!Number.isFinite(p.y))
        ||Math.hypot(wake[0].x-points[0].x,wake[0].y-points[0].y)>width*1e-10
        ||Math.abs(wake[0].x-xmax)>width*1e-10||wake.some((p,i)=>i>0&&p.x-wake[i-1].x<=width*1e-12)){
        throw new Error('The cross-line mesh requires a downstream x-monotone wake starting at the sharp trailing edge.');
      }
    }
    const wakeY=x=>{
      if(!wake)return hi.y;
      let j=1;while(j<wake.length-1&&wake[j].x<x)j++;
      const a=wake[j-1],b=wake[j];return a.y+(b.y-a.y)*(x-a.x)/(b.x-a.x);
    };
    for (let k = 0; k < points.length - 1; k++) if (Math.abs(points[k + 1].x - points[k].x) < width * 1e-13) {
      throw new Error('Vertical body edges are outside the reference mesh topology.');
    }
    const sample = x => {
      if (x <= xmin) return [lo.y, lo.y];
      if (x >= xmax) {const y=wakeY(x);return[y,y];}
      const values = [];
      for (let k = 0; k < points.length - 1; k++) {
        const a = points[k]; const b = points[k + 1];
        if (Math.abs(b.x - a.x) < width * 1e-13) {
          if (Math.abs(x - a.x) < width * 1e-13) throw new Error('Vertical body edges are outside the reference mesh topology.');
          continue;
        }
        if (x >= Math.min(a.x, b.x) && x <= Math.max(a.x, b.x)) values.push(a.y + (b.y - a.y) * (x - a.x) / (b.x - a.x));
      }
      values.sort((a, b) => a - b);
      const distinct = values.filter((v, k) => !k || Math.abs(v - values[k - 1]) > width * 1e-12);
      if (distinct.length !== 2) throw new Error('Body must have exactly two intersections with every interior cross-line.');
      return distinct;
    };
    return { element, points, xmin, xmax, sample, center: sample((xmin + xmax) / 2).reduce((a, b) => a + b) / 2 };
  }).sort((a, b) => a.center - b.center);
  const xValues = prepared.flatMap(c => c.map(p => p.x)).sort((a, b) => a - b);
  const chord = xValues.at(-1) - xValues[0];
  const unique = xValues.filter((x, i) => !i || Math.abs(x - xValues[i - 1]) > chord * 1e-12);
  // Validate geometry at interval midpoints too, before any mesh allocation.
  for (const body of bodies) for (let i = 1; i < unique.length; i++) body.sample((unique[i - 1] + unique[i]) / 2);
  const extension=(start,spacing,direction)=>{
    const points=[];let distance=0,step=spacing;
    while(distance<padding*chord){distance=Math.min(padding*chord,distance+step);points.push(start+direction*distance);step*=outerGrowth;}
    return points;
  };
  let xbase = outerGrowth===undefined?[unique[0] - padding * chord, ...unique, unique.at(-1) + padding * chord]
    :[...extension(unique[0],unique[1]-unique[0],-1).reverse(),...unique,...extension(unique.at(-1),unique.at(-1)-unique.at(-2),1)];
  if(wakePaths){
    const wakeX=wakePaths.flatMap(w=>w.map(p=>p.x));
    if(wakeX.some(x=>x>=xbase.at(-1)))throw new Error('The complete modeled wake must fit inside the mesh farfield.');
    const sorted=[...xbase,...wakeX].sort((a,b)=>a-b);
    xbase=sorted.filter((x,i)=>!i||x-sorted[i-1]>chord*1e-12);
  }
  const xs = [xbase[0]];
  for (let i = 1; i < xbase.length; i++) for (let k = 1; k <= subdivisions; k++) xs.push(k === subdivisions ? xbase[i] : xbase[i - 1] + (xbase[i] - xbase[i - 1]) * k / subdivisions);
  const allY = prepared.flatMap(c => c.map(p => p.y));
  const ymin = Math.min(...allY) - padding * chord; const ymax = Math.max(...allY) + padding * chord;
  const vertices = []; const connectivity = []; const boundaries = []; const blocks = [];
  const cuts = [];
  const add = (x, y) => { vertices.push({ x, y }); return vertices.length - 1; };
  // Each cross-line has a pair of IDs on every body. They collapse to one ID
  // upstream/downstream of the body; this enforces conforming face connectivity.
  const columns = xs.map(x => {
    const surfaces = bodies.map(b => b.sample(x));
    if(surfaces[0][0]<=ymin||surfaces.at(-1)[1]>=ymax)throw new Error('Body or wake cut crosses the outer boundary.');
    for (let e = 0; e < surfaces.length - 1; e++) if (!(surfaces[e + 1][0] - surfaces[e][1] > chord * 1e-10)) throw new Error('Body/dividing cuts cross: this geometry needs a general MSET grid.');
    const ids = surfaces.map(([lo, hi]) => { const a = add(x, lo); return [a, hi === lo ? a : add(x, hi)]; });
    const bottom = add(x, ymin); const top = add(x, ymax);
    const gaps = [];
    for (let g = 0; g <= bodies.length; g++) {
      const low = g === 0 ? bottom : ids[g - 1][1]; const high = g === bodies.length ? top : ids[g][0];
      const column = [low];
      for (let j = 1; j < rows; j++) {
        const eta = j / rows;
        // Resolve body walls while keeping the outer domain affordable for the
        // dense verification solve. Gap blocks cluster toward both bodies.
        const downstream=body=>wakeNormalRelaxation===undefined?1:Math.exp(-Math.max(0,x-body.xmax)/(wakeNormalRelaxation*(body.xmax-body.xmin)));
        const relaxation=g===0?downstream(bodies[0]):g===bodies.length?downstream(bodies.at(-1)):Math.max(downstream(bodies[g-1]),downstream(bodies[g]));
        const growth=normalGrowth===undefined?undefined:1+(normalGrowth-1)*relaxation;
        const stretch=t=>growth===undefined?(1-relaxation)*t+relaxation*t**wallStretch:Math.abs(growth-1)<1e-10?t:Math.expm1(rows*t*Math.log(growth))/Math.expm1(rows*Math.log(growth));
        const fraction = g === 0 ? 1 - stretch(1 - eta) : g === bodies.length ? stretch(eta)
          : (1-relaxation)*eta+relaxation*0.5 * (1 + Math.tanh(2 * (2 * eta - 1)) / Math.tanh(2));
        column.push(add(x, vertices[low].y + fraction * (vertices[high].y - vertices[low].y)));
      }
      column.push(high); gaps.push(column);
    }
    return { ids, gaps };
  });
  const wall = (a, b, element, side) => boundaries.push({ a, b, type: 'wall', element, side });
  for (let g = 0; g <= bodies.length; g++) {
    const cellIndices = [];
    for (let i = 0; i < xs.length - 1; i++) for (let j = 0; j < rows; j++) {
      const left = columns[i].gaps[g]; const right = columns[i + 1].gaps[g];
      cellIndices.push(connectivity.length); connectivity.push([left[j], right[j], right[j + 1], left[j + 1]]);
    }
    blocks.push({ gap: g, cells: cellIndices });
    for (let j = 0; j < rows; j++) {
      boundaries.push({ a: columns[0].gaps[g][j], b: columns[0].gaps[g][j + 1], type: outerBoundary === 'channel' ? 'inlet' : 'farfield' });
      boundaries.push({ a: columns.at(-1).gaps[g][j], b: columns.at(-1).gaps[g][j + 1], type: outerBoundary === 'channel' ? 'outlet' : 'farfield' });
    }
  }
  for (let i = 0; i < xs.length - 1; i++) {
    const left = columns[i]; const right = columns[i + 1]; const xm = (xs[i] + xs[i + 1]) / 2;
    for (const g of [0, bodies.length]) {
      const j = g === 0 ? 0 : rows;
      boundaries.push({ a: left.gaps[g][j], b: right.gaps[g][j], type: outerBoundary === 'channel' ? 'wall' : 'farfield', side: g === 0 ? 'lower' : 'upper' });
    }
    bodies.forEach((body, b) => {
      if (xm > body.xmin && xm < body.xmax) {
        wall(left.ids[b][0], right.ids[b][0], body.element, 'lower');
        wall(left.ids[b][1], right.ids[b][1], body.element, 'upper');
      } else cuts.push({ a: left.ids[b][0], b: right.ids[b][0], element: body.element, type: xm < body.xmin ? 'inlet-cut' : 'wake-cut' });
    });
  }
  const mesh = buildMesh(vertices, connectivity, boundaries);
  const byEdge = new Map(mesh.faces.map((f, i) => [f.a < f.b ? `${f.a}:${f.b}` : `${f.b}:${f.a}`, i]));
  for (const cut of cuts) { cut.face = byEdge.get(cut.a < cut.b ? `${cut.a}:${cut.b}` : `${cut.b}:${cut.a}`); if (mesh.faces[cut.face]?.neighbor === null) throw new Error('Disconnected fluid cut.'); }
  return { ...mesh, blocks, cuts, contours: prepared, topology: 'fixed-cross-line-multielement',
    controls: { rows, padding, subdivisions, outerBoundary,outerGrowth,wallStretch,wakePaths,normalGrowth,wakeNormalRelaxation }, bounds: { xmin: xs[0], xmax: xs.at(-1), ymin, ymax } };
}
