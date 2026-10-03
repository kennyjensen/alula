// SPDX-License-Identifier: GPL-2.0-or-later
// Conservative intrinsic-grid cell, following Drela's thesis (1986),
// sections 2.1–2.4. This local kernel is not yet a body/grid Euler solver.
// Mass and total enthalpy are eliminated exactly along each streamtube;
// the two interface pressures are recovered from normal momentum and the
// auxiliary pressure relation. Optional shared-section transportSpeeds
// supply the speed bias in BOTH momentum components (thesis 2.22–2.27).
// The default remains undissipated. See MSES manual 1.2.4 for the entropy row.
import { requirePositiveSimpleQuad } from '../geometry/simple-quadrilateral.js';
import { requirePositiveSimplePolygon } from '../geometry/simple-polygon.js';
const add=(a,b)=>({x:a.x+b.x,y:a.y+b.y});
const sub=(a,b)=>({x:a.x-b.x,y:a.y-b.y});
const scale=(a,k)=>({x:k*a.x,y:k*a.y});
const mean=(a,b)=>scale(add(a,b),.5);
const cross=(a,b)=>a.x*b.y-a.y*b.x;
const dot=(a,b)=>a.x*b.x+a.y*b.y;
const unit=a=>{const length=Math.hypot(a.x,a.y);if(!(length>0))throw new Error('Zero streamwise interval.');return scale(a,1/length);};

export function streamtubeCellGeometry(lower,upper,{geometryDomain='convex'}={}){
  if(!['convex','positive-simple'].includes(geometryDomain))throw new Error('Unknown streamtube geometry domain.');
  if(![lower,upper].every(row=>Array.isArray(row)&&row.length===3&&row.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y))))
    throw new Error('Supply three finite geometry nodes on each streamline.');
  // Convexity remains the default. The explicit research domain permits
  // concave intermediate finite-volume polygons, never crossed edges or
  // reversed areas. No absolute area substitutes for signed geometry.
  for(let i=0;i<2;i++){
    const points=[lower[i],lower[i+1],upper[i+1],upper[i]];
    if(geometryDomain==='positive-simple')requirePositiveSimpleQuad(points);
    else for(let j=0;j<4;j++)if(!(cross(sub(points[(j+1)%4],points[j]),sub(points[(j+2)%4],points[(j+1)%4]))>0))
      throw new Error('Folded or degenerate streamtube geometry.');
  }
  if (geometryDomain === 'positive-simple') {
    const ml = [mean(lower[0], lower[1]), mean(lower[1], lower[2])];
    const mu = [mean(upper[0], upper[1]), mean(upper[1], upper[2])];
    // The displayed quadrilaterals are not the staggered conservation
    // volumes. A simple concave quad can hide a crossed midpoint half-cell.
    for (let i = 0; i < 2; i++) {
      requirePositiveSimplePolygon([lower[i], ml[i], mu[i], upper[i]], 'streamtube half-volume');
      requirePositiveSimplePolygon([ml[i], lower[i + 1], upper[i + 1], mu[i]], 'streamtube half-volume');
    }
    requirePositiveSimplePolygon([ml[0], lower[1], ml[1], mu[1], upper[1], mu[0]], 'streamtube conservation-volume');
  }
  // Form relative vectors before averaging. Subtracting two rounded absolute
  // midpoints can lose a large fraction of a thin tube's normal width.
  const gaps=lower.map((p,i)=>sub(upper[i],p));
  const segments=[0,1].map(i=>mean(sub(lower[i+1],lower[i]),sub(upper[i+1],upper[i])));
  const directions=segments.map(unit);
  const sections=[mean(gaps[0],gaps[1]),mean(gaps[1],gaps[2])];
  const sides={lower:scale(sub(lower[2],lower[0]),.5),upper:scale(sub(upper[2],upper[0]),.5)};
  const streamwise=mean(sides.lower,sides.upper),transverse=mean(...sections),area=cross(streamwise,transverse);
  const normalAreas=sections.map((a,i)=>cross(directions[i],a));
  if(!(area>0)||!Number.isFinite(area)||normalAreas.some(a=>!(a>0)||!Number.isFinite(a)))throw new Error('Invalid signed streamtube area.');
  // B(lower)-B(upper) = cross(mean(s1),s2L-s2U)
  //                     + cross(s1L-s1U,mean(s2)).
  // Gap differences retain this numerator when the bank bends nearly match.
  // This is the same signed full-node correction, with no coefficient change.
  const bendDifference=cross(segments[0],sub(gaps[1],gaps[2]))+cross(sub(gaps[0],gaps[1]),segments[1]);
  const streamwiseLengths=segments.map(v=>Math.hypot(v.x,v.y));
  return{sections,directions,sides,streamwise,transverse,area,normalAreas,streamwiseLengths,
    pressureCurvature:bendDifference/(2*area)};
}

// Directional derivative of the signed intrinsic geometry. Both Euler
// linearization variants use this local cell kernel, so it lives beside the
// geometry it differentiates rather than in a separate forwarding module.
const dcross = (a, b, da, db) => cross(da, b) + cross(a, db);
const zeroGeometryRow = Object.freeze(Array.from({ length: 3 }, () => Object.freeze({ x: 0, y: 0 })));
function finiteGeometryRow(row) {
  if (!Array.isArray(row) || row.length !== 3) return false;
  for (let i = 0; i < 3; i++) if (!Number.isFinite(row[i].x) || !Number.isFinite(row[i].y)) return false;
  return true;
}

export function linearizeStreamtubeGeometry(lower, upper, value = streamtubeCellGeometry(lower, upper)) {
  lower = lower.map(p => ({ ...p })); upper = upper.map(p => ({ ...p }));
  const gaps = lower.map((p, i) => sub(upper[i], p));
  const segments = [0, 1].map(i => mean(sub(lower[i + 1], lower[i]), sub(upper[i + 1], upper[i])));
  const lengths = segments.map(v => Math.hypot(v.x, v.y));
  // Private synchronous scratch. None of these points escape in a returned
  // tangent, so successive column derivatives still have independent values.
  const dgaps = Array.from({ length: 3 }, () => ({ x: 0, y: 0 }));
  const dsegments = Array.from({ length: 2 }, () => ({ x: 0, y: 0 }));
  const gap12 = sub(gaps[1], gaps[2]), gap01 = sub(gaps[0], gaps[1]);
  const apply = ({ lower: dl = zeroGeometryRow, upper: du = zeroGeometryRow } = {}) => {
    if (!finiteGeometryRow(dl) || !finiteGeometryRow(du))
      throw new Error('Invalid intrinsic geometry tangent.');
    for (let i = 0; i < 3; i++) {
      dgaps[i].x = du[i].x - dl[i].x; dgaps[i].y = du[i].y - dl[i].y;
    }
    const streamwiseLengths = new Array(2), directions = new Array(2);
    for (let i = 0; i < 2; i++) {
      dsegments[i].x = ((dl[i + 1].x - dl[i].x) + (du[i + 1].x - du[i].x)) * .5;
      dsegments[i].y = ((dl[i + 1].y - dl[i].y) + (du[i + 1].y - du[i].y)) * .5;
      const dv = dsegments[i], t = value.directions[i];
      const projection = dot(t, dv), inverseLength = 1 / lengths[i];
      streamwiseLengths[i] = projection;
      directions[i] = { x: (dv.x - t.x * projection) * inverseLength, y: (dv.y - t.y * projection) * inverseLength };
    }
    const sections = [mean(dgaps[0], dgaps[1]), mean(dgaps[1], dgaps[2])];
    const sides = { lower: { x: (dl[2].x - dl[0].x) * .5, y: (dl[2].y - dl[0].y) * .5 },
      upper: { x: (du[2].x - du[0].x) * .5, y: (du[2].y - du[0].y) * .5 } };
    const streamwise = mean(sides.lower, sides.upper), transverse = mean(...sections);
    const area = dcross(value.streamwise, value.transverse, streamwise, transverse);
    const normalAreas = value.sections.map((v, i) => dcross(value.directions[i], v, directions[i], sections[i]));
    const dbendDifference = dcross(segments[0], gap12, dsegments[0], sub(dgaps[1], dgaps[2]))
      + dcross(gap01, segments[1], sub(dgaps[0], dgaps[1]), dsegments[1]);
    const pressureCurvature = (.5 * dbendDifference - value.pressureCurvature * area) / value.area;
    return { sections, directions, sides, streamwise, transverse, area, normalAreas, pressureCurvature, streamwiseLengths };
  };
  return { value, apply };
}

export function streamtubeSection({density,massFlow,normalArea,stagnationEnthalpy,gamma=1.4}){
  if(![density,massFlow,normalArea,stagnationEnthalpy,gamma].every(Number.isFinite)
    ||Math.min(density,massFlow,normalArea,stagnationEnthalpy)<=0||gamma<=1)
    throw new Error('Invalid streamtube thermodynamic inputs.');
  const q=massFlow/(density*normalArea),enthalpy=stagnationEnthalpy-.5*q*q,k=(gamma-1)/gamma;
  if(!(enthalpy>0))throw Object.assign(new Error('Nonpositive streamtube static enthalpy.'), {
    code: 'streamtube-static-enthalpy', diagnostics: { rho: density, massFlow, normalArea,
      stagnationEnthalpy, gamma, q, enthalpy, maximumPhysicalSpeed: Math.sqrt(2 * stagnationEnthalpy) },
  });
  const p=k*density*enthalpy,machSquared=q*q/((gamma-1)*enthalpy);
  if(![q,p,machSquared].every(Number.isFinite)||p<=0)throw new Error('Nonfinite streamtube thermodynamic state.');
  return{rho:density,q,p,machSquared,enthalpy,
    derivatives:{q:{density:-q/density,massFlow:q/massFlow,normalArea:-q/normalArea,stagnationEnthalpy:0},
      p:{density:k*(stagnationEnthalpy+.5*q*q),massFlow:-k*q/normalArea,
        normalArea:k*density*q*q/normalArea,stagnationEnthalpy:k*density}}};
}

export function evaluateStreamtubeCell({lower,upper,densities,massFlow,stagnationEnthalpy,gamma=1.4,pressureCorrectionFactor=.1,geometryDomain='convex',transportSpeeds}){
  if(!Array.isArray(densities)||densities.length!==2)throw new Error('Supply a density at each cross-section.');
  if(!Number.isFinite(pressureCorrectionFactor)||pressureCorrectionFactor<0)throw new Error('Invalid auxiliary pressure correction factor.');
  const geometry=streamtubeCellGeometry(lower,upper,{geometryDomain}),{directions,sections,sides,streamwise,transverse,area}=geometry;
  const states=densities.map((density,i)=>streamtubeSection({density,massFlow,normalArea:geometry.normalAreas[i],stagnationEnthalpy,gamma}));
  const [a,b]=states,pMean=.5*(a.p+b.p),machSquared=.5*(a.machSquared+b.machSquared);
  const biased=transportSpeeds!==undefined;
  if(biased&&(!Array.isArray(transportSpeeds)||transportSpeeds.length!==2
    ||!transportSpeeds.every(q=>Number.isFinite(q)&&q>=0)))throw new Error('Invalid streamtube transport speeds.');
  const [qa,qb]=biased?transportSpeeds:[a.q,b.q];
  // Eq. 2.16 uses signed cross products, not absolute values. This
  // correction controls the geometry sawtooth; it is not shock dissipation.
  const pressureCorrection=machSquared<1?pressureCorrectionFactor*gamma*pMean*machSquared*(1-machSquared)*geometry.pressureCurvature:0;
  const normalInertia=massFlow*(qa*dot(directions[0],transverse)-qb*dot(directions[1],transverse))/area;
  const pressureDifference=normalInertia+pressureCorrection*cross(sections[0],sections[1])/area;
  const pressureSum=a.p+b.p+2*pressureCorrection;
  const interfacePressure={lower:.5*(pressureSum-pressureDifference),upper:.5*(pressureSum+pressureDifference)};
  if(!Object.values(interfacePressure).every(p=>Number.isFinite(p)&&p>0))throw Object.assign(new Error('Nonpositive or nonfinite streamline interface pressure.'), {
    code: 'streamtube-interface-pressure', diagnostics: { interfacePressure, pressureSum, pressureDifference, normalInertia },
  });
  // Outward conservative flux sign (the negative of thesis Eq. 2.18).
  // Matching adjacent cells' interface pressures supplies reduced N-momentum.
  const streamwiseResidual=massFlow*(qb*dot(directions[1],streamwise)-qa*dot(directions[0],streamwise))/area
    +b.p-a.p-pressureCorrection*cross(sides.lower,sides.upper)/area;
  // S = log((h/h0)^(gamma/(gamma-1)) / (p/p0)). With a perfect gas,
  // Delta S = log(h_b/h_a)/(gamma-1) - log(rho_b/rho_a). Stable relative
  // differences retain small smooth-flow changes. This row cannot replace
  // conservative momentum across a shock. Preserve the physical entropy
  // jump for diagnostics even when the optional biased row is requested.
  const entropyJump=Math.log1p((b.enthalpy-a.enthalpy)/a.enthalpy)/(gamma-1)-Math.log1p((b.rho-a.rho)/a.rho);
  if(biased){
    const artificialEnthalpies=transportSpeeds.map(q=>stagnationEnthalpy-.5*q*q);
    if(!artificialEnthalpies.every(h=>Number.isFinite(h)&&h>0))throw new Error('Nonpositive artificial streamtube enthalpy.');
    // S_tilde uses h_tilde but PHYSICAL p. The log(h)/(gamma-1)-log(rho)
    // simplification is valid only for the physical state. Relative changes
    // give exact default-row parity when supplied speeds equal actual q.
    const logRatios=states.map((state,i)=>Math.log1p(.5*(state.q-transportSpeeds[i])*(state.q+transportSpeeds[i])/state.enthalpy));
    const artificialEntropyJump=entropyJump+gamma/(gamma-1)*(logRatios[1]-logRatios[0]);
    const isentropicResidual=-pMean*artificialEntropyJump;
    if(![artificialEntropyJump,isentropicResidual].every(Number.isFinite))throw new Error('Nonfinite artificial streamtube entropy.');
    return{geometry,states,pressureCorrection,interfacePressure,streamwiseResidual,entropyJump,isentropicResidual,
      transportSpeeds:[...transportSpeeds],artificialEnthalpies,artificialEntropyJump};
  }
  const isentropicResidual=-pMean*entropyJump;
  return{geometry,states,pressureCorrection,interfacePressure,streamwiseResidual,entropyJump,isentropicResidual};
}
