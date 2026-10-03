// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour } from '../geometry/airfoil.js';
import { createContourTopology } from '../geometry/contour-topology.js';
import { solveMultielementViscous } from './multielement.js';
import { viscousObservables } from './observables.js';
import { provisionalViscousCoefficients } from './provisional-coefficients.js';
import { prepareAirfoilElement } from '../geometry/airfoil-element.js';

export const multielementLimitations='Incompressible, simultaneous integral boundary layers with shared displacement interaction and a curved wake from every element. Native XFOIL kernels and single-element results are verified. Multielement refinement is checked; experimental accuracy, wake confluence, strong separation and compressible Euler coupling remain unvalidated.';

// Serializable, dimensional browser/API result. Internal Newton contexts never
// cross the worker boundary. Unfinished loads carry explicit provisional status.
export function solveViscousAssembly({elements,alpha=0,mach=0,reynolds=1e6,ncrit=9,trips=[1,1],
  referenceChord=1,momentReference={x:referenceChord/4,y:0},maxIterations=60,
  wakeLength=2,wakeLengths,wakeCount=24,wakePaths,...options}={},controls={}){
  if(!Array.isArray(elements)||elements.length<1||elements.length>6)throw new Error('Supply between one and six lifting elements.');
  if(mach!==0)throw new Error('Coupled multielement flow currently requires Mach 0. Compressible Euler coupling is not yet validated.');
  if(!Number.isFinite(alpha)||Math.abs(alpha)>20)throw new Error('Enter incidence between −20° and 20°.');
  if(!Number.isFinite(reynolds)||reynolds<5e4||reynolds>5e7)throw new Error('Enter Reynolds number 50,000–50,000,000, based on reference chord.');
  if(!Number.isFinite(ncrit)||ncrit<0||ncrit>14)throw new Error('Ncrit must be between 0 and 14.');
  if(!Array.isArray(trips)||trips.length!==2||trips.some(x=>!Number.isFinite(x)||x<=0||x>1))throw new Error('Upper and lower trips must be in (0, 1] of each element chord.');
  if(!Number.isInteger(maxIterations)||maxIterations<1||maxIterations>200)throw new Error('Use 1–200 Newton iterations.');
  if(!Number.isFinite(referenceChord)||referenceChord<=0||!Number.isFinite(momentReference?.x)||!Number.isFinite(momentReference?.y))throw new Error('Invalid reference chord or moment reference.');
  if(wakeLengths&&(!Array.isArray(wakeLengths)||wakeLengths.length!==elements.length||wakeLengths.some(v=>!Number.isFinite(v)||v<=0)))throw new Error('Supply one positive absolute wake length per element.');
  elements=elements.map(prepareAirfoilElement);
  const topologies=elements.map(e=>e.trailingEdge?.kind==='finite-base'?createContourTopology(e.points,e):null);
  const contours=elements.map((e,k)=>topologies[k]?.points??prepareContour(e.points));
  const panelCount=contours.reduce((n,p)=>n+p.length-1,0);
  const wettedPanelCount=contours.reduce((n,p,k)=>n+(topologies[k]?.surface.panels.length??p.length-1),0);
  if(contours.some((p,k)=>(topologies[k]?.surface.points.length??p.length)<41)||wettedPanelCount>700)throw new Error('Use at least 40 wetted panels per element and at most 700 total. Refine transition regions carefully.');
  const normalized=contours.map((p,e)=>({name:elements[e].name,trips:elements[e].trips,
    ...(topologies[e]?{trailingEdge:structuredClone(elements[e].trailingEdge)}:{}),
    points:p.map(v=>({x:v.x/referenceChord,y:v.y/referenceChord}))}));
  const normalizedWakes=wakePaths?.map(path=>path.map(p=>({x:p.x/referenceChord,y:p.y/referenceChord})));
  const solverSettings={maxIterations,wakeLength,wakeLengths:wakeLengths??null,wakeCount,elementTrips:elements.map(e=>e.trips??trips),equationTolerance:controls.tolerance??1e-8,
    wakeTolerance:controls.wakeTolerance??1e-6,maxWakeIterations:controls.maxWakeIterations??40,
    wakeRelaxation:controls.wakeRelaxation??.8,initialization:options.initialization??(elements.length>1?'auto':'native'),
    jacobianMode:options.jacobianMode??'analytic',initialWakePaths:wakePaths??null};
  const raw=solveMultielementViscous({...options,elements:normalized,alpha,mach,reynolds,ncrit,trips,wakeLength,wakeLengths:wakeLengths?.map(v=>v/referenceChord),wakeCount,wakePaths:normalizedWakes},
    {wakeRelaxation:.8,maxIterations,...controls});
  const {system}=raw;
  if(system.initialization)solverSettings.actualInitialization=system.initialization.method;
  const residual=system.residual(raw.x),families=[0,0,0,0];
  residual.forEach((v,i)=>{families[i%4]=Math.max(families[i%4],Math.abs(v));});
  if(!raw.converged)return{model:'multielement-coupled-incompressible',status:'unconverged',
    alpha,mach,reynolds,ncrit,trips,referenceChord,momentReference,panelCount,solverSettings,
    ...provisionalViscousCoefficients(raw,{elements,alpha,mach,referenceChord,momentReference}),
    elements:contours.map((points,e)=>({name:elements[e].name??`Element ${e+1}`,points,cp:[]})),
    history:raw.history,wakeHistory:raw.wakeHistory,
    diagnostics:{equationResidual:Math.max(...families),wakeResidual:raw.wakeResidual,
      iterations:raw.history.filter(h=>h.iteration>0).length,reason:raw.reason},
    warnings:[`Coupled solve did not converge: ${raw.reason}. Any displayed coefficients are provisional values from the last iterate.`],limitations:multielementLimitations};
  const {surfaces,wakes,outputElements,cl,cm,cd,cdf,pressureIntegralDrag,maxRelativeSurfaceDisplacement,maxRelativeWakeDisplacement}=viscousObservables(raw,{elements,alpha,mach,referenceChord,momentReference});
  const warnings=[];
  const equationResidual=Math.max(...families),converged=raw.converged&&raw.wakeConverged&&equationResidual<=(controls.tolerance??1e-8);
  // A deliberately broad rejection threshold for the thin-layer asymptotic
  // model, not an assertion of experimental accuracy below that threshold.
  const maxRelativeDisplacement=Math.max(maxRelativeSurfaceDisplacement,maxRelativeWakeDisplacement),outsideModel=maxRelativeDisplacement>.1;
  if(!converged)warnings.push(`Coupled solve did not converge: ${raw.reason}. Coefficients withheld.`);
  if(outsideModel)warnings.push('Boundary-layer or wake displacement thickness exceeds 10% of its element chord. The thin-layer model is outside its supported range; coefficients withheld.');
  if(contours.some(p=>p.length<161))warnings.push('Refine to at least 160 panels per element to check transition and separation bubbles.');
  if(surfaces.some(s=>s.stations.some(p=>p.cf<0)))warnings.push('Local separation is predicted. Check boundary-layer thickness and panel refinement.');
  if(elements.length>1)warnings.push('Multielement coupling passes numerical refinement checks; comparison with independent multielement viscous reference data is still required.');
  if(reynolds<5e5||reynolds>3e6||Math.abs(alpha)>4)warnings.push('Outside the initial single-element Fortran validation range (Re 0.5–3 million, incidence ±4°).');
  const solved=converged&&!outsideModel;
  if(!solved){for(const e of outputElements){e.cl=null;e.cm=null;}for(const w of wakes)w.drag=null;}
  return{model:'multielement-coupled-incompressible',status:!converged?'unconverged':outsideModel?'outside-model':'solved',
    alpha,mach,reynolds,ncrit,trips,referenceChord,momentReference,panelCount,solverSettings,
    cl:solved?cl:null,cm:solved?cm:null,cd:solved?cd:null,cdf:solved?cdf:null,
    elements:outputElements,boundaryLayer:{surfaces,wakes},history:raw.history,wakeHistory:raw.wakeHistory,
    diagnostics:{equationResidual,equationTolerance:controls.tolerance??1e-8,couplingResidual:families[3],
      blResidual:Math.max(...families.slice(0,3)),wakeResidual:raw.wakeResidual,wakeTolerance:controls.wakeTolerance??1e-6,
      residualFamilies:families,iterations:raw.history.filter(h=>h.iteration>0).length,
      wakeIterations:raw.wakeHistory.length,transition:surfaces.map(s=>s.transition),maxRelativeDisplacement,maxRelativeSurfaceDisplacement,maxRelativeWakeDisplacement,
      wakeStations:wakes.reduce((n,w)=>n+w.stations.length,0),pressureIntegralDrag,reason:raw.reason},
    warnings,limitations:multielementLimitations};
}
