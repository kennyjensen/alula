// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour } from '../geometry/airfoil.js';
import { solveMultielementPotentialViscous } from './coupled.js';
import { viscousObservables } from '../viscous/observables.js';
import { provisionalViscousCoefficients } from '../viscous/provisional-coefficients.js';
import { displacementMass } from '../viscous/mass.js';
import { potentialMeshSnapshot, subcriticalMeshOptions } from './mesh-preview.js';

export const subcriticalLimitations='Simultaneous multielement integral boundary layers and displacement wakes coupled to compressible, uniform-entropy full-potential flow. All local flow must remain subsonic. Native BL kernels, single-element comparisons, conservation and refinement checks are available. Experimental multielement accuracy, wake confluence, strong separation and shock waves are not validated.';

// Serializable browser/API adapter. One common reference chord defines Re,
// coefficients and the requested distance beyond the assembly for all wakes.
export function solveSubcriticalAssembly({elements,alpha=0,mach=.2,reynolds=1e6,ncrit=9,trips=[1,1],
  referenceChord=1,momentReference={x:referenceChord/4,y:0},maxIterations=80,
  wakeLength=2,wakeLengths,wakeCount=48,wakePaths,initialization,jacobianMode='analytic'}={},controls={}){
  if(!Array.isArray(elements)||elements.length<1||elements.length>6)throw new Error('Supply between one and six lifting elements.');
  if(!Number.isFinite(mach)||mach<0||mach>=1)throw new Error('Subcritical flow requires freestream Mach in [0, 1) and subsonic local velocities.');
  if(!Number.isFinite(alpha)||Math.abs(alpha)>20)throw new Error('Enter incidence between −20° and 20°.');
  if(!Number.isFinite(reynolds)||reynolds<5e4||reynolds>5e7)throw new Error('Enter Reynolds number 50,000–50,000,000, based on reference chord.');
  if(!Number.isFinite(ncrit)||ncrit<0||ncrit>14)throw new Error('Ncrit must be between 0 and 14.');
  if(!Array.isArray(trips)||trips.length!==2||trips.some(x=>!Number.isFinite(x)||x<=0||x>1))throw new Error('Upper and lower trips must be in (0, 1] of each element chord.');
  maxIterations=controls.maxIterations??maxIterations;
  if(!Number.isInteger(maxIterations)||maxIterations<1||maxIterations>200)throw new Error('Use 1–200 Newton iterations.');
  if(!Number.isFinite(referenceChord)||referenceChord<=0||!Number.isFinite(momentReference?.x)||!Number.isFinite(momentReference?.y))throw new Error('Invalid reference chord or moment reference.');
  if(!Number.isFinite(wakeLength)||wakeLength<=0||!Number.isInteger(wakeCount)||wakeCount<4||wakeCount>96)throw new Error('Use a positive wake extent and 4–96 wake stations per element.');
  if(wakeLengths&&(!Array.isArray(wakeLengths)||wakeLengths.length!==elements.length||wakeLengths.some(v=>!Number.isFinite(v)||v<=0)))throw new Error('Supply one positive absolute wake length per element.');
  const contours=elements.map(e=>prepareContour(e.points)),panelCount=contours.reduce((n,p)=>n+p.length-1,0);
  if(contours.some(p=>p.length<41)||panelCount>700)throw new Error('Use at least 40 panels per element and at most 700 total.');
  const normalized=contours.map((p,e)=>({name:elements[e].name,trips:elements[e].trips,points:p.map(v=>({x:v.x/referenceChord,y:v.y/referenceChord}))}));
  const target=Math.max(...normalized.flatMap(e=>e.points.map(p=>p.x)))+wakeLength;
  const lengths=wakeLengths?wakeLengths.map(v=>v/referenceChord):normalized.map(e=>target-e.points[0].x);
  const normalizedWakes=wakePaths?.map(path=>path.map(p=>({x:p.x/referenceChord,y:p.y/referenceChord})));
  const mesh=subcriticalMeshOptions(controls.mesh);
  const settings={farfield:'multipole',edgeVelocitySampling:'station-average',bodyMassInterpolation:'hermite',...controls,maxIterations,mesh,
    fillLevel:controls.fillLevel??2,wakeRelaxation:controls.wakeRelaxation??.5,
    onMesh:controls.onMesh?(grid=>controls.onMesh(potentialMeshSnapshot(grid,referenceChord))):undefined};
  initialization??=elements.length>1?'auto':'native';
  const solverSettings={maxIterations,wakeLength,wakeLengths:lengths.map(v=>v*referenceChord),wakeCount,elementTrips:elements.map(e=>e.trips??trips),
    wakeExtentConvention:'arc lengths to the initial common downstream extent, in input coordinate units',
    equationTolerance:controls.tolerance??1e-8,wakeTolerance:controls.wakeTolerance??1e-6,
    maxWakeIterations:controls.maxWakeIterations??40,wakeRelaxation:settings.wakeRelaxation,
    mesh,fluxQuadrature:controls.fluxQuadrature??3,fillLevel:settings.fillLevel,linearTolerance:controls.linearTolerance??1e-7,
    eventStepFraction:controls.eventStepFraction??.25,machContinuation:controls.machContinuation??elements.length>=3,
    assemblyInitialization:controls.assemblyInitialization??elements.length>=3,
    edgeVelocitySampling:settings.edgeVelocitySampling,
    bodyMassInterpolation:settings.bodyMassInterpolation,
    edgeVelocityContinuation:controls.edgeVelocityContinuation??(elements.length>=3&&settings.edgeVelocitySampling!=='vertex'),
    farfield:settings.farfield,
    linearBackend:controls.linearBackend??'klu',coupledPreconditioner:controls.linearBackend==='gmres'?'full-Jacobian ILU':null,
    initialization,jacobianMode,initialWakePaths:wakePaths??null};
  const raw=solveMultielementPotentialViscous({elements:normalized,alpha,mach,reynolds,ncrit,trips,
    wakeLengths:lengths,wakeCount,wakePaths:normalizedWakes,initialization,jacobianMode},settings);
  const common={model:'multielement-coupled-subcritical',alpha,mach,reynolds,ncrit,trips,referenceChord,momentReference,
    panelCount,solverSettings,history:raw.history,wakeHistory:raw.wakeHistory,assemblyInitialization:raw.assemblyInitialization,limitations:subcriticalLimitations,
    // The final solve's actual deformed grid, in input coordinate units.
    // Omit operators/flow states from this compact display/export snapshot.
    mesh:potentialMeshSnapshot(raw.system.mesh,referenceChord)};
  const diagnostics={equationResidual:raw.diagnostics?.residual??Infinity,equationTolerance:solverSettings.equationTolerance,
    wakeResidual:raw.wakeResidual,wakeTolerance:solverSettings.wakeTolerance,
    iterations:raw.history.filter(h=>h.iteration>0).length,wakeIterations:raw.wakeHistory.length,
    cells:raw.system.mesh.cells.length,unknowns:raw.system.n,maxMach:raw.stage==='mach-initialization'?null:raw.diagnostics?.maxMach??null,reason:raw.reason};
  if(!raw.converged||!raw.wakeConverged||diagnostics.equationResidual>solverSettings.equationTolerance)return{...common,status:'unconverged',
    ...provisionalViscousCoefficients(()=>({system:raw.system.bl,...raw.system.bl.decode(raw.x.subarray(raw.system.blOffset)),x:raw.x.subarray(raw.system.blOffset)}),
      {elements,alpha,mach:raw.system.conditions.mach,referenceChord,momentReference,inviscidReference:false}),
    elements:contours.map((points,e)=>({name:elements[e].name??`Element ${e+1}`,points,cp:[]})),
    diagnostics,warnings:[`Coupled subcritical flow did not converge: ${raw.reason}. Any displayed coefficients are provisional values from the last iterate at Mach ${raw.system.conditions.mach}.`]};
  const {bl,outer,blOffset}=raw.system,x=raw.x.subarray(blOffset);
  const o=viscousObservables({system:bl,...raw.boundaryLayer,x},{elements,alpha,mach,referenceChord,momentReference,inviscidReference:false});
  const families=[0,0,0,0];raw.residual.subarray(blOffset).forEach((v,i)=>{families[i%4]=Math.max(families[i%4],Math.abs(v));});
  let farfieldMass=0;raw.system.mesh.faces.forEach((f,i)=>{if(f.boundary?.type==='farfield')farfieldMass+=raw.fluxes[i];});
  const downstreamMass=outer.wakes.reduce((sum,w)=>sum+displacementMass(x[4*w.end+3],x[4*w.end+2]*bl.thicknessScale,{mach}).value,0);
  Object.assign(diagnostics,{couplingResidual:families[3],blResidual:Math.max(...families.slice(0,3)),residualFamilies:families,
    massResidual:raw.diagnostics.massResidual,massBudgetError:Math.abs(farfieldMass-downstreamMass),farfieldMass,downstreamWakeMass:downstreamMass,
    transition:o.surfaces.map(s=>s.transition),wakeStations:o.wakes.reduce((n,w)=>n+w.stations.length,0),
    maxRelativeSurfaceDisplacement:o.maxRelativeSurfaceDisplacement,maxRelativeWakeDisplacement:o.maxRelativeWakeDisplacement,
    maxRelativeDisplacement:Math.max(o.maxRelativeSurfaceDisplacement,o.maxRelativeWakeDisplacement),pressureIntegralDrag:o.pressureIntegralDrag});
  const outside=diagnostics.maxRelativeDisplacement>.1,warnings=[];
  if(outside)warnings.push('Displacement thickness exceeds 10% of its element chord. The thin-layer model is outside its supported range; coefficients withheld.');
  if(contours.some(p=>p.length<161))warnings.push('Refine to at least 160 panels per element and compare Cp, transition and drag.');
  if(o.surfaces.some(s=>s.stations.some(p=>p.cf<0)))warnings.push('Local separation is predicted. Check boundary-layer thickness and mesh refinement.');
  warnings.push('Check surface spacing, outer-mesh spacing, farfield distance and wake extent. Numerical convergence alone does not establish experimental accuracy.');
  if(elements.length>1)warnings.push('Comparison with independent multielement viscous reference data is still required.');
  if(mach>.2||reynolds<5e5||reynolds>3e6||Math.abs(alpha)>4)warnings.push('Outside the initial single-element comparison range: Mach 0–0.2, Re 0.5–3 million, incidence ±4°.');
  if(outside){for(const e of o.outputElements){e.cl=null;e.cm=null;}for(const w of o.wakes)w.drag=null;}
  return{...common,status:outside?'outside-model':'solved',cl:outside?null:o.cl,cm:outside?null:o.cm,cd:outside?null:o.cd,cdf:outside?null:o.cdf,
    elements:o.outputElements,boundaryLayer:{surfaces:o.surfaces,wakes:o.wakes},diagnostics,warnings};
}
