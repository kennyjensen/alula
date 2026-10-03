// SPDX-License-Identifier: GPL-2.0-or-later
// Research observables only. Integrate pressure on the solid contour, never
// the displaced inviscid bank. Wake drag already includes the viscous loss;
// do not add the pressure integral or skin friction to it again.
import {pressureForces} from '../inviscid/pressure-forces.js';
import {wakeMomentumDrag} from '../viscous/wake-drag.js';
import {isentropicState} from '../potential/isentropic.js';
import {createContourTopology} from '../geometry/contour-topology.js';

export function quadCoupledCoefficients({surfaces,wakes,stagnation,bodies,alpha,mach,gamma=1.4,referenceChord,momentReference}){
 if(!(referenceChord>0)||![referenceChord,alpha,mach,momentReference?.x,momentReference?.y].every(Number.isFinite))
  throw new Error('Invalid coupled coefficient reference conditions.');
 const origin={x:momentReference.x/referenceChord,y:momentReference.y/referenceChord};
 const coefficients={cl:0,cm:0,cd:0,cx:0,cy:0,pressureIntegralDrag:0,elements:[],wakes:[],warnings:[]};
 for(let element=0;element<stagnation.length;element++){
  const upper=surfaces.find(s=>s.element===element&&s.side==='upper');
  const lower=surfaces.find(s=>s.element===element&&s.side==='lower');
  if(!upper||!lower)throw new Error('Missing coupled surface for pressure integration.');
  const contour=[...upper.stations.slice().reverse(),{...stagnation[element],cp:isentropicState(0,0,{mach,gamma}).cp},...lower.stations];
  const body=bodies?.find(b=>b.element===element);
  const base=body?.trailingEdge?.kind==='finite-base'?createContourTopology(body.points,body).base:null;
  let basePressureModel;
  if(base){
   const upperPoint=contour[0],lowerPoint=contour.at(-1),first=base.points[0],last=base.points.at(-1);
   if(upperPoint.x!==last.x||upperPoint.y!==last.y||lowerPoint.x!==first.x||lowerPoint.y!==first.y)
    throw new Error('Finite-base pressure loads must retain both original solid TE corners.');
   let arc=0;
   for(let k=1;k<base.points.length;k++){
    arc+=base.panels[k-1].length;
    const fraction=k===base.points.length-1?1:arc/base.length;
    contour.push({...base.points[k],cp:(1-fraction)*lowerPoint.cp+fraction*upperPoint.cp});
   }
   basePressureModel={kind:'TE-pressure interpolation on retained base',panelCount:base.panels.length,
    upperCp:upperPoint.cp,lowerCp:lowerPoint.cp};
  }
  const points=contour.map(p=>({x:p.x/referenceChord,y:p.y/referenceChord}));
  if(!points.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)))throw new Error('Nonfinite solid contour for pressure integration.');
  const force=pressureForces(points,contour.map(p=>p.cp),{alpha,momentOrigin:origin});
  for(const key of ['cl','cm','cx','cy'])coefficients[key]+=force[key];
  coefficients.pressureIntegralDrag+=force.cd;
  coefficients.elements.push({element,...force,pressureIntegralDrag:force.cd,...(basePressureModel?{basePressureModel}:{})});
 }
 for(const wake of wakes){
  try{
   const end=wake.stations.at(-1);
   if(end.wakeGap!==undefined&&(!Number.isFinite(end.wakeGap)||end.wakeGap<0||end.wakeGap>0))
    throw new Error('Wake exit has a nonzero or invalid prescribed dead-air gap; extend the wake before drag extrapolation.');
   const drag=wakeMomentumDrag({theta:end.theta/referenceChord,deltaStar:end.deltaStar/referenceChord,ue:end.ue},{mach,gamma});
   coefficients.cd+=drag.cd;coefficients.wakes.push({element:wake.element,station:end.index,...drag});
  }catch(error){coefficients.cd=null;coefficients.warnings.push(`Wake ${wake.element+1} drag unavailable: ${error.message}`);}
 }
 // An invalid wake must not silently turn total drag into a partial sum.
 if(coefficients.wakes.length!==stagnation.length)coefficients.cd=null;
 if(bodies?.some(b=>b.trailingEdge?.kind==='finite-base'))coefficients.warnings.push('Finite-base pressure interpolates the two TE pressures along every retained base segment; the separated base cavity is not resolved.');
 return {...coefficients,referenceChord,momentReference:{...momentReference},
  method:{lift:'BL-edge pressure integrated on the solid contour; shear contribution to lift omitted',
   moment:'BL-edge pressure integrated on the solid contour; positive nose-up; shear moment omitted',
   drag:'Sum of compressible Squire–Young extrapolations from each final wake station; includes viscous losses',
   pressure:'Isentropic pressure from solved physical edge speed; stagnation uses total pressure'},
  physicalAcceptance:false};
}
