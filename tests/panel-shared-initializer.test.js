import test from'node:test';import assert from'node:assert/strict';import fs from'node:fs';
import{initializePrescribedPanelProfiles}from'../src/viscous/shared-field-initializer.js';
import{createPanelContext}from'../src/viscous/context.js';
import{createCoupledAssembly}from'../src/viscous/assembly.js';
import{splind,seval}from'../src/viscous/xfoil/spline.js';
import{naca4,transform}from'../src/geometry/airfoil.js';
const close=(a,b,t=1e-10)=>assert.ok(Math.abs(a-b)<=t*Math.max(1e-6,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
test('shared prescribed MRCHUE adapter retains original Fortran profile states and transition including warning cases',t=>{
 const f=JSON.parse(fs.readFileSync('docs/current-quad-automatic-mrchue-before.json'));let count=0,max=0,warnings=0;
 for(const c of f.cases){const r=initializePrescribedPanelProfiles([c.profile,c.profile],c.parameters);warnings+=r.localConvergenceWarnings.length;
 for(const s of r.surfaces){assert.equal(s.transition,c.native.transition);assert.equal(s.forced,c.native.forced);close(s.s,c.native.s);
 s.states.forEach((p,i)=>{for(const k of ['ue','aux','theta','deltaStar']){close(p[k],c.native.states[i][k]);max=Math.max(max,Math.abs(p[k]-c.native.states[i][k])/Math.max(1e-6,Math.abs(c.native.states[i][k])));}count++;});}}
 assert.ok(warnings>0);t.diagnostic(JSON.stringify({stations:count,max,warnings}));
});
test('geometry-only context preserves native geometry while omitting dense panel systems',()=>{
 const p=naca4('2412',40),a=createPanelContext(p,3),b=createPanelContext(p,3,{geometryOnly:true});
 for(const k of ['X','Y','S','XP','YP','XLE','YLE','XTE','YTE','CHORD','ANTE','ASTE','DSTE','SHARP'])assert.deepEqual(b[k],a[k]);
 assert.equal(b.AIJ.length,0);assert.equal(b.BIJ.length,0);assert.ok(Number.isFinite(b.SLE));
});
const finite=()=>{const p=naca4('0012',40).map((v,i,a)=>({x:v.x,y:v.y+(i<(a.length-1)/2?.002:-.002)*v.x}));return{points:[...p,{...p[0]}],trailingEdge:{kind:'finite-base',upperIndex:0,lowerIndex:p.length-1}};};
test('explicit shared initialization maps actual signs and material trips without isolated solves or geometry identity assumptions',()=>{
 const e=finite(),input={elements:[e,{points:transform(naca4('0012',40),{chord:.3,x:1.1,y:-.2})}],alpha:4,trips:[.1,.2],wakeCount:8,initialization:'shared-mrchue'};
 const s=createCoupledAssembly(input),d=s.initialization;assert.equal(d.method,'shared-field MRCHUE');assert.equal(d.isolatedSolves,0);
 for(let id=0;id<s.total;id++)assert.equal(Math.sign(s.initial[4*id+3]),Math.sign(s.outer.q0[id]));
 assert.equal(s.surfaces.length,4);assert.equal(s.admissible(s.initial),true);
 for(const e of d.elements)for(const surf of e.surfaces){assert.ok(surf.tripS>0);assert.ok(surf.tripS<surf.requestedProfile.at(-1).s);}
 const seed=s.exportSeed(s.initial);assert.equal(seed.initialization.method,'shared-field MRCHUE');
 const resumed=createCoupledAssembly({...input,seed});assert.deepEqual(resumed.initial,s.initial);assert.equal(resumed.initialization.resumed,true);
});
test('finite multielement auto chooses shared startup, while sharp single-element native remains unchanged',()=>{
 const input={elements:[finite(),{points:transform(naca4('0012',40),{chord:.3,x:1.1,y:-.2})}],alpha:4,trips:[1,1],wakeCount:8};
 const a=createCoupledAssembly({...input,initialization:'auto'}),b=createCoupledAssembly({...input,initialization:'shared-mrchue'});
 assert.deepEqual(a.initial,b.initial);assert.deepEqual(a.initialization,b.initialization);
 const sharp=createCoupledAssembly({elements:[{points:naca4('0012',40)}],alpha:4,wakeCount:8,initialization:'march'});assert.equal(sharp.initialization,undefined);
});

test('a material trip before the first physical station is preserved by an upstream virtual similarity station',t=>{
 const input={elements:[finite(),{points:transform(naca4('0012',40),{chord:.3,x:1.1,y:-.2})}],alpha:4,trips:[1,1],wakeCount:8,initialization:'shared-mrchue'};
 const base=createCoupledAssembly(input),b=base.outer.bodies[0],meta=base.initialization.elements[0],first=meta.surfaces[1].requestedProfile[0];
 const panel=createPanelContext(b.points,0,{geometryOnly:true}),cx=panel.XTE-panel.XLE,cy=panel.YTE-panel.YLE,c2=cx*cx+cy*cy;
 const values=Float64Array.from(b.points,p=>((p.x-panel.XLE)*cx+(p.y-panel.YLE)*cy)/c2),derivative=new Float64Array(values.length);
 splind(values,derivative,b.s,b.s.length,-999,-999);
 const targetArc=meta.stagnation.s+.5*first.s,fraction=seval(targetArc,values,derivative,b.s,b.s.length);
 assert.ok(fraction>0&&fraction<1);
 const elements=structuredClone(input.elements);elements[0].trips=[.2,fraction];
 const result=createCoupledAssembly({...input,elements}),surface=result.initialization.elements[0].surfaces[1];
 assert.equal(surface.virtualStation,true);assert.equal(surface.transition,0);assert.equal(surface.requestedProfile[0].id,null);
 assert.equal(surface.requestedProfile[1].id,first.id);assert.equal(surface.requestedProfile[1].s,first.s);
 assert.equal(surface.requestedProfile[0].s,.5*surface.tripS);assert.ok(surface.tripS<first.s);
 const reconstructed=seval(surface.tripArc,values,derivative,b.s,b.s.length),error=Math.abs(reconstructed-fraction);
 assert.ok(error<2e-10,JSON.stringify({error,reconstructed,fraction}));
 for(let id=0;id<result.total;id++)assert.equal(Math.sign(result.initial[4*id+3]),Math.sign(result.outer.q0[id]));
 t.diagnostic(JSON.stringify({fraction,tripArc:surface.tripArc,tripS:surface.tripS,firstPhysicalS:first.s,virtualS:surface.requestedProfile[0].s,reconstructedFractionError:error}));
});
