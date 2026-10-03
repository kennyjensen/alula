import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { panelTrailingEdgeGeometry, runCoupled } from '../src/viscous/context.js';
import { createBoundaryLayerAssembly } from '../src/viscous/assembly.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';
import { naca4 } from '../src/geometry/airfoil.js';
const native=JSON.parse(fs.readFileSync('tests/fixtures/fortran/finite-base-panel.json'));
const close=(a,b,tol=2e-12)=>assert.ok(Math.abs(a-b)<=tol*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
const element=surface=>({points:[...structuredClone(surface),{...surface[0]}],trailingEdge:{kind:'finite-base',upperIndex:0,lowerIndex:surface.length-1}});
let system;
function fixture(){
 return system??=createBoundaryLayerAssembly({elements:[element(native.cases[1].points)],alpha:4,mach:0,
  reynolds:1e6,ncrit:9,trips:[1,1],wakeCount:12,wakeLength:1,initialization:'march'});
}

test('geometry-only native spline/TECALC contract matches original Fortran on all three finite-base fixtures',()=>{
 for(const c of native.cases){
  const before=structuredClone(c.points),g=panelTrailingEdgeGeometry(c.points);
  close(g.width,c.expected.normalGap);close(g.tangentialProjection,c.expected.tangentialGap);close(g.magnitude,c.expected.gap);
  assert.equal(g.sharp,c.expected.sharp);
  close(g.upperDerivative.x,c.expected.nodes[0].dxds);close(g.upperDerivative.y,c.expected.nodes[0].dyds);
  close(g.lowerDerivative.x,c.expected.nodes.at(-1).dxds);close(g.lowerDerivative.y,c.expected.nodes.at(-1).dyds);
  assert.deepEqual(c.points,before);g.center.x=999;assert.notEqual(panelTrailingEdgeGeometry(c.points).center.x,999);
 }
 assert.throws(()=>panelTrailingEdgeGeometry([{x:0,y:0},{x:0,y:0},{x:1,y:0}]),/distinct/);
});

test('finite-gap wake and TE kernels replay unchanged original Fortran blocks at Mach0 and .3',()=>{
 const fixture=JSON.parse(fs.readFileSync('tests/fixtures/fortran/integral.json'));
 const cases=fixture.cases.filter(c=>c.input.name.startsWith('finite-'));assert.equal(cases.length,4);
 for(const {input:c,expected} of cases){
  const k=createIntegralKernel({...c.parameters,exactJacobian:true});
  const v=c.regime==='te'?k.trailingEdge(...c.surfaceStates,c.states[1],c.gap):k.interval({regime:c.regime,upstream:c.states[0],downstream:c.states[1]});
  v.residual.forEach((r,i)=>close(r,expected.residual[i],3e-10));
 }
});

test('finite assembly uses total displacement, native gap at fixed wake distance and gap-aware TE rows',()=>{
 const s=fixture(),x=s.initial.slice(),w=s.outer.wakes[0],base=w.body.baseGeometry;
 const model=createXfoilDeadAirGap({normalGap:base.width,upperDerivative:base.upperDerivative,lowerDerivative:base.lowerDerivative,sharp:base.sharp});
 const {states}=s.decode(x);assert.ok(states[w.start].wakeGap>0);assert.equal(states.at(-1).wakeGap,0);
 assert.ok(w.s.some((arc,i)=>i>0&&model.at(arc).gap>0),'fixture resolves the prescribed finite-gap tail');
 for(let j=0;j<w.s.length;j++)close(states[w.start+j].wakeGap,model.at(w.s[j]).gap,1e-16);
 const q=s.residual(x),v=s.kernel.trailingEdge(states[w.body.start],states[w.body.end],states[w.start],base.width);
 for(let k=0;k<3;k++)close(q[4*w.start+k],v.residual[k]*(k===0?20:1/s.thicknessScale),1e-14);
 const oldWrong=s.kernel.trailingEdge(states[w.body.start],states[w.body.end],states[w.start]);
 close((oldWrong.residual[2]-v.residual[2])/s.thicknessScale,base.width/s.thicknessScale,1e-12);
 for(let id=w.start+1;id<=w.end;id++){
  const exact=s.kernel.interval({regime:'wake',upstream:states[id-1],downstream:states[id]});
  exact.residual.forEach((r,k)=>close(q[4*id+k],r*(k===0?20:1),1e-14));
 }
 // A stagnation perturbation changes absolute BL s but cannot change WGAP.
 const t=s.decode(x).st[0],moved=x.slice();moved[4*t.left+3]*=1.002;
 const other=s.decode(moved).states;assert.notEqual(other[w.start].s,states[w.start].s);
 for(let id=w.start;id<=w.end;id++)assert.equal(other[id].wakeGap,states[id].wakeGap);
 const bad=x.slice();bad[4*w.start+2]=bad[4*w.start+1]+.5*base.width/s.thicknessScale;
 assert.equal(s.admissible(bad),false,'positive total H must not conceal negative fluid thickness');
 const d=new Float64Array(x.length);d[4*w.start+2]=-100;
 const trial=s.guardedTrialState(x,d,1),gap=states[w.start].wakeGap/s.thicknessScale;
 close(trial[4*w.start+2]-trial[4*w.start+1]-gap,.2*(x[4*w.start+2]-x[4*w.start+1]-gap));
});

test('native wake seed interpolation preserves fluid DSTR and restores the current prescribed gap once',()=>{
 const s=fixture(),w=s.outer.wakes[0],raw=runCoupled(w.body.points,{alpha:4,reynolds:1e6,ncrit:9,trips:[1,1],maxIterations:1});
 const b=raw.bl,te=b.IBLTE[2],nativeS=[];for(let i=te+1;i<=b.NBL[2];i++)nativeS.push(b.XSSI[i][2]-b.XSSI[te+1][2]);
 for(let j=0;j<w.s.length;j++){
  const id=w.start+j,arc=w.s[j];let k=0;while(k<nativeS.length-2&&nativeS[k+1]<arc)k++;
  const f=Math.min(1,(arc-nativeS[k])/(nativeS[k+1]-nativeS[k]));
  const fluid=(b.DSTR[te+1+k][2]-b.WGAP[k+1])*(1-f)+(b.DSTR[te+2+k][2]-b.WGAP[k+2])*f;
  close(s.initial[4*id+2]*s.thicknessScale-s.wakeGaps[id],fluid,2e-15);
 }
});

test('finite-base warm wake transfer preserves fluid thickness, detached gaps and full solid identity',()=>{
 const s=fixture(),seed=s.exportSeed(s.initial),before=structuredClone(seed),w=s.outer.wakes[0];
 const wakePaths=[structuredClone(w.points)];wakePaths[0][1].y+=.0003;
 const input={elements:[element(native.cases[1].points)],alpha:4,mach:0,reynolds:1e6,ncrit:9,
  trips:[1,1],wakeCount:12,wakeLength:1,initialization:'march',seed,wakePaths};
 const next=createBoundaryLayerAssembly(input);
 assert.notEqual(next.wakeGaps[w.start+1],s.wakeGaps[w.start+1]);
 for(let id=w.start;id<=w.end;id++)close(next.initial[4*id+2]*next.thicknessScale-next.wakeGaps[id],
  s.initial[4*id+2]*s.thicknessScale-s.wakeGaps[id],2e-15);
 assert.deepEqual(seed,before);next.exportSeed(next.initial).wakeGaps.fill(999);assert.notEqual(next.wakeGaps[0],999);
 const missing={...seed};delete missing.wakeGaps;
 assert.throws(()=>createBoundaryLayerAssembly({...input,seed:missing}),/previous prescribed wake gaps/);
 const changed=structuredClone(input.elements);changed[0].points.splice(-1,0,{x:1.0016,y:0});
 assert.throws(()=>createBoundaryLayerAssembly({...input,elements:changed}),/warm-start geometry or conditions changed/);
});

test('complete finite-base BL Jacobian matches independent FD4 including wake, TE, and stagnation chains',t=>{
 const s=fixture(),x=s.initial.slice();s.updateActive(x);const before=x.slice(),n=x.length,j=s.jacobian(x);
 let maximum=0,worst=null;
 for(const phase of [.271,.713]){
  const d=x.map((v,i)=>Math.sin(phase*(i+1))*Math.max(Math.abs(v),i%4===0?.003:.01)),h=2e-6;
  const at=k=>s.residual(x.map((v,i)=>v+k*h*d[i])),[mm,m,p,pp]=[-2,-1,1,2].map(at);
  for(let row=0;row<n;row++)if(row%4!==3){
   let a=0;for(let col=0;col<n;col++)a+=j[row*n+col]*d[col];
   const b=(mm[row]-8*m[row]+8*p[row]-pp[row])/(12*h),e=Math.abs(a-b)/Math.max(1,Math.abs(a),Math.abs(b));
   if(e>maximum){maximum=e;worst={row,analytic:a,fd:b,phase};}
  }
 }
 assert.deepEqual(x,before);assert.ok(maximum<3e-6,JSON.stringify({maximum,worst}));
 t.diagnostic(JSON.stringify({unknowns:n,maximum,worst}));
});

test('omitted sharp path preserves archived initialized states, residuals and all BL Jacobian entries exactly',async()=>{
 const rewrite=(source,dir)=>source.replace(/from\s+(['"])(\.[^'"]+)\1/g,(_,q,p)=>`from ${q}${pathToFileURL(resolve(dir,p)).href}${q}`);
 const oldContext=fs.readFileSync('docs/nlr-panel-bl/src-viscous-context.js.before.txt','utf8');
 const contextURL='data:text/javascript;base64,'+Buffer.from(rewrite(oldContext,'src/viscous')).toString('base64');
 let oldAssembly=fs.readFileSync('docs/nlr-panel-bl/src-viscous-assembly.js.before.txt','utf8');
 oldAssembly=oldAssembly.replace("from './context.js'",`from '${contextURL}'`);
 const old=await import('data:text/javascript;base64,'+Buffer.from(rewrite(oldAssembly,'src/viscous')).toString('base64'));
 const input={elements:[{points:naca4('0012',40)}],alpha:4,trips:[.1,.1],wakeCount:8,initialization:'march'};
 const a=createBoundaryLayerAssembly(input),b=old.createBoundaryLayerAssembly(input);
 assert.deepEqual(a.initial,b.initial);assert.deepEqual(a.decode(a.initial),b.decode(b.initial));
 assert.deepEqual(a.residual(a.initial),b.residual(b.initial));assert.deepEqual(a.jacobian(a.initial),b.jacobian(b.initial));
 assert.equal(a.wakeGaps,undefined);assert.equal(a.baseGeometry,undefined);
 const d=a.initial.map((v,i)=>Math.sin(i*.3));assert.deepEqual(a.guardedTrialState(a.initial,d,.1),b.guardedTrialState(b.initial,d,.1));
});
