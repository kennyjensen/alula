import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { makePanel } from '../src/inviscid/panel.js';
import { viscousObservables } from '../src/viscous/observables.js';
import { provisionalViscousCoefficients } from '../src/viscous/provisional-coefficients.js';

const close=(a,b,t=2e-13)=>assert.ok(Math.abs(a-b)<t,`${a} versus ${b}`);
function fixture({finite=true,count=1,referenceChord=1,normalizedScale=1}={}){
  const surface=[{x:1,y:.005},{x:.7,y:.15},{x:.3,y:.18},{x:.05,y:.08},{x:0,y:0},
    {x:.05,y:-.07},{x:.3,y:-.1},{x:.7,y:-.08},{x:1,y:finite?-.005:.005}];
  const base=[{x:1.0015,y:-.0025},{x:1.002,y:0},{x:1.0015,y:.0025},{...surface[0]}];
  const bodies=[],panels=[],elements=[],states=[],st=[],wakes=[],surfaces=[];
  for(let e=0;e<count;e++){
    const point=p=>({x:(p.x+1.5*e)*normalizedScale,y:(p.y-.2*e)*normalizedScale});
    const points=surface.map(point),solidPoints=finite?[...points,...base.map(point)]:points;
    const start=e*surface.length,first=panels.length,s=[0];
    for(let j=1;j<points.length;j++){const p=makePanel(points[j-1],points[j],e);panels.push(p);s.push(s.at(-1)+p.length);}
    const body={points,start,end:start+points.length-1,first,last:panels.length-1,s,
      ...(finite?{solidPoints,trailingEdge:{kind:'finite-base',upperIndex:0,lowerIndex:8},
        baseGeometry:{width:.01*normalizedScale,center:point({x:1,y:0})}}:{})};
    bodies.push(body);elements.push({name:`Body ${e+1}`,points:solidPoints.map(p=>({x:p.x*referenceChord,y:p.y*referenceChord}))});
    st.push({s:.5*(s[4]+s[5])});
    const q=[-.9,-1.2,-1.05,-.6,-.1,.6,.85,1.1,.9];
    q.forEach((v,j)=>states[start+j]={s:Math.abs(s[j]-st[e].s),ue:Math.abs(v),theta:.001*normalizedScale,deltaStar:.0024*normalizedScale,aux:.03,signedUe:v});
    for(const [side,ids] of [[1,[4,3,2,1,0]],[2,[5,6,7,8]]]){
      const direction=side===1?-1:1,tripArc=st[e].s+direction*.4*(side===1?st[e].s:s.at(-1)-st[e].s);
      surfaces.push({body,side,ids:ids.map(i=>start+i),transition:1,transitionId:start+ids[1],tripArc});
    }
  }
  const surfaceCount=count*surface.length;
  bodies.forEach((body,e)=>{
    const start=surfaceCount+3*e,points=[0,.1,1].map(d=>({x:(1+1.5*e+d)*normalizedScale,y:-.2*e*normalizedScale}));
    wakes.push({body,start,end:start+2,points,s:[0,.1*normalizedScale,normalizedScale]});
    [finite?.01:0,finite?.003:0,0].forEach((gap,j)=>states[start+j]={s:body.s.at(-1)-st[e].s+j*.5*normalizedScale,ue:.9,
      theta:.002*normalizedScale,deltaStar:(.0048+gap)*normalizedScale,aux:.03,...(finite?{wakeGap:gap*normalizedScale}:{})});
  });
  const x=new Float64Array(4*states.length),q0=new Float64Array(states.length);
  states.forEach((s,i)=>{x[4*i+3]=s.signedUe??s.ue;q0[i]=x[4*i+3]*1.03;});
  const system={outer:{bodies,panels,wakes,q0},zeroNodes:[],activeSurfaces:()=>surfaces,
    kernel:{station:(s,regime)=>({cf:regime==='wake'?0:.004,rho:1,hk:(s.deltaStar-(s.wakeGap??0))/s.theta}),
      interval:({tripS})=>({transition:{s:tripS,forced:true}})}};
  return{state:{system,states,st,x},conditions:{elements,alpha:4,mach:0,referenceChord,momentReference:{x:.25*referenceChord,y:0}},surfaces};
}
// Independent midpoint/Gauss edge pressure traction and moment integral.
function forceOracle(f){
  let cx=0,cy=0,cm=0;const alpha=f.conditions.alpha*Math.PI/180;
  const integrate=(a,b,pa,pb)=>{
    const dx=b.x-a.x,dy=b.y-a.y;
    for(const s of [.5-.5/Math.sqrt(3),.5+.5/Math.sqrt(3)]){
      const cp=pa+s*(pb-pa),fx=-.5*cp*dy,fy=.5*cp*dx;
      cx+=fx;cy+=fy;cm-= (a.x+s*dx-.25)*fy-(a.y+s*dy)*fx;
    }
  };
  for(const b of f.state.system.outer.bodies){
    const cp=b.points.map((_,j)=>1-f.state.x[4*(b.start+j)+3]**2);
    for(let j=1;j<b.points.length;j++)integrate(b.points[j-1],b.points[j],cp[j-1],cp[j]);
    if(b.solidPoints){const base=b.solidPoints.slice(b.trailingEdge.lowerIndex),p=.5*(cp[0]+cp.at(-1));
      for(let j=1;j<base.length;j++)integrate(base[j-1],base[j],p,p);}
  }
  return{cl:cy*Math.cos(alpha)-cx*Math.sin(alpha),pressureIntegralDrag:cx*Math.cos(alpha)+cy*Math.sin(alpha),cm};
}

test('Finite-base loads integrate all retained edges with modeled TE pressure and no base skin friction',()=>{
  const f=fixture({count:2}),o=viscousObservables(f.state,f.conditions),expected=forceOracle(f);
  for(const key of ['cl','cm','pressureIntegralDrag'])close(o[key],expected[key]);
  assert.equal(o.surfaces.length,4);assert.equal(o.wakes.length,2);
  o.outputElements.forEach((e,i)=>{
    assert.deepEqual(e.points,f.conditions.elements[i].points);assert.equal(e.cp.length,e.points.length);
    assert.equal(e.finiteBase.points.length,5);assert.match(e.finiteBase.pressureModel,/no base skin friction/);
    close(e.finiteBase.pressure,.19);
  });
  const angle=f.conditions.alpha*Math.PI/180;let cdf=0;
  for(const b of f.state.system.outer.bodies)for(let j=1;j<b.points.length;j++){
    const a=f.state.states[b.start+j-1],z=f.state.states[b.start+j],p=f.state.system.outer.panels[b.first+j-1];
    cdf+=.5*.004*(a.ue**2+z.ue**2)*Math.sign((a.signedUe??a.ue)+(z.signedUe??z.ue))
      *((p.b.x-p.a.x)*Math.cos(angle)+(p.b.y-p.a.y)*Math.sin(angle));
  }
  close(o.cdf,cdf);
});

test('Wake total displacement, effective shape and geometric matching remain distinct, with closed-exit Squire–Young unchanged',()=>{
  const f=fixture(),o=viscousObservables(f.state,f.conditions),w=o.wakes[0],p=w.stations[0];
  close(p.deltaStar,.0148);close(p.wakeGap,.01);close(p.viscousDeltaStar,.0048);close(p.h,2.4);close(p.totalH,7.4);
  close(w.matching.momentum,0);close(w.matching.displacement,0);close(w.matching.massFlux,0);
  close(w.matching.geometricMassFlux,.009);close(w.matching.totalMassFluxDifference,.009);
  close(w.drag,2*.002*.9**3.7);close(o.cd,w.drag);
});

test('Reference-length scaling preserves dimensional output geometry and scales dimensionless loads consistently',()=>{
  const a=fixture(),b=fixture({referenceChord:2,normalizedScale:.5});
  // Same physical moment origin for both reference lengths.
  b.conditions.momentReference={x:.25,y:0};
  const oa=viscousObservables(a.state,a.conditions),ob=viscousObservables(b.state,b.conditions);
  assert.deepEqual(ob.outputElements[0].points,oa.outputElements[0].points);
  close(ob.cl,oa.cl/2);close(ob.cm,oa.cm/4);close(ob.cd,oa.cd/2);close(ob.cdf,oa.cdf/2);
  close(ob.wakes[0].stations[0].wakeGap,oa.wakes[0].stations[0].wakeGap);
  close(ob.wakes[0].stations[0].deltaStar,oa.wakes[0].stations[0].deltaStar);
  assert.deepEqual(ob.surfaces.map(s=>s.transition),oa.surfaces.map(s=>s.transition));
});

test('An exit inside the dead-air gap is explicitly unavailable instead of extrapolating a fictitious viscous shape',()=>{
  const f=fixture(),end=f.state.states.at(-1);end.wakeGap=.001;end.deltaStar+=.001;
  assert.throws(()=>viscousObservables(f.state,f.conditions),e=>e.code==='WAKE_EXIT_GAP_OPEN');
  const p=provisionalViscousCoefficients(f.state,f.conditions);
  assert.equal(p.coefficientStatus,'unavailable');assert.equal(p.coefficientCode,'WAKE_EXIT_GAP_OPEN');
  for(const key of ['cl','cd','cm','cdf'])assert.equal(p[key],null);
});

test('Finite-base reporting preserves the input and works on provisional states; missing geometry is rejected',()=>{
  const f=fixture(),saved={states:structuredClone(f.state.states),x:f.state.x.slice(),elements:structuredClone(f.conditions.elements)};
  const p=provisionalViscousCoefficients(f.state,f.conditions);assert.equal(p.coefficientStatus,'unconverged');
  assert.ok(['cl','cd','cm','cdf'].every(k=>Number.isFinite(p[k])));
  assert.deepEqual(f.state.states,saved.states);assert.deepEqual(f.state.x,saved.x);assert.deepEqual(f.conditions.elements,saved.elements);
  delete f.state.system.outer.bodies[0].solidPoints;
  assert.equal(provisionalViscousCoefficients(f.state,f.conditions).coefficientStatus,'unavailable');
});

test('The complete sharp-only observable and provisional objects remain exactly equal to the archived implementation',async()=>{
  const source=fs.readFileSync('docs/nlr-panel-bl/observables.js.before.txt','utf8')
    .replace(/from '(\.[^']+)'/g,(_,p)=>`from '${pathToFileURL(resolve('src/viscous',p)).href}'`);
  const old=await import('data:text/javascript;base64,'+Buffer.from(source).toString('base64'));
  const f=fixture({finite:false,count:2});
  assert.deepEqual(viscousObservables(f.state,f.conditions),old.viscousObservables(f.state,f.conditions));
  const oldProvisionalSource=fs.readFileSync('docs/nlr-panel-bl/provisional-coefficients.js.before.txt','utf8')
    .replace("'./observables.js'",`'data:text/javascript;base64,${Buffer.from(source).toString('base64')}'`);
  const oldProvisional=await import('data:text/javascript;base64,'+Buffer.from(oldProvisionalSource).toString('base64'));
  assert.deepEqual(provisionalViscousCoefficients(f.state,f.conditions),oldProvisional.provisionalViscousCoefficients(f.state,f.conditions));
});
