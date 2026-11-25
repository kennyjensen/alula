// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { naca4, naca4Standard, transform } from '../src/geometry/airfoil.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { createDisplacementOperator } from '../src/inviscid/displacement.js';
import { pressureForces } from '../src/inviscid/pressure-forces.js';
import { solveViscousAssembly } from '../src/viscous/result.js';
import { wakeMomentumDrag } from '../src/viscous/wake-drag.js';
import { createPanelContext } from '../src/viscous/context.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateLeadingTransitionInterval } from '../src/viscous/leading-transition-interval.js';
import { geometricTripArc, tripDistance } from '../src/viscous/trips.js';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';

const exe = fileURLToPath(new URL('./check-displacement', import.meta.url));
const sharp = {points:naca4('2412',48)};
function finite(code, panels, options) {
  const p = transform(naca4Standard(code,panels),options);
  return {points:[...p,p[0]], trailingEdge:{kind:'finite-base',upperIndex:0,lowerIndex:p.length-1}};
}
const main = finite('2412',60,{});
const flap = finite('0012',32,{chord:.3,x:.94,y:-.08,angle:-15});
const nlr = getBenchmarkAirfoil('nlr7301').elements;
function native(elements, alpha, nw, mode, full, paths, mass, parameters=[]) {
  const data = [elements.length];
  for (const e of elements) {
    const te = e.trailingEdge;
    data.push(e.points.length, te ? `F ${te.upperIndex} ${te.lowerIndex}` : 'S');
    for (const p of e.points) data.push(p.x,p.y);
  }
  data.push(alpha,nw,2,mode,+full);
  for (const path of paths ?? []) for (const p of path) data.push(p.x,p.y);
  data.push(...mass,...parameters);
  return spawnSync(exe, {input:data.join(' '), encoding:'utf8', maxBuffer:32*1024*1024});
}

function checkProfiles(ref, rows, trip, reynolds, close, solved=false) {
  const [interval,outer] = rows[0];
  assert.ok(interval >= 0 && interval <= 1e-9);
  if (solved) assert.ok(outer <= 1e-9);
  const regimes=['similarity','laminar','turbulent','wake','transition'];
  const states=rows.slice(1+ref.bodies.length).map(([regime,speed,s,ue,theta,deltaStar,aux,wakeGap])=>
    ({regime:regimes[regime],speed,s,ue,theta,deltaStar,aux,wakeGap}));
  const matched=ref.evaluate(states.map(p=>p.speed*p.deltaStar));
  close(outer,Math.max(...matched.map((v,i)=>Math.abs(v-states[i].speed))),2e-8,'seed edge mismatch');
  const kernel=createIntegralKernel({reynolds,ncrit:9,mach:0,transitionTolerance:1e-12});
  const balance=(block,id)=>block.residual.forEach((v,j)=>close(v,0,2e-7,`BL station ${id}/${j}`));
  ref.bodies.forEach((body,e)=>{
    const [left,right,stagnation,...arcs]=rows[e+1], w=ref.wakes[e];
    const geo=createPanelContext(body.points,0,{geometryOnly:true});
    const ctx={...geo,XSTRIP:[0,...trip]};
    for(const key of ['X','Y','S']) ctx[key]=Float64Array.from([0,...geo[key].subarray(0,geo.N)]);
    const upperS=body.s[left-body.start],lowerS=body.s[right-body.start];
    const q=solved?states.map(s=>s.speed):ref.q0;
    const stag=right>left+1?body.s[left+1-body.start]:upperS-q[left]*(lowerS-upperS)/(q[right]-q[left]);
    close(stagnation,stag,2e-11,'stagnation');
    for(const [side,sign] of [[1,-1],[2,1]]) {
      const arc=geometricTripArc(ctx,side);
      close(arcs[side-1],arc,1e-8,'material trip');
      const tripS=tripDistance({body,side,tripArc:arc},stagnation);
      let previous, transitioned=false;
      for(let id=side===1?left:right;id>=body.start && id<=body.end;id+=sign) {
        const p=states[id];
        close(p.s,Math.abs(body.s[id-body.start]-stagnation),2e-11,'surface distance');
        close(p.speed,sign*p.ue,1e-14,'signed speed');
        assert.ok(p.ue>0 && p.deltaStar>p.theta && p.theta>0);
        if(p.regime==='transition') {
          assert.equal(transitioned,false);
          transitioned=true;
        } else assert.equal(p.regime,previous?(transitioned?'turbulent':'laminar'):'similarity');
        const block=!previous && p.regime==='transition'?
          evaluateLeadingTransitionInterval(kernel,{downstream:p,tripS}):
          kernel.interval({regime:p.regime,upstream:previous,downstream:p,tripS});
        balance(block,`${id} (${p.regime})`);
        previous=p;
      }
      assert.ok(transitioned);
    }
    const base=body.baseGeometry;
    const deadAir=base && createXfoilDeadAirGap({...base,normalGap:base.width});
    for(let id=w.start;id<=w.end;++id) {
      const p=states[id],j=id-w.start;
      assert.equal(p.regime,'wake');
      close(p.s,body.s.at(-1)-stagnation+w.s[j],2e-11,'wake distance');
      close(p.wakeGap,deadAir?.at(w.s[j]).gap??0,2e-11,'wake gap');
      assert.ok(p.ue>0 && p.deltaStar-p.wakeGap>p.theta && p.theta>0);
      balance(j?kernel.interval({regime:'wake',upstream:states[id-1],downstream:p}):
        kernel.trailingEdge(states[body.start],states[body.end],p,base?.width??0),id);
    }
  });
  return states;
}

function check(label,elements,nw,mode,full=true,trip=null,alpha=4,reynolds=1e6) {
  const options={elements,alpha,wakeCount:nw,wakeInitialization:mode===1?'inviscid':'straight'};
  const before=performance.now();
  const seed=createDisplacementOperator(options);
  let ms=performance.now()-before;
  const paths=mode===2?seed.wakes.map(w=>w.points.map((p,j)=>({...p,y:p.y+.01*(1-Math.exp(-w.s[j]))}))):null;
  const pathStart=performance.now();
  const ref=paths?createDisplacementOperator({...options,wakePaths:paths}):seed;
  if (paths) ms=performance.now()-pathStart;
  const mass=Float64Array.from({length:ref.total},(_,i)=>.00012*Math.sin(i*.63)+.00005*Math.cos(i*.11));
  const parameters=trip?[reynolds,9,...elements.flatMap(()=>trip)]:[];
  const out=native(elements,alpha,nw,mode,full,paths,mass,parameters);
  if (process.env.ALULA_PANEL_BENCH) console.log(label,out.stderr.trim(),`JS assembly ${ms.toFixed(2)} ms`);
  assert.equal(out.status,0,`${label} assembly failed ${out.stderr}`);
  const rows=out.stdout.trim().split('\n').map(s=>s.trim().split(/\s+/).map(Number));
  const speed=ref.evaluate(mass),velocity=ref.velocityField(mass);
  const stations=[...ref.bodies,...ref.wakes].flatMap(b=>b.points.map((p,j)=>[p.x,p.y,b.s[j]]));
  const close=(a,b,tol,note)=>assert.ok(Number.isFinite(a)&&Math.abs(a-b)<=tol*(1+Math.abs(b)),
    `${label} ${note}: C ${a}, JS ${b}`);
  assert.equal(rows.length,ref.total+elements.length*(nw-1)+(full?ref.total:0)+(trip?1+elements.length+ref.total:0));
  for (let i=0;i<ref.total;++i) {
    stations[i].forEach((p,k)=>close(rows[i][k],p,2e-11,`station ${i}/${k}`));
    close(rows[i][3],ref.q0[i],2e-9,`inviscid ${i}`);
    close(rows[i][4],speed[i],2e-9,`displaced ${i}`);
  }
  let index=ref.total;
  for (const wake of ref.wakes) for (const p of wake.segments) {
    const v=velocity(p);
    close(rows[index][0],v.u,2e-9,'midpoint u'); close(rows[index++][1],v.v,2e-9,'midpoint v');
  }
  if (full) for (let i=0;i<ref.total;++i) for (let j=0;j<ref.total;++j)
    close(rows[index+i][j],ref.influence[i*ref.total+j],1e-8,`influence ${i}/${j}`);
  if (trip) checkProfiles(ref,rows.slice(index+(full?ref.total:0)),trip,reynolds,close);
  console.log(`Native displacement parity passed for ${label} (${ref.total} stations).`);
  return ref;
}

check('sharp', [sharp],9,0,true,[.08,.13]);
check('trip before first lower station', [sharp],9,0,false,[.08,.005]);
check('inviscid wake', [sharp],9,1,true,[1,1]);
const curve=check('prescribed curved wake', [main,flap],8,2,true,[.15,.18]);
check('mixed trailing edges', [sharp,flap],8,1);
check('NLR 7301 segmented bases', nlr,8,0,false,[1,1],6,2510000);
const paths=curve.wakes.map(w=>w.points.map(p=>({...p})));
paths[0][2]={x:.5,y:.01};
assert.equal(native([main,flap],4,8,2,false,paths,new Float64Array(curve.total)).status,2,
  'A prescribed wake crossing the main element must be rejected.');

const viscousExe=fileURLToPath(new URL('./alula-viscous',import.meta.url));
function solve(elements,alpha,nw,trips,extra=[],progress=true) {
  const input=[elements.length,...elements.flatMap(({points,trailingEdge:te})=>[
    points.length,te?`F ${te.upperIndex} ${te.lowerIndex}`:'S',...points.flatMap(p=>[p.x,p.y])])].join(' ');
  const out=spawnSync(viscousExe,['--topology','--wake-count',String(nw),...(progress?['--progress']:[]),
    ...elements.flatMap(()=>['--trip',...trips.map(String)]),...extra,String(alpha),'1000000','9'],
    {input,encoding:'utf8',maxBuffer:32*1024*1024});
  const events=progress?out.stderr.trim().split('\n').map(line=>JSON.parse(line)):[];
  return {code:out.status,final:JSON.parse(out.stdout),events};
}

for (const [label,elements,alpha,nw,trip] of [
  ['early material trip',[sharp],4,9,[.08,.005]],
  ['symmetric stagnation',[{points:naca4('0012',64)}],0,9,[.15,.15]],
  ['two finite bases',[main,flap],4,8,[.15,.18]]
]) {
  const {code,final,events}=solve(elements,alpha,nw,trip);
  assert.equal(code,0,`${label} failed: ${final.reason}`);
  assert.equal(final.converged,true);
  assert.ok(final.equationResidual<=1e-8 && final.wakeResidual<=1e-6);
  assert.ok(events.some(e=>e.stage==='panel/BL' && e.iteration>0) && events.at(-1).stage==='wake');
  assert.ok(events.every(e=>e.boundaryLayerCurves.length===2*elements.length && Number.isFinite(e.cl)));
  const paths=final.streamlines.map(w=>w.map(([x,y])=>({x,y})));
  const ref=createDisplacementOperator({elements,alpha,wakeCount:nw,wakePaths:paths});
  const mass=new Float64Array(ref.total),params=[1e6,9,...elements.flatMap(()=>trip),80];
  const out=native(elements,alpha,nw,2,false,paths,mass,params);
  assert.equal(out.status,0,out.stderr);
  const rows=out.stdout.trim().split('\n').slice(ref.total+elements.length*(nw-1)).map(s=>s.trim().split(/\s+/).map(Number));
  const close=(a,b,tol,note)=>assert.ok(Math.abs(a-b)<=tol*(1+Math.abs(b)),`${label} ${note}: ${a}, ${b}`);
  const states=checkProfiles(ref,rows,trip,1e6,close,true);
  const field=ref.velocityField(states.map(s=>s.speed*s.deltaStar));
  let normal=0;
  for (const w of ref.wakes) for (const seg of w.segments) {
    const v=field(seg);
    assert.ok(v.u*seg.tx+v.v*seg.ty>0);
    normal=Math.max(normal,Math.abs(v.u*seg.nx+v.v*seg.ny));
  }
  close(final.wakeResidual,normal,4e-8,'wake normal velocity');
  let cl=0,cm=0,cd=0;
  for (const [e,b] of ref.bodies.entries()) {
    const cp=Array.from({length:b.points.length},(_,j)=>1-states[b.start+j].speed**2);
    for (let j=0;j<cp.length;++j) close(final.cp[b.start+j].cp,cp[j],5e-8,'surface pressure');
    let points=b.points;
    if (b.baseGeometry) {
      const n=b.solidPoints.length-1,{lowerIndex:lo,upperIndex:hi}=b.trailingEdge;
      const base=Array.from({length:(hi-lo+n)%n+1},(_,j)=>b.solidPoints[(lo+j)%n]);
      const baseCp=.5*(cp[0]+cp.at(-1));points=[...points,...base];cp.push(...base.map(()=>baseCp));
    }
    const force=pressureForces(points,cp,{alpha});cl+=force.cl;cm+=force.cm;
    cd+=wakeMomentumDrag(states[ref.wakes[e].end],{mach:0}).cd;
  }
  close(final.cl,cl,5e-8,'lift');close(final.cm,cm,5e-8,'pitching moment');close(final.viscousWakeDrag,cd,5e-8,'drag');
  console.log(`Native viscous solve passed for ${label} with wake residual ${normal.toExponential(2)}.`);
}
for (const [flag,expected] of [['--iterations','panel/BL Newton solve stopped'],['--wake-limit','wake iteration limit']]) {
  const result=solve([sharp],4,9,[.08,.13],[flag,'0']);
  assert.equal(result.code,1);
  assert.equal(result.final.reason,expected);
  assert.equal(result.final.converged,false);
}

const nativeSingle=solve([sharp],4,9,[.08,.13]).final;
const webSingle=solveViscousAssembly({elements:[sharp],alpha:4,reynolds:1e6,ncrit:9,trips:[.08,.13],wakeCount:9});
assert.equal(webSingle.status,'solved');
for (const [key,value] of [['cl',webSingle.cl],['cm',webSingle.cm],['viscousWakeDrag',webSingle.cd]])
  assert.ok(Math.abs(nativeSingle[key]-value)<1e-8,`Independent viscous solve differs for ${key}`);
for (const [i,p] of webSingle.elements[0].cp.entries()) assert.ok(Math.abs(nativeSingle.cp[i].cp-p.cp)<1e-7);
if (process.env.ALULA_FINE_COUPLED==='1') {
  const input=[nlr.length,...nlr.flatMap(e=>[e.points.length,'F',e.trailingEdge.upperIndex,
    e.trailingEdge.lowerIndex,...e.points.flatMap(p=>[p.x,p.y])])].join(' ');
  const out=spawnSync(viscousExe,['--topology','--live-progress','6','2510000','9'],
    {input,encoding:'utf8',maxBuffer:32*1024*1024});
  assert.equal(out.status,0,out.stdout.slice(0,400));
  const result=JSON.parse(out.stdout),events=out.stderr.trim().split('\n').map(line=>JSON.parse(line));
  assert.ok(events.some(e=>e.actualAlpha<6));
  assert.equal(result.actualAlpha,6);
  assert.ok(result.converged && result.equationResidual<=1e-8 && result.wakeResidual<=1e-6);
  assert.ok(result.viscousWakeDrag>0 && result.boundaryLayerCurves.length===4);
  console.log('Native incompressible NLR 7301 panel/BL and both wakes converged at 6 degrees.');
}

if (process.env.ALULA_PANEL_BENCH) {
  const elements=[{points:naca4('2412',96)}], trips=[.08,.13];
  const milliseconds=[[],[]];
  for(let run=0;run<=5;++run) {
    const result=[];
    for(const solver of run%2?[0,1]:[1,0]) {
      const start=performance.now();
      result[solver]=solver===0?solve(elements,4,12,trips,[],false).final:
        solveViscousAssembly({elements,alpha:4,trips,reynolds:1e6,ncrit:9,wakeCount:12});
      if(run) milliseconds[solver].push(performance.now()-start);
    }
    const [c,js]=result;
    assert.ok(c.converged && js.status==='solved' && Math.abs(c.cl-js.cl)<1e-7 &&
      Math.abs(c.cm-js.cm)<1e-7 && Math.abs(c.viscousWakeDrag-js.cd)<1e-8);
  }
  const [c,js]=milliseconds.map(samples=>samples.sort((a,b)=>a-b)[2]);
  console.log(`Complete 96-panel viscous solve, 12 wake stations, five-run median: C CLI ${c.toFixed(1)} ms, JS in Node ${js.toFixed(1)} ms.`);
}
