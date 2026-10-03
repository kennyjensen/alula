import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { finiteBaseWakePotentialDifference } from '../src/inviscid/finite-base-wake-potential.js';
import { potentialDifference } from '../src/inviscid/streamfunction.js';
import { velocityAt } from '../src/inviscid/linear-vortex.js';
import { makePanel } from '../src/inviscid/panel.js';
const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) <= tolerance, `scalar difference ${a - b} exceeds ${tolerance}`);
const fixture = rounded => {
  const surface = [{x:0,y:.2},{x:-1,y:.2},{x:-1,y:-.2},{x:0,y:-.2}];
  const base = rounded ? [surface.at(-1),{x:.08,y:-.1},{x:.1,y:0},{x:.08,y:.1},surface[0]] : [surface.at(-1),surface[0]];
  return { a:surface[0],field:{u:1,v:.1,gamma:Float64Array.of(.2,.6,.8,-.2),
    panels:surface.slice(1).map((p,i)=>({...makePanel(surface[i],p,0),node:i})),
    basePanels:base.slice(1).map((p,i)=>({...makePanel(base[i],p,0),sourceStrength:.3,vortexStrength:.23,
      cutOrigin:{x:0,y:0},cutDirection:{x:1,y:0}}))} };
};
const pathValue = (path,field) => path.slice(1).reduce((sum,p,i)=>sum+potentialDifference(path[i],p,field),0);
const options = (a,b,field,direction={x:1,y:0})=>({trailingEdge:a,wakeSeed:b,direction,element:0,field});

const saved=JSON.parse(fs.readFileSync(new URL('../docs/nlr-finite-base/finer-potential-path/local-before.json',import.meta.url)));
saved.field.gamma=Float64Array.from(Object.values(saved.field.gamma));
const finer=JSON.parse(fs.readFileSync(new URL('../docs/nlr-finite-base/finer-potential-path/finer-local.json',import.meta.url)));

test('exact NLR64/128 crossing paths continue through the exterior with the same farther-seed gauge',()=>{
  const before=JSON.stringify(saved),far=saved.profiles[0].seed,a=saved.profiles[0].te;
  for(const c of finer.cases.filter(c=>c.body===0)){
    assert.throws(()=>potentialDifference(c.te,c.seed,saved.field),/crosses a finite-base sheet/);
    const r=finiteBaseWakePotentialDifference(options(a,c.seed,saved.field,saved.profiles[0].tangent));
    assert.equal(r.detour,true);assert.equal(r.path.length,4);
    assert.deepEqual(r.path[0],a);assert.deepEqual(r.path.at(-1),c.seed);
    close(r.value,potentialDifference(a,far,saved.field)+potentialDifference(far,c.seed,saved.field));
    const subdivided=r.path.flatMap((p,i)=>i?[{x:.5*(p.x+r.path[i-1].x),y:.5*(p.y+r.path[i-1].y)},p]:[p]);
    close(r.value,pathValue(subdivided,saved.field));
  }
  assert.equal(JSON.stringify(saved),before);
});

test('valid exact NLR32 and straight finite-base direct paths retain bit-identical increments',()=>{
  for(const p of saved.profiles){
    const r=finiteBaseWakePotentialDifference({...options(p.te,p.seed,saved.field,p.tangent),element:p.body});
    assert.equal(r.detour,false);assert.equal(r.value,potentialDifference(p.te,p.seed,saved.field));
  }
  const {a,field}=fixture(false),b={x:.1,y:0};
  const r=finiteBaseWakePotentialDifference(options(a,b,field));
  assert.equal(r.detour,false);assert.equal(r.value,potentialDifference(a,b,field));
});

test('rounded retained cap has path-independent continuation and the physical directional derivative',()=>{
  const {a,field}=fixture(true),b={x:.11,y:.01},far={x:.8,y:.01};
  assert.throws(()=>potentialDifference(a,b,field),/finite-base sheet/);
  const r=finiteBaseWakePotentialDifference(options(a,b,field));assert.equal(r.detour,true);
  close(r.value,potentialDifference(a,far,field)+potentialDifference(far,b,field));
  const h=1e-5,f=x=>finiteBaseWakePotentialDifference(options(a,{x,y:b.y},field)).value;
  const fd=(-f(b.x+2*h)+8*f(b.x+h)-8*f(b.x-h)+f(b.x-2*h))/(12*h);
  close(fd,velocityAt(b,field).u,3e-10);
  // A complete exterior winding keeps the original total sheet circulation.
  const loop=[{x:-2,y:-1},{x:2,y:-1},{x:2,y:1},{x:-2,y:1},{x:-2,y:-1}];
  const gamma=field.panels.reduce((v,p)=>v+.5*p.length*(field.gamma[p.node]+field.gamma[p.node+1]),0)
    +field.basePanels.reduce((v,p)=>v+p.length*p.vortexStrength,0);
  close(pathValue(loop,field),gamma);
});

test('exterior continuation is covariant under rigid motion and physical length scaling',()=>{
  const {a,field}=fixture(true),b={x:.11,y:.01},reference=finiteBaseWakePotentialDifference(options(a,b,field));
  const c=Math.cos(.713),s=Math.sin(.713),L=3.7;
  const rotate=p=>({x:c*p.x-s*p.y,y:s*p.x+c*p.y});
  const transform=p=>{const q=rotate(p);return{x:8+L*q.x,y:-3+L*q.y};};
  const v=rotate({x:field.u,y:field.v}),next={...field,u:v.x,v:v.y,
    panels:field.panels.map(p=>({...makePanel(transform(p.a),transform(p.b),0),node:p.node})),
    basePanels:field.basePanels.map(p=>({...makePanel(transform(p.a),transform(p.b),0),sourceStrength:p.sourceStrength,vortexStrength:p.vortexStrength,
      cutOrigin:transform(p.cutOrigin),cutDirection:rotate(p.cutDirection)}))};
  const result=finiteBaseWakePotentialDifference(options(transform(a),transform(b),next,rotate({x:1,y:0})));
  assert.equal(result.detour,true);close(result.value,L*reference.value,2e-11);
});

test('physical sheet and interior guards remain active; a second body obstructing the continuation is rejected',()=>{
  const {a,field}=fixture(true);
  assert.throws(()=>finiteBaseWakePotentialDifference(options(a,{x:-.5,y:0},field)),/solid body|physical sheet|vortex sheet/);
  assert.throws(()=>potentialDifference({x:.08,y:-.01},{x:.15,y:-.01},field),/finite-base sheet/);
  const b={x:.11,y:.01},good=finiteBaseWakePotentialDifference(options(a,b,field));
  const center={x:.5*(good.path[1].x+good.path[2].x),y:.5*(good.path[1].y+good.path[2].y)};
  const d=.002,obstacle=[{x:center.x-d,y:center.y-d},{x:center.x+d,y:center.y-d},{x:center.x+d,y:center.y+d},{x:center.x-d,y:center.y+d}];
  const blocked={...field,gamma:Float64Array.from([...field.gamma,0,0,0,0,0]),panels:[...field.panels,
    ...obstacle.map((p,i)=>({...makePanel(p,obstacle[(i+1)%4],1),node:4+i}))]};
  assert.throws(()=>finiteBaseWakePotentialDifference(options(a,b,blocked)),/physical sheet|vortex sheet/);
  assert.throws(()=>finiteBaseWakePotentialDifference(options(a,b,field,{x:-1,y:0})),/does not face downstream/);
});
