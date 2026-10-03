import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4 } from '../src/geometry/airfoil.js';
import { createBoundaryLayerAssembly } from '../src/viscous/assembly.js';
import { tripDistance,surfaceTransition } from '../src/viscous/trips.js';

test('a trip before the first BL node uses a physical leading-edge interval with the complete moving-trip Jacobian',()=>{
  for(const mach of [0,.2]){
    const s=createBoundaryLayerAssembly({elements:[{points:naca4('0012',40)}],alpha:.01,mach,
      trips:[.0001,.0001],initialization:'march',wakeCount:12});
    const x=s.initial.slice();s.updateActive(x);const {states,st}=s.decode(x);
    const surface=s.activeSurfaces(st).find(s=>s.transition===0);assert.ok(surface);
    const trip=tripDistance(surface,st[0].s),block=s.leadingTransitionInput(states[surface.ids[0]],trip);
    assert.ok(block.upstream.s>0&&block.upstream.s<trip&&trip<block.downstream.s);
    assert.ok(Math.abs(block.upstream.ue/block.upstream.s-block.downstream.ue/block.downstream.s)<1e-12);
    const tr=surfaceTransition(s,surface,states,st[0].s);
    assert.equal(tr.forced,true);assert.ok(Math.abs(tr.s-trip)<1e-12);
    const j=s.jacobian(x),n=x.length,d=x.map((v,i)=>Math.sin(.37*i)*Math.max(1e-5,Math.abs(v))*.01),h=1e-5;
    const p=s.residual(x.map((v,i)=>v+h*d[i])),m=s.residual(x.map((v,i)=>v-h*d[i]));
    for(let i=0;i<n;i++){
      let value=0;for(let k=0;k<n;k++)value+=j[i*n+k]*d[k];
      assert.ok(Math.abs(value-(p[i]-m[i])/(2*h))/Math.max(1,Math.abs(value))<1e-6,`M=${mach}, row ${i}`);
    }
  }
});
