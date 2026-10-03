import test from 'node:test';
import assert from 'node:assert/strict';
import {multipoleBasis,createMultipoleFarfield} from '../src/potential/farfield.js';
import {affine,evaluateAffine} from '../src/potential/reconstruction.js';
import {buildMesh} from '../src/euler/mesh.js';
import {createPotentialSystem,solvePotential} from '../src/potential/system.js';

test('MSES farfield expansion differentiates rotation, PG scaling and its second-order circulation term',()=>{
  const settings={center:{x:.25,y:-.1},alpha:13};
  for(const mach of [0,.2,.6])for(const p of [{x:-3,y:1},{x:4,y:-2},{x:-3,y:-1}]){
    const a=multipoleBasis(p,{...settings,mach}),h=1e-5;
    for(const [k,key] of ['x','y'].entries()){
      const plus=multipoleBasis({...p,[key]:p[key]+h},{...settings,mach});
      const minus=multipoleBasis({...p,[key]:p[key]-h},{...settings,mach});
      for(let i=0;i<5;i++)assert.ok(Math.abs(a.velocity[i][k]-(plus.potential[i]-minus.potential[i])/(2*h))<2e-11,`M=${mach} term ${i}, ${key}`);
    }
  }
});

function annulus(angular,radial){
  const vertices=[],cells=[],boundaries=[],id=(i,j)=>(i%angular)*(radial+1)+j;
  for(let i=0;i<angular;i++)for(let j=0;j<=radial;j++){
    const angle=2*Math.PI*i/angular,r=1+2*j/radial;
    vertices.push({x:r*Math.cos(angle),y:r*Math.sin(angle)});
  }
  for(let i=0;i<angular;i++)for(let j=0;j<radial;j++)cells.push([id(i,j),id(i,j+1),id(i+1,j+1),id(i+1,j)]);
  for(let i=0;i<angular;i++)for(const j of [0,radial])boundaries.push({a:id(i,j),b:id(i+1,j),type:j===0?'wall':'farfield'});
  return buildMesh(vertices,cells,boundaries);
}

test('unknown farfield source and doublets recover an independent exterior harmonic solution',()=>{
  // Phi = A log(r) + B x/r² + C y/r², prescribed only by its
  // inner-boundary normal derivative; the outer strengths are unknown.
  const a=.03,b=.06,c=-.04,errors=[];
  const phi=p=>a*Math.log(Math.hypot(p.x,p.y))+(b*p.x+c*p.y)/(p.x*p.x+p.y*p.y);
  const velocity=p=>{
    const r2=p.x*p.x+p.y*p.y,r4=r2*r2;
    return[1+a*p.x/r2+b*(p.y*p.y-p.x*p.x)/r4-2*c*p.x*p.y/r4,
      a*p.y/r2-2*b*p.x*p.y/r4+c*(p.x*p.x-p.y*p.y)/r4];
  };
  for(const [nt,nr] of [[24,4],[48,8],[96,16]]){
    const mesh=annulus(nt,nr),count=mesh.cells.length;
    const ff=createMultipoleFarfield(mesh,{circulation:{centers:[]},offset:count,center:{x:0,y:0}});
    const s=createPotentialSystem(mesh,{sparse:true,boundaryPotential:ff.potential,
      boundaryNormalVelocity:f=>{const [u,v]=velocity(f);return affine(u*f.nx+v*f.ny);},constraints:ff.constraints});
    const r=solvePotential(s);assert.equal(r.converged,true,r.reason);
    const expected=[a,b,c].map(v=>2*Math.PI*v);
    const error=Math.max(...expected.map((v,i)=>Math.abs(r.x[count+i]-v)));errors.push(error);
    assert.ok(r.diagnostics.residual<1e-10);
    assert.ok(Math.abs(r.x[count+3])<1e-14);
    if(nt===96){
      const reference=Float64Array.from(mesh.cells,phi),shift=r.x[0]-reference[0];
      assert.ok(Math.max(...reference.map((v,i)=>Math.abs(r.x[i]-v-shift)))<.002);
    }
  }
  for(let i=1;i<errors.length;i++)assert.ok(errors[i]<.4*errors[i-1],String(errors));
  assert.ok(errors.at(-1)<.002,String(errors));
});

test('farfield potential correction stays continuous across the analytic vortex branch cut',()=>{
  const mesh={cells:Array.from({length:2}),faces:Array.from({length:4},(_,i)=>({boundary:{type:'farfield'}}))};
  const ff=createMultipoleFarfield(mesh,{circulation:{centers:[{x:0,y:0}]},offset:3,center:{x:.25,y:0},mach:.4});
  const x=[0,0,-.2,.003,.02,-.01,.04];
  const a=evaluateAffine(ff.potential({x:-4,y:1e-9}),x),b=evaluateAffine(ff.potential({x:-4,y:-1e-9}),x);
  assert.ok(Math.abs(a-b)<1e-10);
});
