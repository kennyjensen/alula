import test from 'node:test';
import assert from 'node:assert/strict';
import {interpolatePanelAverages,integratedPanelVelocity} from '../src/potential/edge-sampling.js';
import {affine,evaluateAffine} from '../src/potential/reconstruction.js';
import {channelMesh} from '../src/euler/mesh.js';
import {createPotentialSystem,solvePotential} from '../src/potential/system.js';

function samples(positions,subdivisions,velocityAt,halves){
  return positions.slice(1).flatMap((b,i)=>{
    const a=positions[i],panel={a:{x:a,y:0},length:b-a,tx:1,ty:0};
    const segments=Array.from({length:subdivisions},(_,j)=>({a:{x:a+(b-a)*j/subdivisions,y:0},b:{x:a+(b-a)*(j+1)/subdivisions,y:0}}));
    return Array.from({length:halves},(_,j)=>({s:a+(b-a)*(j+.5)/halves,
      velocity:integratedPanelVelocity(panel,segments,velocityAt,{lo:j/halves,hi:(j+1)/halves})}));
  });
}

test('station-centered averages retain the alternating mode of an exact harmonic outer flow',()=>{
  // phi = a sin(pi*x) exp(-pi*y)/pi gives wall u = 1+a cos(pi*x)
  // and v = -a sin(pi*x), the displacement derivative for m=a cos(pi*x)/pi.
  // Whole-panel means annihilate this mode. Each centered half-panel pair
  // must retain 2/pi of it, including with an odd outer-face subdivision.
  const positions=Array.from({length:9},(_,i)=>i),velocityAt=(_,p)=>[affine(1,[[0,Math.cos(Math.PI*p.x)]]),affine()];
  for(const subdivisions of [4,5]){
    const old=interpolatePanelAverages(samples(positions,subdivisions,velocityAt,1),positions);
    const centered=interpolatePanelAverages(samples(positions,subdivisions,velocityAt,2),positions);
    for(let i=1;i<positions.length-1;i++){
      assert.ok(Math.abs(evaluateAffine(old[i],[.1])-1)<1e-13);
      assert.ok(Math.abs(evaluateAffine(centered[i],[.1])-(1+.1*2/Math.PI*(-1)**i))<1e-8);
    }
  }
});

test('half-panel integration recovers linear nodal velocity on nonuniform intervals and at both TEs',()=>{
  const positions=[0,.001,.01,.07,.13,.45,.7,1],velocityAt=(_,p)=>[affine(2,[[0,p.x]]),affine()];
  for(const subdivisions of [3,4,5]){
    const forms=interpolatePanelAverages(samples(positions,subdivisions,velocityAt,2),positions);
    forms.forEach((f,i)=>assert.ok(Math.abs(evaluateAffine(f,[3])-(2+3*positions[i]))<1e-13));
  }
});

test('the solved harmonic Neumann response converges to its nonzero station-centered velocity',()=>{
  const amplitude=.1,k=Math.PI,expectedAmplitude=2*amplitude/Math.PI,errors=[];
  // Fixed BL stations at integer x. Refine only the outer grid. Normal wall
  // flux is the exact integral of the analytic displacement derivative;
  // potential and all recovered velocities must come from the PDE solve.
  for(const subdivisions of [4,8,16,32]){
    const mesh=channelMesh({nx:8*subdivisions,ny:subdivisions,length:8});
    const system=createPotentialSystem(mesh,{sparse:true,
      boundaryPotential:p=>affine(amplitude/k*Math.sin(k*p.x)*Math.exp(-k*p.y)),
      boundaryNormalVelocity:f=>{
        const a=mesh.vertices[f.a],b=mesh.vertices[f.b],lo=Math.min(a.x,b.x),hi=Math.max(a.x,b.x);
        return affine(-f.ny*amplitude*Math.exp(-k*f.y)*(Math.cos(k*lo)-Math.cos(k*hi))/(k*f.length));
      }});
    const result=solvePotential(system);assert.equal(result.converged,true,result.reason);
    const samples=[];
    for(let j=0;j<8;j++){
      const panel={a:{x:j,y:0},length:1,tx:1,ty:0};
      const segments=mesh.faces.filter(f=>f.boundary?.side==='lower'&&f.x>j&&f.x<j+1)
        .map(f=>({...f,a:mesh.vertices[f.a],b:mesh.vertices[f.b]}));
      for(let half=0;half<2;half++)samples.push({s:j+.25+.5*half,
        velocity:integratedPanelVelocity(panel,segments,(f,p)=>system.reconstruction.velocity(f.owner,p),
          {lo:.5*half,hi:.5*(half+1)})});
    }
    const velocities=interpolatePanelAverages(samples,Array.from({length:9},(_,i)=>i)).map(f=>evaluateAffine(f,result.x));
    // Endpoints use one-sided extrapolation, not centered integration.
    errors.push(Math.max(...velocities.slice(1,-1).map((v,i)=>Math.abs(v-(1+expectedAmplitude*(-1)**(i+1))))));
  }
  for(let i=1;i<errors.length;i++)assert.ok(errors[i]<.35*errors[i-1],String(errors));
  assert.ok(errors.at(-1)<.01*expectedAmplitude,String(errors));
});

test('station-centered recovery retains second-order accuracy on nonuniform quadratic data',()=>{
  const errors=[20,40,80].map(n=>{
    const positions=Array.from({length:n+1},(_,i)=>(i/n)**1.5);
    const forms=interpolatePanelAverages(samples(positions,5,(_,p)=>[affine(p.x*p.x),affine()],2),positions);
    return Math.max(...forms.map((f,i)=>Math.abs(f.constant-positions[i]**2)));
  });
  assert.ok(errors[1]<.28*errors[0]&&errors[2]<.28*errors[1]);
});

test('nonuniform interval averages recover linear edge speed and its state derivative at all nodes including both TEs',()=>{
  const positions=[0,.01,.07,.13,.45,.7,1];
  const samples=positions.slice(1).map((b,i)=>({s:(positions[i]+b)/2,velocity:affine(2,[[0,.5*(positions[i]+b)]])}));
  const forms=interpolatePanelAverages(samples,positions);
  forms.forEach((f,i)=>{assert.ok(Math.abs(evaluateAffine(f,[3])-(2+3*positions[i]))<1e-14);assert.ok(Math.abs((f.coefficients.get(0)??0)-positions[i])<1e-14);});
});

test('nodal velocity recovery from quadratic interval averages has second-order truncation error',()=>{
  const errors=[20,40,80].map(n=>{
    const positions=Array.from({length:n+1},(_,i)=>(i/n)**1.5);
    const samples=positions.slice(1).map((b,i)=>{const a=positions[i];return{s:(a+b)/2,velocity:affine((a*a+a*b+b*b)/3)};});
    return Math.max(...interpolatePanelAverages(samples,positions).map((f,i)=>Math.abs(f.constant-positions[i]**2)));
  });
  assert.ok(errors[1]<.28*errors[0]&&errors[2]<.28*errors[1]);
});
