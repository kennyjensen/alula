import test from 'node:test';
import assert from 'node:assert/strict';
import { createDisplacementOperator, polynomialPotential, polynomialSourceBasis } from '../src/inviscid/displacement.js';
import { makePanel } from '../src/inviscid/panel.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { factorLinear, solveLinear } from '../src/numerics/linear.js';

test('reused LU retains lower-triangle pivots for independent right-hand sides',()=>{
  const a=[0,2e-14,3e-14,4e12,5e12,6e12,7,8,10],solve=factorLinear(a);
  for(const b of[[13e-14,32e12,53],[3,2,1],[1,0,0]]){
    const exact=solveLinear(a,b),actual=solve(b);
    actual.forEach((v,i)=>assert.ok(Math.abs(v-exact[i])<1e-12*Math.max(1,Math.abs(v))));
  }
});

test('polynomial source potentials differentiate to independently integrated velocities',()=>{
  const p=makePanel({x:.2,y:-.3},{x:1.1,y:.4}),point={x:.35,y:.8},h=1e-5;
  const b=polynomialSourceBasis(point,p),xp=polynomialPotential({...point,x:point.x+h},p),xm=polynomialPotential({...point,x:point.x-h},p);
  const yp=polynomialPotential({...point,y:point.y+h},p),ym=polynomialPotential({...point,y:point.y-h},p);
  for(let k=0;k<3;k++){
    let u=0,v=0;const n=4000;
    for(let j=0;j<n;j++){
      const xi=(j+.5)/n,dx=point.x-p.a.x-xi*p.length*p.tx,dy=point.y-p.a.y-xi*p.length*p.ty;
      const f=xi**k*p.length/(2*Math.PI*(dx*dx+dy*dy)*n);u+=f*dx;v+=f*dy;
    }
    assert.ok(Math.abs(u-b[k].u)<1e-8);assert.ok(Math.abs(v-b[k].v)<1e-8);
    assert.ok(Math.abs((xp[k].real-xm[k].real)/(2*h)-u)<1e-8);
    assert.ok(Math.abs((yp[k].imag-ym[k].imag)/(2*h)-u)<1e-8);
    assert.ok(Math.abs((yp[k].real-ym[k].real)/(2*h)-v)<1e-8);
    assert.ok(Math.abs((xp[k].imag-xm[k].imag)/(2*h)+v)<1e-8);
  }
});

test('small source panels retain far-field accuracy without polynomial-recurrence cancellation',()=>{
  const p=makePanel({x:.7,y:-.2},{x:.700006,y:-.199992});
  for(const point of[{x:1.3,y:.6},{x:1000,y:-200}]){
    const b=polynomialSourceBasis(point,p),f=polynomialPotential(point,p);
    for(let k=0;k<3;k++){
      let u=0,v=0,real=0,imag=0;const n=10000;
      for(let j=0;j<n;j++){
        const xi=(j+.5)/n,dx=point.x-p.a.x-xi*p.length*p.tx,dy=point.y-p.a.y-xi*p.length*p.ty;
        const weight=xi**k*p.length/(2*Math.PI*n),rr=dx*dx+dy*dy;
        u+=weight*dx/rr;v+=weight*dy/rr;real+=weight*.5*Math.log(rr);imag+=weight*Math.atan2(dy,dx);
      }
      // Midpoint quadrature of xi^2 has O(n^-2) relative error; the analytic
      // formulas should be substantially more accurate than this oracle.
      assert.ok(Math.abs(b[k].u-u)<1e-8*Math.abs(u));assert.ok(Math.abs(b[k].v-v)<1e-8*Math.abs(v));
      assert.ok(Math.abs(f[k].real-real)<1e-8*Math.max(Math.abs(real),p.length));
      assert.ok(Math.abs(f[k].imag-imag)<1e-8*Math.max(Math.abs(imag),p.length));
    }
  }
});

test('a known source-sheet midpoint has zero principal normal velocity after large translations',()=>{
  for(const y of [0,16,256]){
    const p=makePanel({x:1,y},{x:1.0015,y:y+.00013});
    const b=polynomialSourceBasis(p,p,{principal:true,selfMidpoint:true});
    for(let k=0;k<3;k++){
      assert.ok(Math.abs(b[k].u*p.nx+b[k].v*p.ny)<1e-15);
      assert.ok(Math.abs(b[k].u*p.tx+b[k].v*p.ty-(k===0?0:-1/(2*Math.PI)))<1e-15);
    }
  }
});

test('every wake segment conserves displacement flux on nonuniform spacing',()=>{
  const o=createDisplacementOperator({elements:[{points:naca4('0012',40)}],alpha:3,wakeCount:16});
  const mass=Float64Array.from({length:o.total},(_,i)=>.002*Math.sin(i*.31)),sigma=o.sources(mass);
  for(const w of o.wakes)for(let j=0;j<w.points.length-1;j++){
    const row=3*(w.sourceStart+j),h=w.s[j+1]-w.s[j];
    const flux=h*(sigma[row]+sigma[row+1]/2+sigma[row+2]/3);
    assert.ok(Math.abs(flux-(mass[w.start+j+1]-mass[w.start+j]))<1e-15);
    if(j>0){const prev=row-3;assert.ok(Math.abs(sigma[row]-sigma[prev]-sigma[prev+1]-sigma[prev+2])<1e-12);}
  }
});

test('shared displacement and Kutta equations retain influences between elements',()=>{
  const o=createDisplacementOperator({elements:[{points:naca4('0012',40)},{points:transform(naca4('0012',40),{x:1.1,y:-.2,chord:.4})}],alpha:3,wakeCount:16});
  const mass=new Float64Array(o.total),b=o.bodies[1];for(let i=b.start;i<=b.end;i++)mass[i]=.001*Math.sin((i-b.start)*.2);
  const q=o.evaluate(mass);let effect=0;for(let i=0;i<=o.bodies[0].end;i++)effect=Math.max(effect,Math.abs(q[i]-o.q0[i]));
  assert.ok(effect>1e-4&&effect<.1,`cross-element effect ${effect}`);
  for(const b of o.bodies)assert.ok(Math.abs(q[b.start]+q[b.end])<1e-12);
});

test('multielement source angle charts preserve translation, rotation and scale invariance',()=>{
  const elements=[{points:naca4('2412',40)},{points:transform(naca4('0012',40),{x:1.1,y:-.2,chord:.4})}];
  const o=createDisplacementOperator({elements,alpha:3,wakeCount:16});
  const scale=1.7,moved=createDisplacementOperator({elements:elements.map(e=>({points:transform(e.points,{x:3.2,y:-2.4,angle:31,chord:scale})})),alpha:34,wakeCount:16});
  const m=Float64Array.from({length:o.total},(_,i)=>.001*Math.sin(.37*i)),q=o.evaluate(m),qm=moved.evaluate(m.map(v=>v*scale));
  q.forEach((v,i)=>assert.ok(Math.abs(v-qm[i])<1e-8,`${i}: ${v} vs ${qm[i]}`));
});
