import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4 } from '../src/geometry/airfoil.js';
import { createBoundaryLayerAssembly } from '../src/viscous/assembly.js';

test('compressible similarity retains its finite stagnation limit and velocity derivatives near a contour node',()=>{
  const s=createBoundaryLayerAssembly({elements:[{points:naca4('0012',40)}],alpha:2,mach:.2,trips:[.1,.1]});
  const x=s.initial.slice();s.releaseNonzeroStagnation(x);s.updateActive(x);
  const {left,right}=s.decode(x).st[0];
  for(const speed of [1e-5,1e-9,1e-12]){
    x[4*left+3]=-speed;
    const r=s.residual(x),j=s.jacobian(x),n=x.length;
    assert.ok(r.every(Number.isFinite));
    for(const col of [4*left+3,4*right+3]){
      const h=Math.abs(x[col])*1e-3,p=x.slice(),m=x.slice();p[col]+=h;m[col]-=h;
      const rp=s.residual(p),rm=s.residual(m);
      for(const id of [left,right])for(let k=1;k<=2;k++){
        const i=4*id+k,fd=(rp[i]-rm[i])/(2*h),exact=j[i*n+col];
        // The tiny own-velocity perturbation reaches roundoff at 1e-12.
        // Its finite derivative is instead checked with a one-sided step
        // that keeps this node on the same side of stagnation.
        if(col===4*left+3&&speed<1e-5)continue;
        assert.ok(Math.abs(fd-exact)/Math.max(1,Math.abs(exact))<2e-5,`${speed}, row ${i}, col ${col}`);
      }
    }
    if(speed<1e-5){
      const h=1e-6,m=x.slice();m[4*left+3]-=h;const rm=s.residual(m);
      for(const id of [left,right])for(let k=1;k<=2;k++){
        const i=4*id+k,fd=(r[i]-rm[i])/h,exact=j[i*n+4*left+3];
        assert.ok(Math.abs(fd-exact)/Math.max(1,Math.abs(exact))<2e-5);
      }
    }
  }
});
