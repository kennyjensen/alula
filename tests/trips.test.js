import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform } from '../src/geometry/airfoil.js';
import { runCoupled } from '../src/viscous/context.js';
import { geometricTripArc,tripDistance } from '../src/viscous/trips.js';
import { xifset } from '../src/viscous/xfoil/xbl.js';

test('material trip coordinates survive an upstream/downstream stagnation crossing and match native XIFSET',()=>{
  const points=transform(naca4('2412',80),{chord:.4,angle:17,x:.9,y:-.2});
  const {bl}=runCoupled(points,{alpha:3,trips:[.05,.1],maxIterations:1});
  const body={s:Array.from(bl.S.subarray(1,bl.N+1))};
  for(const side of [1,2]){
    const arc=geometricTripArc(bl,side),surface={side,body,tripArc:arc};
    for(const displacement of [-.01,.01,-.005]){
      bl.SST=arc+displacement;
      bl.XSSI[bl.IBLTE[1]][1]=bl.SST-body.s[0];
      bl.XSSI[bl.IBLTE[2]][2]=body.s.at(-1)-bl.SST;
      xifset(bl,side);
      assert.ok(Math.abs(tripDistance(surface,bl.SST)-bl.XIFORC)<1e-13);
      assert.equal(geometricTripArc(bl,side),arc,'a physical trip cannot move or become the TE when stagnation crosses it');
    }
  }
});
