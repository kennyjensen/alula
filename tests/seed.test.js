import test from 'node:test';
import assert from 'node:assert/strict';
import { remapInitialMach,remapSeedMach } from '../src/viscous/seed.js';
import { createIntegralKernel } from '../src/viscous/integral.js';

test('Mach continuation initializes a nearly mixed wake with the same kinematic shape factor',()=>{
  const x=Float64Array.of(.02,.001,.00100001,1),y=remapInitialMach(x,0,.2);
  assert.deepEqual(x,Float64Array.of(.02,.001,.00100001,1));
  const kernel=createIntegralKernel({mach:.2,reynolds:1e6});
  assert.doesNotThrow(()=>kernel.station({s:2,aux:y[0],theta:y[1],deltaStar:y[2],ue:y[3]},'wake'));
  const z=remapInitialMach(y,.2,0);for(let i=0;i<4;i++)assert.ok(Math.abs(z[i]-x[i])<1e-15);
  const seed={key:JSON.stringify({mach:0,alpha:3}),x};const changed=remapSeedMach(seed,.2);
  assert.equal(JSON.parse(changed.key).mach,.2);assert.equal(JSON.parse(seed.key).mach,0);
});
