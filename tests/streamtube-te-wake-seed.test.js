import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { initializeCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-initializer.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/main-flap-te-wake-seed.json', import.meta.url)));
const turn = (a, b, c) => {
  const ux = b.x-a.x, uy = b.y-a.y, vx = c.x-b.x, vy = c.y-b.y;
  return Math.atan2(ux*vy-uy*vx, ux*vx+uy*vy)*180/Math.PI;
};

test('cold main/flap seeds a gradual wake without changing surface profiles or requested conditions', () => {
  const {input, options} = structuredClone(fixture), before = JSON.stringify({input, options});
  const result = initializeCoupledStreamtubeBody(input, options), {system} = result;
  const value = system.evaluate(system.initial), x = system.initial.slice(system.ne);
  assert.equal(result.initialization.wakeInitialization, 'iset-linear-shape');
  assert.equal(result.initialization.thicknessFactor, 1);
  assert.equal(result.initialization.equationsChanged, false);
  assert.equal(result.initialization.flowSolved, false);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(system.conditions.ncrit, 9);
  assert.equal(system.conditions.reynolds, options.reynolds);
  assert.deepEqual(system.bl.trips, options.tripFractions);
  assert.ok(Array.from(value.residual).every(Number.isFinite));
  for (const surface of system.bl.surfaces) for (const id of surface.ids)
    for (let k=0;k<4;k++) assert.equal(x[4*id+k],fixture.originalBL[4*id+k],
      'Changing the cold wake guess must retain every MRCHUE surface variable.');
  const thickness = system.bl.thicknesses(x);
  for (const wake of system.bl.wakes) {
    const first = wake.ids[0], theta = x[4*first+1], count = wake.ids.length-1;
    const h0 = (x[4*first+2] - (system.bl.initialWakeGaps?.[first] ?? 0)/system.bl.scale)/theta;
    for (let k=1; k<wake.ids.length; k++) {
      const id=wake.ids[k];
      assert.equal(x[4*id+1],theta,'Constant momentum thickness in the source-style wake guess.');
      const h=(x[4*id+2]-(system.bl.initialWakeGaps?.[id]??0)/system.bl.scale)/theta;
      assert.ok(Math.abs(h-(h0+(1.1-h0)*k/count))<2e-14);
    }
  }
  const main=input.bodies.findIndex(b=>b.element===0), flap=input.bodies.findIndex(b=>b.element===1);
  const angle=(b,upper)=>{
    const i=input.bodies[b].trailingIndex,g=upper?b+1:b,j=upper?0:value.outer.nodes[g][0].length-1;
    return Math.abs(turn(...[i-1,i,i+1].map(k=>value.outer.nodes[g][k][j])));
  };
  // Measured original startup: main upper50.24°, flap upper70.90°.
  // Keep finite airfoil wedge turns; this is a seed-regression bound,
  // not a new admissibility tolerance for converged physical solutions.
  assert.ok(angle(main,true)<25);
  assert.ok(angle(main,false)<12);
  assert.ok(angle(flap,true)<5);
  assert.ok(angle(flap,false)<8);
  const total=thickness.surfaces[main].upper.at(-1)+thickness.surfaces[main].lower.at(-1);
  assert.ok(1-thickness.wakes[main][0]/total<.02,'No original31% first-interval wake-gap collapse.');
  assert.equal(JSON.stringify({input,options}),before);
});
