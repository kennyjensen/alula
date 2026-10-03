import test from 'node:test';
import assert from 'node:assert/strict';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { distributeSourceSurface } from '../src/geometry/curvature-surface-spacing.js';
const points = [[1,0],[.8,.05],[.2,.15],[0,0],[.1,-.12],[.65,-.08],[1,0]].map(([x,y])=>({x,y}));
const curve = createContourCurve(points), stagnation = curve.knots[3];
const close = (a,b) => assert.ok(Math.abs(a-b) < 1e-13, `${a} != ${b}`);
for (const side of ['upper','lower']) {
  test(`source ${side} spacing preserves original nodes at native branch count`,()=>{
    const f = distributeSourceSurface({curve,side,stagnation,count:4});
    const expected = side === 'upper' ? points.slice(0,4).reverse() : points.slice(3);
    f.forEach((fraction,i)=>{const p=curve.branch(side,fraction,stagnation).point;close(p.x,expected[i].x);close(p.y,expected[i].y);});
  });
  test(`source ${side} refinement does not add quadratic stagnation clustering`,()=>{
    const coarse=distributeSourceSurface({curve,side,stagnation,count:7});
    const fine=distributeSourceSurface({curve,side,stagnation,count:13});
    close(fine[1],coarse[1]/2);
    coarse.forEach((f,i)=>close(f,fine[2*i]));
    assert.equal(fine[0],0);assert.equal(fine.at(-1),1);
    assert.ok(fine.every((f,i)=>!i || f>fine[i-1]));
  });
}
test('source spacing handles an interior stagnation point and rigid transforms without changing geometry',()=>{
  const before=structuredClone(points), angle=.4,scale=2.7;
  const transformed=createContourCurve(points.map(p=>({x:scale*(p.x*Math.cos(angle)-p.y*Math.sin(angle))+3,y:scale*(p.x*Math.sin(angle)+p.y*Math.cos(angle))-2})));
  const stag=(curve.knots[2]+curve.knots[3])/2;
  for(const side of ['upper','lower']) {
    const a=distributeSourceSurface({curve,side,stagnation:stag,count:129});
    const b=distributeSourceSurface({curve:transformed,side,stagnation:stag*scale,count:129});
    a.forEach((v,i)=>close(v,b[i]));
    assert.ok(a.every((v,i)=>Number.isFinite(v)&&(!i||v>a[i-1])));
  }
  assert.deepEqual(points,before);
});
test('source spacing rejects degenerate source coordinates and invalid requests',()=>{
  const valid={curve,side:'upper',stagnation,count:5};
  for(const change of [{side:'wake'},{stagnation:0},{count:2},{curve:{...curve,knots:[0,1,1,curve.length]}}])
    assert.throws(()=>distributeSourceSurface({...valid,...change}),/Invalid source/);
});
