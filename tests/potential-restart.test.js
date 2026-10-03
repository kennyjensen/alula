import test from 'node:test';
import assert from 'node:assert/strict';
import {assertPotentialMeshSettings,assertPotentialRestart} from '../scripts/validation/potential-restart.js';

const points=[{x:1,y:0},{x:0,y:0},{x:1,y:0}];
const mesh={topology:'constrained-triangular-multielement',contours:[points],vertices:new Array(100),
  controls:{padding:8,surfaceScale:.2,boundaryScale:.2,growth:.35,farSpacing:8/3}};
const controls={mesh:{type:'triangular',padding:8,surfaceScale:.2,growth:.35},farfield:'multipole'};

test('refinement cannot silently reuse a mesh from a different spacing, boundary subdivision or domain',()=>{
  assert.doesNotThrow(()=>assertPotentialMeshSettings(mesh,controls.mesh));
  for(const change of [{surfaceScale:.15},{boundaryScale:.3},{padding:16},{growth:.4},{farSpacing:2}])
    assert.throws(()=>assertPotentialMeshSettings(mesh,{...controls.mesh,...change}),/differs from requested/);
  assert.throws(()=>assertPotentialMeshSettings(mesh,{...controls.mesh,maxVertices:99}),/vertex budget/);
});

test('full-state restarts preserve the requested geometry and wake discretization',()=>{
  const input={elements:[{points}],wakeCount:48,wakeLength:2};
  const saved={config:{input,controls},debug:{originalMesh:mesh},wakePaths:[[{x:1,y:0},{x:3,y:0}]]};
  assert.doesNotThrow(()=>assertPotentialRestart(input,controls,saved));
  for(const change of [{wakeCount:72},{wakeLength:4},{wakePaths:[[{x:1,y:0},{x:5,y:0}]]},
    {elements:[{points:points.map(p=>({...p,x:p.x+1}))}]}])
    assert.throws(()=>assertPotentialRestart({...input,...change},controls,saved));
  assert.throws(()=>assertPotentialRestart(input,{...controls,farfield:'fixed'},saved),/unknown count/);
});
