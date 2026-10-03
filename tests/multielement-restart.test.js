import test from 'node:test';
import assert from 'node:assert/strict';
import {naca4} from '../src/geometry/airfoil.js';
import {createCoupledAssembly} from '../src/viscous/assembly.js';
import {solveMultielementViscous} from '../src/viscous/multielement.js';
const input=()=>({elements:[{points:naca4('0012',40)}],alpha:4,wakeCount:8,initialization:'march',wakeInitialization:'inviscid'});
test('public wake-coupling wrapper retains supplied seed instead of silently cold-starting',()=>{
 const data=input(),system=createCoupledAssembly(data),seed=system.exportSeed(system.initial);
 seed.x[5]*=1.001;
 const raw=solveMultielementViscous({...data,seed},{maxIterations:0,maxWakeIterations:0});
 assert.deepEqual(raw.system.initial,seed.x);
 assert.equal(raw.history.filter(h=>h.iteration>0).length,0);
});
test('public wake-coupling wrapper enforces the saved seed geometry/condition identity',()=>{
 const data=input(),system=createCoupledAssembly(data),seed=system.exportSeed(system.initial);
 assert.throws(()=>solveMultielementViscous({...data,alpha:5,seed},{maxIterations:0,maxWakeIterations:0}),/warm-start geometry or conditions changed/);
});
