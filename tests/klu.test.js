import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {solveSparseDirect,solveSparseDirectAligned} from '../src/numerics/klu.js';
import {sparseMatrix,sparseProduct,sparseDense,solveSparse} from '../src/numerics/sparse.js';
import {solveLinear} from '../src/numerics/linear.js';

test('WASM sparse LU agrees with independent dense LU on a nonsymmetric matrix requiring pivoting',()=>{
  const n=31,a=sparseMatrix(Array.from({length:n},(_,i)=>[i,(i+1)%n,(i+3)%n,(i+19)%n]));
  for(let i=0;i<n;i++)for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)a.values[p]=a.colIndex[p]===i?1e-9:Math.sin(.3*i+.7*a.colIndex[p]);
  const b=Float64Array.from({length:n},(_,i)=>Math.cos(i)),expected=solveLinear(sparseDense(a),b);
  for(let repeat=0;repeat<3;repeat++){
    const r=solveSparseDirect(a,b);assert.ok(r.relativeResidual<1e-10);
    for(let i=0;i<n;i++)assert.ok(Math.abs(r.x[i]-expected[i])<2e-10*Math.max(1,Math.abs(expected[i])));
  }
});

test('WASM sparse LU handles a scaled nonsymmetric PDE stencil and preserves the original matrix',()=>{
  const width=25,n=width**2,a=sparseMatrix(Array.from({length:n},(_,i)=>[i,...(i%width?[i-1]:[]),...(i%width<width-1?[i+1]:[]),...(i>=width?[i-width]:[]),...(i<n-width?[i+width]:[])]));
  for(let i=0;i<n;i++)for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)a.values[p]=(a.colIndex[p]===i?4.1:a.colIndex[p]===i-1?-1.2:a.colIndex[p]===i+1?-.8:-1)*10**(i%13-6);
  const original=a.values.slice(),exact=Float64Array.from({length:n},(_,i)=>Math.sin(.37*i)),b=sparseProduct(a,exact),rhs=b.slice();
  const r=solveSparseDirect(a,b),iterative=solveSparse(a,b);
  assert.ok(r.relativeResidual<1e-10);assert.ok(r.factorNonzeros>n);
  assert.ok(Math.max(...r.x.map((v,i)=>Math.abs(v-exact[i])))<2e-12);
  assert.ok(Math.max(...r.x.map((v,i)=>Math.abs(v-iterative.x[i])))<2e-8);
  assert.deepEqual(a.values,original);assert.deepEqual(b,rhs);
});

test('WASM sparse LU refuses singular, malformed, or nonfinite systems and recovers for the next solve',()=>{
  const a=sparseMatrix([[0,1],[0,1]]);a.values.set([1,1,2,2]);
  assert.throws(()=>solveSparseDirect(a,Float64Array.of(1,0)),/factorization failed/);
  a.values[0]=NaN;assert.throws(()=>solveSparseDirect(a,Float64Array.of(1,0)),/matrix entry/);
  a.values.set([0,2,3,0]);
  assert.deepEqual(Array.from(solveSparseDirect(a,Float64Array.of(4,9)).x),[3,2]);
  assert.throws(()=>solveSparseDirect(a,Float64Array.of(1,Infinity)),/right-hand side/);
  assert.deepEqual(solveSparseDirect(a,new Float64Array(2)).x,new Float64Array(2));
  assert.throws(()=>solveSparseDirect(a,Float64Array.of(1,0),{ordering:'invalid'}),/controls/);
  for(const pivotTolerance of [0,-1,1.01,NaN,Infinity])
    assert.throws(()=>solveSparseDirect(a,Float64Array.of(1,0),{pivotTolerance}),/controls/);
  assert.throws(()=>solveSparseDirect(a,Float64Array.of(1,0),{pivotFallback:'yes'}),/controls/);
});

test('failed AMD factor accuracy retries COLAMD on the same frozen refined coupled matrix',t=>{
  const fixture=JSON.parse(readFileSync(new URL('./fixtures/klu-refined-coupled.json',import.meta.url)));
  const a={...fixture.matrix,rowPtr:Int32Array.from(fixture.matrix.rowPtr),colIndex:Int32Array.from(fixture.matrix.colIndex),values:Float64Array.from(fixture.matrix.values)};
  const b=Float64Array.from(fixture.rhs),values=a.values.slice(),rhs=b.slice();
  const r=solveSparseDirect(a,b);
  assert.equal(r.ordering,'colamd');assert.equal(r.attempts.length,2);
  assert.ok(r.attempts[0].relativeResidual>1e-4);
  assert.ok(r.relativeResidual<1e-10);
  const ax=sparseProduct(a,r.x),norm=v=>Math.hypot(...v);
  assert.ok(norm(b.map((v,i)=>v-ax[i]))/norm(b)<1e-10);
  assert.deepEqual(a.values,values);assert.deepEqual(b,rhs);
  // A known solution supplies an independent forward-error check without
  // requiring a dense 7919-by-7919 factorization or another airfoil solve.
  const exact=Float64Array.from({length:a.n},(_,i)=>Math.sin(.37*i));
  const manufactured=solveSparseDirect(a,sparseProduct(a,exact),{preferredOrdering:r.ordering});
  assert.equal(manufactured.attempts.length,1);
  const error=Math.max(...exact.map((v,i)=>Math.abs(v-manufactured.x[i])));
  assert.ok(error<2e-7,`Forward error ${error}`);
  t.diagnostic(JSON.stringify({attempts:r.attempts,manufacturedForwardError:error}));
});

test('shipped WASM has no host imports and matches its source/build provenance',async()=>{
  const base=new URL('../third_party/klu/',import.meta.url),build=JSON.parse(readFileSync(new URL('BUILD.json',base)));
  const hash=path=>createHash('sha256').update(readFileSync(new URL(path,base))).digest('hex');
  assert.equal(hash('klu.wasm'),build.wasmSha256);assert.equal(hash('bridge.c'),build.bridgeSha256);assert.equal(hash('build.js'),build.buildScriptSha256);
  const provenance=JSON.parse(readFileSync(new URL('PROVENANCE.json',base)));
  for(const [path,expected]of Object.entries(provenance.sourceHashes))assert.equal(hash(path),expected,path);
  const module=await WebAssembly.compile(readFileSync(new URL('klu.wasm',base)));
  assert.deepEqual(WebAssembly.Module.imports(module),[]);
});

test('an impossible residual target retains accurate residual and every factor attempt in failure diagnostics',()=>{
  const a={n:3,rowPtr:Int32Array.of(0,3,6,9),colIndex:Int32Array.of(0,1,2,0,1,2,0,1,2),
    values:Float64Array.of(.11,.27,-.31,.23,-.17,.71,.61,.19,.43)},b=Float64Array.of(.73,-.91,.37);
  assert.throws(()=>solveSparseDirect(a,b,{ordering:'amd',pivotFallback:false,maxRefinements:0,tolerance:1e-30}),error=>{
    assert.equal(error.code,'KLU_RESIDUAL_LIMIT');
    assert.equal(error.diagnostics.relativeResidual,error.relativeResidual);
    assert.equal(error.diagnostics.tolerance,1e-30);
    assert.equal(error.diagnostics.residualEvaluation,'compensated-original-system');
    assert.deepEqual(error.diagnostics.attempts,error.attempts);
    assert.notEqual(error.diagnostics.attempts,error.attempts);
    assert(error.relativeResidual>1e-30);return true;
  });
});

test('equation alignment preserves the original solution, residual gate and inputs',()=>{
  const a=sparseMatrix([[0,1],[1,2],[0,2]]);a.values.set([2,3,4,5,6,7]);
  const exact=Float64Array.of(1,-2,3),b=sparseProduct(a,exact),before=structuredClone(a),rhs=b.slice();
  const r=solveSparseDirectAligned(a,b,[2,0,1]);
  assert.equal(r.equationOrdering,'aligned');assert.ok(r.relativeResidual<=1e-10);
  for(let i=0;i<3;i++)assert.ok(Math.abs(r.x[i]-exact[i])<1e-12);
  assert.deepEqual(a,before);assert.deepEqual(b,rhs);
  assert.throws(()=>solveSparseDirectAligned(a,b,[0,0,2]),/bijection/);
  assert.throws(()=>solveSparseDirectAligned(a,b,[0,1]),/bijection/);
  for(const corrupt of [m=>{m.rowPtr[1]=-1;},m=>{m.colIndex[0]=3;},m=>{m.values[0]=NaN;}]){
    const bad=structuredClone(a);corrupt(bad);
    assert.throws(()=>solveSparseDirectAligned(bad,b,[2,0,1]),/Invalid sparse/);
  }
  assert.deepEqual(solveSparseDirectAligned(a,new Float64Array(3),[2,0,1]).x,new Float64Array(3));
  const singular=sparseMatrix([[0,1],[0,1]]);singular.values.set([1,1,2,2]);
  assert.throws(()=>solveSparseDirectAligned(singular,Float64Array.of(1,0),[1,0]),/factorization failed/);
});
