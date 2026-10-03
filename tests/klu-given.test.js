import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {solveSparseDirect} from '../src/numerics/klu.js';
import {solveSparseDirect as archivedDefault} from '../docs/linear-ordering-integration/before/src/numerics/klu.js';
import {solveSparseDirectGiven as prototype} from '../scripts/validation/klu-given.js';
const solveSparseDirectGiven=(a,b,options)=>solveSparseDirect(a,b,{ordering:'given',...options});
import {solveLinear} from '../src/numerics/linear.js';

const csr = rows => ({n:rows.length,rowPtr:Int32Array.from({length:rows.length+1},(_,i)=>i*rows.length),
  colIndex:Int32Array.from({length:rows.length**2},(_,i)=>i%rows.length),values:Float64Array.from(rows.flat())});
const product = (a,x) => Float64Array.from({length:a.n},(_,i)=>{
  let sum=0;for(let k=a.rowPtr[i];k<a.rowPtr[i+1];k++)sum+=a.values[k]*x[a.colIndex[k]];return sum;});
const identity=n=>Int32Array.from({length:n},(_,i)=>i);
const relative = (a,b,x)=>Math.hypot(...product(a,x).map((v,i)=>v-b[i]))/Math.hypot(...b);
const options=(p,q,extra={})=>({rowPermutation:p,columnPermutation:q,pivotTolerance:1,btf:false,...extra});
const hash=p=>createHash('sha256').update(fs.readFileSync(p)).digest('hex');

test('independent nonidentity P and Q preserve original RHS/output on a nonsymmetric pivoting system',()=>{
  const rows=[[0,2,0,-1,0],[3,0,5,0,0],[0,-2,0,4,1],[1,0,-3,0,2],[0,1,2,0,-1]],a=csr(rows);
  const x=Float64Array.of(1.25,-.75,2.5,.125,-3),b=product(a,x),p=Int32Array.of(3,0,4,1,2),q=Int32Array.of(2,4,1,0,3);
  const before=structuredClone({a,b,p,q});
  const r=solveSparseDirectGiven(a,b,options(p,q));
  assert(Math.max(...r.x.map((v,i)=>Math.abs(v-x[i])))<1e-12);
  const dense=solveLinear(Float64Array.from(rows.flat()),b);
  assert(Math.max(...r.x.map((v,i)=>Math.abs(v-dense[i])))<1e-12);
  assert(relative(a,b,r.x)<1e-12);assert.equal(r.attempts.length,1);
  assert.equal(r.btf,false);assert.equal(r.symbolicBlocks,1);assert.equal(r.ordering,'given');assert.equal(r.pivotTolerance,1);
  assert.deepEqual({a,b,p,q},before);
});

test('explicit A(P,Q), b(P) solve agrees with KLU internal permutation and inverse Q mapping',()=>{
  const rows=[[1e-12,2,-1,0],[3,-2,0,1],[0,5,4,-1],[2,0,1,3]],a=csr(rows),b=Float64Array.of(.7,-2.1,3.4,.33);
  const p=Int32Array.of(2,0,3,1),q=Int32Array.of(1,3,0,2);
  const r=solveSparseDirectGiven(a,b,options(p,q));
  const ordered=csr(Array.from(p,i=>Array.from(q,j=>rows[i][j]))),rhs=Float64Array.from(p,i=>b[i]);
  const y=solveSparseDirectGiven(ordered,rhs,options(identity(4),identity(4))).x,back=new Float64Array(4);
  q.forEach((original,k)=>{back[original]=y[k];});
  assert(Math.max(...back.map((v,i)=>Math.abs(v-r.x[i])))<2e-12);
  assert(relative(a,b,r.x)<1e-12);
});

test('scaled nonsymmetric adversarial matrix with tiny diagonal passes original CSR gate and forward check',()=>{
  const n=31,rows=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>
    (j===i?1e-12:[(i+1)%n,(i+3)%n,(i+19)%n].includes(j)?Math.sin(.3*i+.7*j):0)*2**(i%13-6)));
  const a=csr(rows),exact=Float64Array.from({length:n},(_,i)=>Math.cos(.17*i)),b=product(a,exact);
  const p=Int32Array.from({length:n},(_,i)=>(7*i+3)%n),q=Int32Array.from({length:n},(_,i)=>(11*i+5)%n),before=structuredClone({a,b,p,q});
  const r=solveSparseDirectGiven(a,b,options(p,q));
  assert(r.relativeResidual<=1e-10);assert(relative(a,b,r.x)<=1e-10);
  assert(Math.max(...r.x.map((v,i)=>Math.abs(v-exact[i])))<1e-10);
  assert.deepEqual({a,b,p,q},before);assert.equal(r.attempts.length,1);
});

test('BTF false is honored by symbolic kernel and true is explicit, never an automatic retry',()=>{
  const a=csr([[2,1,0,0],[0,3,-1,0],[0,0,4,2],[0,0,0,5]]),b=Float64Array.of(1,2,3,4),o=options(identity(4),identity(4));
  const off=solveSparseDirectGiven(a,b,o),on=solveSparseDirectGiven(a,b,{...o,btf:true});
  assert.equal(off.btf,false);assert.equal(off.symbolicBlocks,1);
  assert.equal(on.btf,true);assert.equal(on.symbolicBlocks,4);
  assert.deepEqual(off.x,on.x);assert.equal(off.attempts.length,1);assert.equal(on.attempts.length,1);
});

test('singular systems fail once, including zero RHS, and cleanup permits a subsequent solve',()=>{
  const a=csr([[1,2,3],[2,4,6],[3,6,9]]),o=options(Int32Array.of(2,0,1),Int32Array.of(1,2,0));
  for(const b of [Float64Array.of(1,0,0),new Float64Array(3)])
    assert.throws(()=>solveSparseDirectGiven(a,b,o),e=>e.code==='KLU_GIVEN_FACTORIZATION'&&e.status===1&&e.attempts.length===1);
  const valid=csr([[0,2,0],[3,0,0],[0,0,4]]),r=solveSparseDirectGiven(valid,Float64Array.of(4,9,16),o);
  assert.deepEqual(Array.from(r.x),[3,2,4]);
  assert.deepEqual(solveSparseDirectGiven(valid,new Float64Array(3),o).x,new Float64Array(3));
});

test('permutation/control/CSR validation rejects malformed requests without silently falling back',()=>{
  const a=csr([[2,1],[3,4]]),b=Float64Array.of(1,2),o=options(identity(2),identity(2));
  for(const key of ['rowPermutation','columnPermutation'])for(const bad of [undefined,[0],[0,0],[0,2],[-1,1],[0,.5],[0,NaN]])
    assert.throws(()=>solveSparseDirectGiven(a,b,{...o,[key]:bad}),/bijection/);
  for(const bad of [{btf:1},{maxRefinements:5},{maxRefinements:-1},{pivotTolerance:0},{pivotTolerance:1.1},{tolerance:0},{ordering:'amd'}])
    assert.throws(()=>solveSparseDirectGiven(a,b,{...o,...bad}),/controls/);
  assert.throws(()=>solveSparseDirectGiven({...a,colIndex:Int32Array.of(1,0,0,1)},b,o),/matrix entry/);
  assert.throws(()=>solveSparseDirectGiven(a,Float64Array.of(1,Infinity),o),/right-hand side/);
});

test('impossible precision fails the unchanged original-system gate after one configured factor',()=>{
  const a=csr([[.11,.27,-.31],[.23,-.17,.71],[.61,.19,.43]]),b=Float64Array.of(.73,-.91,.37),before=structuredClone({a,b});
  assert.throws(()=>solveSparseDirectGiven(a,b,options(Int32Array.of(2,0,1),Int32Array.of(1,2,0),{tolerance:1e-30,maxRefinements:0})),
    e=>e.code==='KLU_RESIDUAL_LIMIT'&&e.relativeResidual>1e-30&&e.attempts.length===1&&e.attempts[0].refinements===0);
  assert.deepEqual({a,b},before);
});

test('given path retains prototype factors and direction while certifying the compensated residual',()=>{
  const n=17,rows=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>
    j===i?1e-10:[(i+1)%n,(i+3)%n,(i+11)%n].includes(j)?Math.cos(.31*i+.57*j):0));
  const a=csr(rows),b=Float64Array.from({length:n},(_,i)=>Math.sin(.13*i));
  const p=Int32Array.from({length:n},(_,i)=>(7*i+1)%n),q=Int32Array.from({length:n},(_,i)=>(5*i+3)%n);
  for(const pivotTolerance of [.001,1]){
    const controls=options(p,q,{pivotTolerance}),old=prototype(a,b,controls),current=solveSparseDirectGiven(a,b,controls);
    for(const key of ['x','refinements','factorNonzeros','symbolicBlocks','ordering','pivotTolerance','btf'])
      assert.deepEqual(current[key],old[key],key);
    assert(current.relativeResidual<=1e-10);assert(relative(a,b,current.x)<=1e-10);
    assert.equal(current.residualEvaluation,'compensated-original-system');
    const withoutResidual=attempts=>attempts.map(({relativeResidual,...rest})=>rest);
    assert.deepEqual(withoutResidual(current.attempts),withoutResidual(old.attempts));
    assert.equal(current.backend,'klu-wasm');assert.equal(current.attempts.length,1);
  }
});

test('default and explicit AMD/COLAMD retain archived factors and IEEE directions with accurate residual metadata',()=>{
  const n=31,rows=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>
    j===i?1e-9:[(i+1)%n,(i+3)%n,(i+19)%n].includes(j)?Math.sin(.3*i+.7*j):0));
  const a=csr(rows),b=Float64Array.from({length:n},(_,i)=>Math.cos(i));
  for(const controls of [{},{ordering:'amd',pivotFallback:false},{ordering:'colamd',pivotTolerance:1,pivotFallback:false},{preferredOrdering:'colamd'}]){
    const current=solveSparseDirect(a,b,controls),old=archivedDefault(a,b,controls);
    assert.equal(current.residualEvaluation,'compensated-original-system');
    assert(current.relativeResidual<=1e-10);assert(relative(a,b,current.x)<=1e-10);
    const withoutResidual=({relativeResidual,residualEvaluation,attempts,...rest})=>({...rest,
      attempts:attempts.map(({relativeResidual,...attempt})=>attempt)});
    assert.deepEqual(withoutResidual(current),withoutResidual(old));
  }
  assert.deepEqual(solveSparseDirect(a,new Float64Array(n)),archivedDefault(a,new Float64Array(n)));
  const singular=csr([[1,1],[2,2]]),rhs=Float64Array.of(1,0);
  const fail=f=>{try{f(singular,rhs);assert.fail('Expected singular error.');}catch(e){return{message:e.message,attempts:e.attempts};}};
  assert.deepEqual(fail(solveSparseDirect),fail(archivedDefault));
});

test('given ordering never falls back and non-given modes reject permutation/BTF controls',()=>{
  const a=csr([[.11,.27,-.31],[.23,-.17,.71],[.61,.19,.43]]),b=Float64Array.of(.73,-.91,.37);
  const controls=options(Int32Array.of(2,0,1),Int32Array.of(1,2,0),{tolerance:1e-30,maxRefinements:0,pivotTolerance:.001,pivotFallback:true});
  assert.throws(()=>solveSparseDirectGiven(a,b,controls),e=>e.code==='KLU_RESIDUAL_LIMIT'&&e.attempts.length===1&&e.attempts[0].pivotTolerance===.001);
  for(const bad of [{btf:false},{rowPermutation:identity(3)},{columnPermutation:identity(3)}])
    assert.throws(()=>solveSparseDirect(a,b,bad),/controls/);
});

test('new binary preserves old ABI and immutable before archive, with unmodified upstream sources',async()=>{
  const base='third_party/klu/',build=JSON.parse(fs.readFileSync(base+'BUILD.json'));
  for(const [file,key]of [['klu.wasm','wasmSha256'],['bridge.c','bridgeSha256'],['build.js','buildScriptSha256']])assert.equal(hash(base+file),build[key]);
  const provenance=JSON.parse(fs.readFileSync(base+'PROVENANCE.json'));
  for(const [p,expected]of Object.entries(provenance.sourceHashes))assert.equal(hash(base+p),expected,p);
  const archive='docs/linear-ordering-integration/before/',before=JSON.parse(fs.readFileSync(archive+'manifest.json'));
  for(const [p,expected]of Object.entries(before.sha256))assert.equal(hash(archive+p),expected,p);
  const m=await WebAssembly.compile(fs.readFileSync(base+'klu.wasm'));
  const old=await WebAssembly.compile(fs.readFileSync(archive+'src/numerics/vendor/klu/klu.wasm'));
  assert.deepEqual(WebAssembly.Module.imports(m),[]);
  const exports=WebAssembly.Module.exports(m).map(e=>e.name);
  for(const e of WebAssembly.Module.exports(old))assert(exports.includes(e.name));
  assert(exports.includes('mses_klu_factor_given'));
});
