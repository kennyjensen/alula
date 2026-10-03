// SPDX-License-Identifier: GPL-2.0-or-later
// Pivoted sparse LU using the unmodified SuiteSparse KLU kernel in WebAssembly.
// Every accepted solution is checked against the original JavaScript matrix.
import { compensatedSparseResidual } from './compensated-sparse-residual.js';
const url=new URL('../../third_party/klu/klu.wasm',import.meta.url);
const bytes=url.protocol==='file:'
  ?await(await import('node:fs/promises')).readFile(url)
  :await fetch(url).then(response=>{if(!response.ok)throw new Error(`Cannot load sparse solver: HTTP ${response.status}.`);return response.arrayBuffer();});
const {instance}=await WebAssembly.instantiate(bytes,{}),wasm=instance.exports;
wasm._initialize();

function norm(values){
  let scale=0,sum=0;
  for(const v of values){if(!Number.isFinite(v))return Infinity;const a=Math.abs(v);if(a>scale){sum=1+sum*(scale/a)**2;scale=a;}else if(a)sum+=(a/scale)**2;}
  return scale*Math.sqrt(sum);
}

function validateSparseMatrix(a){
  const {n,rowPtr,colIndex,values}=a;
  if(!Number.isInteger(n)||n<1||n>=2**31||rowPtr.length!==n+1||rowPtr[0]!==0
    ||rowPtr[n]!==values.length||values.length!==colIndex.length)throw new Error('Invalid sparse matrix dimensions.');
  for(let i=0;i<n;i++){
    if(!Number.isInteger(rowPtr[i+1])||rowPtr[i+1]<rowPtr[i])throw new Error('Invalid sparse row pointers.');
    let previous=-1;
    for(let k=rowPtr[i];k<rowPtr[i+1];k++){
      const j=colIndex[k];
      if(!Number.isInteger(j)||j<=previous||j>=n||!Number.isFinite(values[k]))throw new Error('Invalid sparse matrix entry.');
      previous=j;
    }
  }
}

function compressedColumns(a){
  validateSparseMatrix(a);
  const {n,rowPtr,colIndex,values}=a,ap=new Int32Array(n+1);
  // Exact structural zeros need no numerical factorization entry.
  for(let k=0;k<values.length;k++)if(values[k]!==0)ap[colIndex[k]+1]++;
  for(let j=0;j<n;j++)ap[j+1]+=ap[j];
  const ai=new Int32Array(ap[n]),ax=new Float64Array(ap[n]),next=ap.slice();
  for(let i=0;i<n;i++)for(let k=rowPtr[i];k<rowPtr[i+1];k++)if(values[k]!==0){
    const p=next[colIndex[k]]++;ai[p]=i;ax[p]=values[k];
  }
  return{ap,ai,ax};
}

export function solveSparseDirect(a,b,{tolerance=1e-10,maxRefinements=4,ordering='auto',preferredOrdering='amd',pivotTolerance=0.001,pivotFallback=true,
  rowPermutation,columnPermutation,btf}={}){
  if(b.length!==a.n||!b.every(Number.isFinite)||!Number.isFinite(tolerance)||tolerance<=0
    ||!Number.isInteger(maxRefinements)||maxRefinements<0||!['auto','amd','colamd','given'].includes(ordering)
    ||!['amd','colamd'].includes(preferredOrdering)||!Number.isFinite(pivotTolerance)||pivotTolerance<=0||pivotTolerance>1
    ||typeof pivotFallback!=='boolean'
    ||ordering!=='given'&&(rowPermutation!==undefined||columnPermutation!==undefined||btf!==undefined)
    ||ordering==='given'&&(maxRefinements>4||btf!==undefined&&typeof btf!=='boolean'))throw new Error('Invalid sparse solve controls or right-hand side.');
  const {ap,ai,ax}=compressedColumns(a),bnorm=norm(b);
  if(ordering==='given')return solveGiven(a,b,{ap,ai,ax,bnorm,rowPermutation,columnPermutation,
    tolerance,maxRefinements,pivotTolerance,btf:btf??false});
  if(!bnorm)return{x:new Float64Array(a.n),relativeResidual:0,refinements:0,factorNonzeros:0,backend:'klu-wasm'};
  const allocations=[];let factor=0;
  const copy=data=>{
    const p=wasm.malloc(data.byteLength);
    if(!p)throw new Error('WebAssembly sparse solver ran out of memory.');
    allocations.push(p);
    // malloc and factorization may grow memory: never retain a heap view.
    new data.constructor(wasm.memory.buffer,p,data.length).set(data);return p;
  };
  try{
    const apPointer=copy(ap),aiPointer=copy(ai),axPointer=copy(ax),rhsPointer=copy(Float64Array.from(b));
    const orderings=ordering==='auto'?[preferredOrdering,preferredOrdering==='amd'?'colamd':'amd']:[ordering],attempts=[];
    const configurations=orderings.map(ordering=>({ordering,pivotTolerance}));
    // Only after the normal factors miss the original-system accuracy gate,
    // use maximum-magnitude partial pivoting. Start with the last ordering
    // tried; COLAMD recovered the retained refined coupled failure.
    if(pivotFallback&&pivotTolerance<1)for(const ordering of [...orderings].reverse())configurations.push({ordering,pivotTolerance:1});
    for(let attempt=0;attempt<configurations.length;attempt++){
      const {ordering:selected,pivotTolerance:threshold}=configurations[attempt],last=attempt===configurations.length-1;
      // Ordering changes sparse elimination and its rounding error. Retrying
      // must factor the same matrix and pass the same original-system gate.
      // Free the rejected factors before allocating the alternative ordering.
      if(factor){wasm.mses_klu_free(factor);factor=0;}
      factor=wasm.mses_klu_factor_with_pivot(a.n,apPointer,aiPointer,axPointer,selected==='amd'?0:1,threshold);
      if(!factor){
        const status=wasm.mses_klu_status();
        attempts.push({ordering:selected,pivotTolerance:threshold,status});
        if(status===1&&!last)continue;
        const error=new Error(`Sparse LU factorization failed (KLU status ${status}).`);
        error.attempts=attempts;throw error;
      }
      const solve=rhs=>{
        new Float64Array(wasm.memory.buffer,rhsPointer,a.n).set(rhs);
        if(!wasm.mses_klu_solve(factor,rhsPointer))throw new Error(`Sparse LU solve failed (KLU status ${wasm.mses_klu_status()}).`);
        return new Float64Array(wasm.memory.buffer,rhsPointer,a.n).slice();
      };
      const x=solve(b);let relativeResidual=Infinity,refinements=0;
      for(;refinements<=maxRefinements;refinements++){
        const residual=compensatedSparseResidual(a,x,b);
        relativeResidual=norm(residual)/bnorm;
        if(x.every(Number.isFinite)&&relativeResidual<=tolerance)break;
        if(!Number.isFinite(relativeResidual)||refinements===maxRefinements)break;
        const correction=solve(residual);for(let i=0;i<a.n;i++)x[i]+=correction[i];
      }
      const factorNonzeros=wasm.mses_klu_nnz(factor);
      attempts.push({ordering:selected,pivotTolerance:threshold,relativeResidual,refinements,factorNonzeros});
      if(x.every(Number.isFinite)&&relativeResidual<=tolerance)
        return{x,relativeResidual,refinements,factorNonzeros,backend:'klu-wasm',ordering:selected,pivotTolerance:threshold,attempts,
          residualEvaluation:'compensated-original-system'};
      if(last){
        const error=new Error(`Sparse LU missed the requested relative residual (${relativeResidual} > ${tolerance}).`);
        error.code='KLU_RESIDUAL_LIMIT';error.relativeResidual=relativeResidual;error.attempts=attempts;
        error.diagnostics={relativeResidual,tolerance,residualEvaluation:'compensated-original-system',attempts:structuredClone(attempts)};throw error;
      }
    }
  }finally{
    if(factor)wasm.mses_klu_free(factor);
    for(const p of allocations)wasm.free(p);
  }
}

// Row alignment improves AMD's structural graph without permuting the solution.
// Retain the original-system accuracy gate and original ordering as a numerical
// fallback: a fill-reducing ordering must never relax the accepted tolerance.
export function solveSparseDirectAligned(a,b,rowPermutation,options={}) {
  const order=givenPermutation(rowPermutation,a.n,'rowPermutation');
  // Validate the original CSR before indexing into it.
  validateSparseMatrix(a);
  if(b.length!==a.n||!b.every(Number.isFinite))throw new Error('Invalid sparse right-hand side.');
  const rowPtr=new Int32Array(a.n+1),colIndex=new Int32Array(a.colIndex.length),values=new Float64Array(a.values.length);
  const rhs=new Float64Array(a.n);
  let next=0;
  for(let i=0;i<a.n;i++){
    const row=order[i];rhs[i]=b[row];
    for(let k=a.rowPtr[row];k<a.rowPtr[row+1];k++){colIndex[next]=a.colIndex[k];values[next++]=a.values[k];}
    rowPtr[i+1]=next;
  }
  let rejected;
  try {
    const result=solveSparseDirect({n:a.n,rowPtr,colIndex,values},rhs,options);
    const residual=compensatedSparseResidual(a,result.x,b),bnorm=norm(b);
    const relativeResidual=norm(residual)/(bnorm||1);
    if(relativeResidual<=(options.tolerance??1e-10))
      return {...result,relativeResidual,equationOrdering:'aligned'};
    rejected={code:'KLU_RESIDUAL_LIMIT',relativeResidual};
  } catch(error) {
    if(error.code!=='KLU_RESIDUAL_LIMIT'&&!error.attempts)throw error;
    rejected={code:error.code,attempts:error.attempts};
  }
  return {...solveSparseDirect(a,b,options),equationOrdering:'original',alignmentFallback:rejected};
}

function givenPermutation(input,n,name){
  if(!(Array.isArray(input)||ArrayBuffer.isView(input))||input.length!==n)
    throw new Error(`${name} must be a bijection of 0..n-1.`);
  const seen=new Uint8Array(n),result=new Int32Array(n);
  for(let k=0;k<n;k++){
    const i=input[k];
    if(!Number.isInteger(i)||i<0||i>=n||seen[i])throw new Error(`${name} must be a bijection of 0..n-1.`);
    seen[i]=1;result[k]=i;
  }
  return result;
}

// Opt-in symbolic A(P,Q), with position-to-original-index permutations.
// KLU maps the original RHS/solution. Numeric partial pivoting remains enabled;
// there is one factorization attempt, with no ordering or pivot fallback.
function solveGiven(a,b,{ap,ai,ax,bnorm,rowPermutation,columnPermutation,tolerance,maxRefinements,pivotTolerance,btf}){
  const p=givenPermutation(rowPermutation,a.n,'rowPermutation'),q=givenPermutation(columnPermutation,a.n,'columnPermutation');
  const allocations=[],attempts=[];let factor=0;
  const copy=data=>{
    const pointer=wasm.malloc(Math.max(1,data.byteLength));
    if(!pointer)throw new Error('WebAssembly sparse solver ran out of memory.');
    allocations.push(pointer);
    new data.constructor(wasm.memory.buffer,pointer,data.length).set(data);return pointer;
  };
  try{
    const apPointer=copy(ap),aiPointer=copy(ai),axPointer=copy(ax),pPointer=copy(p),qPointer=copy(q),rhsPointer=copy(Float64Array.from(b));
    factor=wasm.mses_klu_factor_given(a.n,apPointer,aiPointer,axPointer,pPointer,qPointer,pivotTolerance,btf?1:0);
    if(!factor){
      const status=wasm.mses_klu_status();attempts.push({ordering:'given',pivotTolerance,btf,status});
      const error=new Error(`Given-order sparse LU factorization failed (KLU status ${status}).`);
      error.code='KLU_GIVEN_FACTORIZATION';error.status=status;error.attempts=attempts;throw error;
    }
    const actualBtf=wasm.mses_klu_btf(factor)!==0,symbolicBlocks=wasm.mses_klu_blocks(factor);
    if(actualBtf!==btf)throw new Error('KLU did not retain the requested BTF control.');
    const solve=rhs=>{
      new Float64Array(wasm.memory.buffer,rhsPointer,a.n).set(rhs);
      if(!wasm.mses_klu_solve(factor,rhsPointer)){
        const error=new Error(`Given-order sparse LU solve failed (KLU status ${wasm.mses_klu_status()}).`);
        error.code='KLU_GIVEN_SOLVE';error.attempts=attempts;throw error;
      }
      return new Float64Array(wasm.memory.buffer,rhsPointer,a.n).slice();
    };
    const x=solve(b);let relativeResidual=Infinity,refinements=0;
    for(;refinements<=maxRefinements;refinements++){
      const residual=compensatedSparseResidual(a,x,b);
      const rnorm=norm(residual);relativeResidual=bnorm?rnorm/bnorm:rnorm;
      if(x.every(Number.isFinite)&&relativeResidual<=tolerance)break;
      if(!Number.isFinite(relativeResidual)||refinements===maxRefinements)break;
      const correction=solve(residual);for(let i=0;i<a.n;i++)x[i]+=correction[i];
    }
    const factorNonzeros=wasm.mses_klu_nnz(factor);
    attempts.push({ordering:'given',pivotTolerance,btf:actualBtf,relativeResidual,refinements,factorNonzeros,symbolicBlocks});
    if(!x.every(Number.isFinite)||relativeResidual>tolerance){
      const error=new Error(`Given-order sparse LU missed the requested relative residual (${relativeResidual} > ${tolerance}).`);
      error.code='KLU_RESIDUAL_LIMIT';error.relativeResidual=relativeResidual;error.attempts=attempts;
      error.diagnostics={relativeResidual,tolerance,residualEvaluation:'compensated-original-system',attempts:structuredClone(attempts)};throw error;
    }
    return{x,relativeResidual,refinements,factorNonzeros,backend:'klu-wasm',ordering:'given',pivotTolerance,btf:actualBtf,symbolicBlocks,attempts,
      residualEvaluation:'compensated-original-system'};
  }finally{
    if(factor)wasm.mses_klu_free(factor);
    for(const pointer of allocations)wasm.free(pointer);
  }
}

// Multiple right-hand sides share one aligned numeric factor. Each answer
// still passes the original-system residual gate, with the ordinary pivoted
// solver as the per-right-hand-side fallback.
export function solveSparseDirectAlignedMany(a, rightHandSides, rowPermutation) {
  validateSparseMatrix(a);
  const order = givenPermutation(rowPermutation, a.n, 'rowPermutation');
  if (!Array.isArray(rightHandSides) || rightHandSides.some(b => !b || b.length !== a.n || typeof b.every !== 'function' || !b.every(Number.isFinite)))
    throw new Error('Invalid sparse right-hand sides.');
  if (!rightHandSides.length) return [];
  if (rightHandSides.every(b => norm(b) === 0)) return rightHandSides.map(() => ({
    x: new Float64Array(a.n), relativeResidual: 0, refinements: 0, factorNonzeros: 0, backend: 'klu-wasm', equationOrdering: 'aligned' }));
  const rowPtr = new Int32Array(a.n + 1), colIndex = new Int32Array(a.colIndex.length), values = new Float64Array(a.values.length);
  let next = 0;
  for (let i = 0; i < a.n; i++) {
    for (let k = a.rowPtr[order[i]]; k < a.rowPtr[order[i] + 1]; k++) {
      colIndex[next] = a.colIndex[k]; values[next++] = a.values[k];
    }
    rowPtr[i + 1] = next;
  }
  const { ap, ai, ax } = compressedColumns({ n: a.n, rowPtr, colIndex, values });
  const allocations = []; let factor = 0;
  const copy = data => {
    const pointer = wasm.malloc(Math.max(1, data.byteLength));
    if (!pointer) throw new Error('WebAssembly sparse solver ran out of memory.');
    allocations.push(pointer);
    new data.constructor(wasm.memory.buffer, pointer, data.length).set(data); return pointer;
  };
  try {
    const apPointer = copy(ap), aiPointer = copy(ai), axPointer = copy(ax), rhsPointer = copy(new Float64Array(a.n));
    factor = wasm.mses_klu_factor_with_pivot(a.n, apPointer, aiPointer, axPointer, 0, .001);
    const solve = b => {
      const heap = new Float64Array(wasm.memory.buffer, rhsPointer, a.n);
      for (let i = 0; i < a.n; i++) heap[i] = b[order[i]];
      if (!wasm.mses_klu_solve(factor, rhsPointer)) throw new Error('Batched sparse solve failed.');
      return new Float64Array(wasm.memory.buffer, rhsPointer, a.n).slice();
    };
    return rightHandSides.map(b => {
      const bnorm = norm(b);
      if (!bnorm) return { x: new Float64Array(a.n), relativeResidual: 0, refinements: 0, factorNonzeros: 0, backend: 'klu-wasm', equationOrdering: 'aligned' };
      if (factor) {
        const x = solve(b);
        for (let refinements = 0; refinements <= 4; refinements++) {
          const residual = compensatedSparseResidual(a, x, b), relativeResidual = norm(residual) / bnorm;
          if (x.every(Number.isFinite) && relativeResidual <= 1e-10)
            return { x, relativeResidual, refinements, factorNonzeros: wasm.mses_klu_nnz(factor), backend: 'klu-wasm', equationOrdering: 'aligned' };
          if (!Number.isFinite(relativeResidual) || refinements === 4) break;
          const dx = solve(residual); for (let i = 0; i < a.n; i++) x[i] += dx[i];
        }
      }
      return solveSparseDirectAligned(a, b, order);
    });
  } finally {
    if (factor) wasm.mses_klu_free(factor);
    for (const pointer of allocations) wasm.free(pointer);
  }
}
