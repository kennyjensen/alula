// SPDX-License-Identifier: GPL-2.0-or-later
// Compressed-row storage and a restarted, right-preconditioned GMRES backend.
// Only the preconditioner is approximate; success tests the original system.
import { factorLinear } from './linear.js';
export function sparseMatrix(columns){
  const n=columns.length,rowPtr=new Int32Array(n+1);
  const sorted=columns.map((row,i)=>[...new Set([...row,i])].sort((a,b)=>a-b));
  for(let i=0;i<n;i++)rowPtr[i+1]=rowPtr[i]+sorted[i].length;
  const colIndex=Int32Array.from(sorted.flat()),values=new Float64Array(colIndex.length);
  if(!n||colIndex.some(j=>j<0||j>=n))throw new Error('Invalid sparse matrix pattern.');
  return{n,rowPtr,colIndex,values};
}
export function sparseIndex(a,row,col){
  let lo=a.rowPtr[row],hi=a.rowPtr[row+1]-1;
  while(lo<=hi){const mid=(lo+hi)>>1,j=a.colIndex[mid];if(j===col)return mid;if(j<col)lo=mid+1;else hi=mid-1;}
  return -1;
}
export function sparseAdd(a,row,col,value){
  const index=sparseIndex(a,row,col);
  if(index<0)throw new Error(`Sparse matrix entry (${row},${col}) is outside the stencil.`);
  a.values[index]+=value;
}
export function sparseProduct(a,x){
  const y=new Float64Array(a.n);
  for(let i=0;i<a.n;i++)for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)y[i]+=a.values[p]*x[a.colIndex[p]];
  return y;
}
export function sparseDense(a){
  const result=new Float64Array(a.n*a.n);
  for(let i=0;i<a.n;i++)for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)result[i*a.n+a.colIndex[p]]=a.values[p];
  return result;
}

export function ilu0(a){
  const values=a.values.slice(),diagonal=new Int32Array(a.n),scales=new Float64Array(a.n);
  for(let i=0;i<a.n;i++){
    diagonal[i]=sparseIndex(a,i,i);
    for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)scales[i]=Math.max(scales[i],Math.abs(values[p]));
    if(!(scales[i]>0)||!Number.isFinite(scales[i]))throw new Error('Invalid sparse preconditioner row.');
    for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)values[p]/=scales[i];
  }
  let shiftedPivots=0;
  for(let i=0;i<a.n;i++){
    for(let p=a.rowPtr[i];p<diagonal[i];p++){
      const j=a.colIndex[p],factor=values[p]/values[diagonal[j]];values[p]=factor;
      for(let q=diagonal[j]+1;q<a.rowPtr[j+1];q++){
        const index=sparseIndex(a,i,a.colIndex[q]);if(index>=0)values[index]-=factor*values[q];
      }
    }
    if(Math.abs(values[diagonal[i]])<1e-10){values[diagonal[i]]=Math.sign(values[diagonal[i]]||1)*1e-10;shiftedPivots++;}
  }
  if(!values.every(Number.isFinite))throw new Error('Nonfinite sparse preconditioner.');
  const solve=rhs=>{
    const x=Float64Array.from(rhs,(v,i)=>v/scales[i]);
    for(let i=0;i<a.n;i++)for(let p=a.rowPtr[i];p<diagonal[i];p++)x[i]-=values[p]*x[a.colIndex[p]];
    for(let i=a.n-1;i>=0;i--){
      for(let p=diagonal[i]+1;p<a.rowPtr[i+1];p++)x[i]-=values[p]*x[a.colIndex[p]];
      x[i]/=values[diagonal[i]];
    }
    return x;
  };
  return{solve,shiftedPivots};
}

// Scalar level-of-fill ILU for the mixed potential / BL unknowns. Dropped
// entries and shifted pivots affect only preconditioning, never the matrix
// whose residual certifies the GMRES solve.
export function scalarIlu(a,{fillLevel=1}={}){
  if(fillLevel===0)return ilu0(a);
  if(!Number.isInteger(fillLevel)||fillLevel<0||fillLevel>3)throw new Error('Invalid scalar ILU fill level.');
  const scales=new Float64Array(a.n),diagonal=new Float64Array(a.n),rows=[],lookup=[];let shiftedPivots=0;
  for(let i=0;i<a.n;i++){
    let scale=0;for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)scale=Math.max(scale,Math.abs(a.values[p]));
    if(!(scale>0)||!Number.isFinite(scale))throw new Error('Invalid sparse preconditioner row.');scales[i]=scale;
    const row=[],map=new Map();
    for(let p=a.rowPtr[i];p<a.rowPtr[i+1];p++)if(a.values[p]!==0||a.colIndex[p]===i){
      const entry={col:a.colIndex[p],value:a.values[p]/scale,level:0};row.push(entry);map.set(entry.col,entry);
    }
    rows.push(row);lookup.push(map);
  }
  for(let i=0;i<a.n;i++){
    const row=rows[i],map=lookup[i];
    for(let index=0;index<row.length;index++){
      const entry=row[index],j=entry.col;if(j>=i)break;
      const factor=entry.value/diagonal[j];entry.value=factor;
      for(const upper of rows[j])if(upper.col>j){
        const level=entry.level+upper.level+1;let target=map.get(upper.col);
        if(!target&&level<=fillLevel){
          target={col:upper.col,value:0,level};map.set(upper.col,target);
          let lo=index+1,hi=row.length;while(lo<hi){const mid=(lo+hi)>>1;if(row[mid].col<upper.col)lo=mid+1;else hi=mid;}row.splice(lo,0,target);
        }
        if(target){target.value-=factor*upper.value;target.level=Math.min(target.level,level);}
      }
    }
    const pivot=map.get(i);
    if(Math.abs(pivot.value)<1e-10){pivot.value=Math.sign(pivot.value||1)*1e-10;shiftedPivots++;}
    diagonal[i]=pivot.value;
    if(row.some(e=>!Number.isFinite(e.value)))throw new Error('Nonfinite scalar ILU preconditioner.');
  }
  const solve=rhs=>{
    const x=Float64Array.from(rhs,(v,i)=>v/scales[i]);
    for(let i=0;i<a.n;i++)for(const entry of rows[i]){if(entry.col>=i)break;x[i]-=entry.value*x[entry.col];}
    for(let i=a.n-1;i>=0;i--){for(const entry of rows[i])if(entry.col>i)x[i]-=entry.value*x[entry.col];x[i]/=diagonal[i];}
    return x;
  };
  return{solve,shiftedPivots};
}

const dot=(a,b)=>a.reduce((sum,v,i)=>sum+v*b[i],0);
const norm=a=>Math.sqrt(dot(a,a));

// Keep the four Euler equations together during incomplete elimination.
// Scalar pivots can be tiny near a wall even when the local acoustic/advective
// block is well conditioned. Invert each diagonal block with pivoted LU.
export function blockIlu(a,{fillLevel=1}={}){
  if(a.n%4!==0)throw new Error('Four-variable ILU needs a multiple of four rows.');
  if(!Number.isInteger(fillLevel)||fillLevel<0||fillLevel>3)throw new Error('Invalid block ILU fill level.');
  const count=a.n/4,rows=Array.from({length:count},(_,i)=>{
    const columns=new Set();
    for(let k=0;k<4;k++)for(let p=a.rowPtr[4*i+k];p<a.rowPtr[4*i+k+1];p++)columns.add(Math.floor(a.colIndex[p]/4));
    return[...columns].sort((a,b)=>a-b).map(j=>{
      const data=new Float64Array(16);
      for(let row=0;row<4;row++)for(let col=0;col<4;col++){const index=sparseIndex(a,4*i+row,4*j+col);if(index>=0)data[4*row+col]=a.values[index];}
      return{col:j,data,level:0};
    });
  });
  const lookup=rows.map(row=>new Map(row.map(b=>[b.col,b]))),inverses=[];
  const product=(a,b)=>{
    const c=new Float64Array(16);
    for(let i=0;i<4;i++)for(let k=0;k<4;k++)for(let j=0;j<4;j++)c[4*i+j]+=a[4*i+k]*b[4*k+j];
    return c;
  };
  let shiftedPivots=0;
  for(let i=0;i<count;i++){
    for(const block of rows[i]){
      const j=block.col;if(j>=i)break;
      const lower=product(block.data,inverses[j]);block.data.set(lower);
      for(const upper of rows[j])if(upper.col>j){
        const level=block.level+upper.level+1;
        let target=lookup[i].get(upper.col);
        if(!target&&level<=fillLevel){
          target={col:upper.col,data:new Float64Array(16),level};lookup[i].set(upper.col,target);rows[i].push(target);
        }
        if(!target)continue;
        target.level=Math.min(target.level,level);
        const correction=product(lower,upper.data);for(let k=0;k<16;k++)target.data[k]-=correction[k];
      }
      // New fill columns are greater than the current elimination column;
      // sorted insertion preserves the current iterator's processed prefix.
      rows[i].sort((a,b)=>a.col-b.col);
    }
    const diagonal=lookup[i].get(i).data;let solve;
    try{solve=factorLinear(diagonal,4);}catch{
      const scale=Math.max(...diagonal.map(Math.abs));
      for(let k=0;k<4;k++)diagonal[5*k]+=Math.max(scale,1)*1e-8;
      solve=factorLinear(diagonal,4);shiftedPivots++;
    }
    const inverse=new Float64Array(16);
    for(let col=0;col<4;col++){const rhs=new Float64Array(4);rhs[col]=1;const x=solve(rhs);for(let row=0;row<4;row++)inverse[4*row+col]=x[row];}
    inverses.push(inverse);
  }
  const solve=rhs=>{
    const x=Float64Array.from(rhs);
    const subtract=(i,b)=>{for(let row=0;row<4;row++)for(let col=0;col<4;col++)x[4*i+row]-=b.data[4*row+col]*x[4*b.col+col];};
    for(let i=0;i<count;i++)for(const b of rows[i]){if(b.col>=i)break;subtract(i,b);}
    for(let i=count-1;i>=0;i--){
      for(const b of rows[i])if(b.col>i)subtract(i,b);
      const value=x.slice(4*i,4*i+4),inverse=inverses[i];
      for(let row=0;row<4;row++){x[4*i+row]=0;for(let col=0;col<4;col++)x[4*i+row]+=inverse[4*row+col]*value[col];}
    }
    return x;
  };
  return{solve,shiftedPivots};
}

export function solveSparse(a,rhs,{tolerance=1e-10,restart=50,maxIterations=1000,blockSize=1,fillLevel=blockSize===4?1:0,preconditionerMatrix=a,initial}={}){
  if(rhs.length!==a.n||!Array.from(rhs).every(Number.isFinite)||!a.values.every(Number.isFinite))throw new Error('Invalid sparse linear system.');
  if(!(tolerance>0)||!Number.isInteger(restart)||restart<1||!Number.isInteger(maxIterations)||maxIterations<1)throw new Error('Invalid GMRES controls.');
  if(initial&&(initial.length!==a.n||!Array.from(initial).every(Number.isFinite)))throw new Error('Invalid GMRES initial state.');
  const x=initial?Float64Array.from(initial):new Float64Array(a.n),rhsNorm=norm(rhs);
  if(rhsNorm===0)return{x:new Float64Array(a.n),iterations:0,relativeResidual:0,shiftedPivots:0};
  if(![1,4].includes(blockSize))throw new Error('Unsupported ILU block size.');
  if(preconditionerMatrix.n!==a.n)throw new Error('Sparse preconditioner dimensions differ.');
  const preconditioner=blockSize===4?blockIlu(preconditionerMatrix,{fillLevel}):scalarIlu(preconditionerMatrix,{fillLevel}),precondition=preconditioner.solve;
  let iterations=0,relativeResidual=Infinity;
  while(iterations<maxIterations){
    const ax=sparseProduct(a,x),r=Float64Array.from(rhs,(v,i)=>v-ax[i]);
    relativeResidual=norm(r)/rhsNorm;
    if(relativeResidual<=tolerance)return{x,iterations,relativeResidual,shiftedPivots:preconditioner.shiftedPivots};
    const beta=norm(r),width=Math.min(restart,maxIterations-iterations);
    if(!(beta>0)||!Number.isFinite(beta))throw new Error('GMRES preconditioner breakdown.');
    const basis=[r.map(v=>v/beta)],directions=[],h=new Float64Array((width+1)*width),g=new Float64Array(width+1),c=new Float64Array(width),s=new Float64Array(width);g[0]=beta;
    let used=0;
    for(let j=0;j<width;j++){
      const direction=precondition(basis[j]),w=sparseProduct(a,direction);directions.push(direction);
      // Reorthogonalize: nonnormal Euler operators otherwise lose Arnoldi
      // orthogonality before a sufficiently accurate Newton step is available.
      for(let pass=0;pass<2;pass++)for(let i=0;i<=j;i++){
        const value=dot(w,basis[i]);h[i*width+j]+=value;
        for(let k=0;k<a.n;k++)w[k]-=value*basis[i][k];
      }
      const length=norm(w);h[(j+1)*width+j]=length;
      if(length>0)basis.push(w.map(v=>v/length));
      for(let i=0;i<j;i++){
        const u=h[i*width+j],v=h[(i+1)*width+j];h[i*width+j]=c[i]*u+s[i]*v;h[(i+1)*width+j]=-s[i]*u+c[i]*v;
      }
      const u=h[j*width+j],v=h[(j+1)*width+j],magnitude=Math.hypot(u,v);
      if(!(magnitude>0))throw new Error('GMRES Arnoldi breakdown.');
      c[j]=u/magnitude;s[j]=v/magnitude;h[j*width+j]=magnitude;h[(j+1)*width+j]=0;
      g[j+1]=-s[j]*g[j];g[j]*=c[j];used=j+1;iterations++;
      if(length===0||Math.abs(g[j+1])<=tolerance*rhsNorm)break;
    }
    const y=g.slice(0,used);
    for(let i=used-1;i>=0;i--){for(let j=i+1;j<used;j++)y[i]-=h[i*width+j]*y[j];y[i]/=h[i*width+i];}
    for(let i=0;i<used;i++)for(let k=0;k<a.n;k++)x[k]+=y[i]*directions[i][k];
    if(!x.every(Number.isFinite))throw new Error('Nonfinite GMRES solution.');
  }
  const ax=sparseProduct(a,x);relativeResidual=norm(Float64Array.from(rhs,(v,i)=>v-ax[i]))/rhsNorm;
  if(relativeResidual>tolerance){
    const error=new Error(`GMRES iteration limit: relative residual ${relativeResidual}.`);
    Object.assign(error,{code:'GMRES_ITERATION_LIMIT',solution:x,relativeResidual,iterations});throw error;
  }
  return{x,iterations,relativeResidual,shiftedPivots:preconditioner.shiftedPivots};
}
