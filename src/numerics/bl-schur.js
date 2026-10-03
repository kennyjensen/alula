// SPDX-License-Identifier: GPL-2.0-or-later
import { factorLinear, solveLinear } from './linear.js';

// Exact elimination of the three BL variables at every station. The remaining
// system contains every edge velocity and every cross-element interaction.
// This is a linear factorization of one simultaneous Newton step, not a
// prescribed-velocity nonlinear BL march. Separated-flow local blocks can be
// singular even when the coupled matrix is invertible: retain a full fallback.
export function solveBlSchur(matrix,rhs,order){
  const size=rhs.length,n=size/4,width=n+1;
  if(!Number.isInteger(n)||matrix.length!==size*size||order.length!==n||new Set(order).size!==n
    ||order.some(i=>!Number.isInteger(i)||i<0||i>=n))throw new Error('Invalid BL Schur dimensions or station order.');
  const backwardError=x=>{
    let error=0;
    for(let i=0;i<size;i++){
      let value=-rhs[i],scale=Math.abs(rhs[i]);
      for(let j=0;j<size;j++){const term=matrix[i*size+j]*x[j];value+=term;scale+=Math.abs(term);}
      error=Math.max(error,Math.abs(value)/Math.max(scale,Number.MIN_VALUE));
    }
    return error;
  };
  let fallbackReason;
  try{
    // z stores [B^-1 C | B^-1 rhs] in station-major BL-variable order.
    const z=new Float64Array(3*n*width),done=new Uint8Array(n);
    for(const id of order){
      const local=new Float64Array(9),work=new Float64Array(3*width);
      for(let row=0;row<3;row++){
        const offset=(4*id+row)*size;
        for(let col=0;col<3;col++)local[3*row+col]=matrix[offset+4*id+col];
        for(let j=0;j<n;j++)work[row*width+j]=matrix[offset+4*j+3];
        for(let j=0;j<n;j++){
          if(j===id)continue;
          for(let k=0;k<3;k++){
            const value=matrix[offset+4*j+k];if(value===0)continue;
            if(!done[j])throw new Error('BL block has a dependency outside its elimination order.');
            const source=(3*j+k)*width;
            for(let col=0;col<width;col++)work[row*width+col]-=value*z[source+col];
          }
        }
        work[row*width+n]+=rhs[4*id+row];
      }
      const solve=factorLinear(local,3);
      for(let col=0;col<width;col++){
        const solution=solve([work[col],work[width+col],work[2*width+col]]);
        for(let row=0;row<3;row++)z[(3*id+row)*width+col]=solution[row];
      }
      done[id]=1;
    }
    const schur=new Float64Array(n*n),b=new Float64Array(n);
    for(let i=0;i<n;i++){
      const offset=(4*i+3)*size;b[i]=rhs[4*i+3];
      for(let j=0;j<n;j++)schur[i*n+j]=matrix[offset+4*j+3];
      for(let j=0;j<n;j++)for(let k=0;k<3;k++){
        const value=matrix[offset+4*j+k];if(value===0)continue;
        const source=(3*j+k)*width;
        b[i]-=value*z[source+n];
        for(let col=0;col<n;col++)schur[i*n+col]-=value*z[source+col];
      }
    }
    const dq=solveLinear(schur,b),x=new Float64Array(size);
    for(let i=0;i<n;i++){
      x[4*i+3]=dq[i];
      for(let k=0;k<3;k++){
        const offset=(3*i+k)*width;let value=z[offset+n];
        for(let j=0;j<n;j++)value-=z[offset+j]*dq[j];
        x[4*i+k]=value;
      }
    }
    const error=backwardError(x);
    if(!x.every(Number.isFinite)||!(error<1e-9))throw new Error(`Schur backward error ${error}.`);
    return{x,method:'schur',backwardError:error};
  }catch(error){fallbackReason=error.message;}
  const solve=factorLinear(matrix,size),x=solve(rhs);
  let error=backwardError(x);
  // Separated intervals can make the Schur reduction ill conditioned. Reuse
  // the full factorization for iterative refinement before rejecting a step;
  // the acceptance threshold is unchanged and uses the original matrix.
  for(let iteration=0;error>=1e-9&&iteration<3;iteration++){
    const residual=rhs.slice();
    for(let i=0;i<size;i++){
      let sum=0,correction=0;
      for(let j=0;j<size;j++){
        const term=matrix[i*size+j]*x[j]-correction,next=sum+term;
        correction=(next-sum)-term;sum=next;
      }
      residual[i]-=sum;
    }
    const delta=solve(residual);
    for(let i=0;i<size;i++)x[i]+=delta[i];
    error=backwardError(x);
  }
  if(!(error<1e-9))throw new Error(`Full Newton linear solve backward error ${error}.`);
  return{x,method:'full',backwardError:error,fallbackReason};
}
