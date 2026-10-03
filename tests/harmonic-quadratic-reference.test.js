import test from 'node:test';
import assert from 'node:assert/strict';
import { solveHarmonicScalarReference } from '../src/geometry/tests/harmonic-scalar-reference.js';

const grid = n => Array.from({length:n+1},(_,i)=>Array.from({length:n+1},(_,j)=>{
  const x=i/n,y=j/n,d=.06*Math.sin(Math.PI*x)*Math.sin(Math.PI*y);
  return {x:x+d,y:y-.4*d};
}));

test('assembled Q2 reference recovers a physical harmonic quadratic on distorted cells', () => {
  const nodes=grid(5),exact=p=>2*p.x*p.y+3*p.x-.4*p.y+7;
  const labels=nodes.map((row,i)=>row.map((p,j)=>!i||i===5||!j||j===5?exact(p):Math.sin(i+j)));
  const original=structuredClone({nodes,labels});
  for(const refinement of [1,2,4]){
    const r=solveHarmonicScalarReference({nodes,labels},{degree:2,refinement});
    r.values.forEach((row,i)=>row.forEach((v,j)=>assert.ok(Math.abs(v-exact(nodes[i][j]))<3e-12)));
    assert.equal(r.unknowns,(10*refinement-1)**2);
  }
  assert.deepEqual({nodes,labels},original);
});

test('Q2 and Q1 refine toward an independent Fourier solution with fixed polygon boundary data', t => {
  // Unit square, zero on three edges, triangular tent on the top edge.
  // Sine coefficients follow from integrating the two linear top segments.
  // The tent's peak is an original boundary node, so subdivision changes
  // neither geometry nor Dirichlet data. Evaluate the series only inside.
  const nodes=grid(6),labels=nodes.map((row,i)=>row.map((p,j)=>j===6?2*Math.min(p.x,1-p.x):0));
  const exact=(p,terms)=>{
    let value=0;
    for(let n=1;n<=terms;n+=2){const k=n*Math.PI;
      value+=8*Math.sin(k/2)/(k*k)*Math.sin(k*p.x)*Math.exp(k*(p.y-1))*(-Math.expm1(-2*k*p.y))/(-Math.expm1(-2*k));
    }
    return value;
  };
  const errors={};
  for(const degree of [1,2]){
    errors[degree]=[1,2,4].map(refinement=>{
      const r=solveHarmonicScalarReference({nodes,labels},{degree,refinement});let error=0;
      for(let i=1;i<6;i++)for(let j=1;j<6;j++){
        assert.ok(Math.abs(exact(nodes[i][j],64)-exact(nodes[i][j],128))<1e-14);
        error=Math.max(error,Math.abs(r.values[i][j]-exact(nodes[i][j],128)));
      }
      return error;
    });
    assert.ok(errors[degree][1]<.35*errors[degree][0]);
    assert.ok(errors[degree][2]<.35*errors[degree][1]);
  }
  assert.ok(errors[2][2]<1e-5);
  assert.ok(errors[2][1]<errors[1][2]);
  t.diagnostic(JSON.stringify(errors));
});

test('Q2 preserves constant fields and rejects unsupported degree, quadrature and budget', () => {
  const nodes=grid(3),labels=nodes.map(row=>row.map(()=>12));
  assert.equal(solveHarmonicScalarReference({nodes,labels},{degree:2}).maximumError,0);
  for(const options of [{degree:0},{degree:3},{degree:2,quadratureOrder:2},{degree:1,quadratureOrder:5},{degree:2,maxUnknowns:24}])
    assert.throws(()=>solveHarmonicScalarReference({nodes,labels},options));
});
