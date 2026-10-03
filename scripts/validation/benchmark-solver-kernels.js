// SPDX-License-Identifier: GPL-2.0-or-later
// Compare an unchanged source snapshot with the current kernels and full NLR
// Jacobian. node scripts/validation/benchmark-solver-kernels.js BASELINE_ROOT
import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { requirePositiveSimplePolygon as polygon } from '../../src/geometry/simple-polygon.js';
import { linearizeStreamtubeCell as cell } from '../../src/euler/streamtube-linearization.js';
import { createCoupledStreamtubeBody as body } from '../../src/euler/streamtube-coupled.js';
import { solveSparseDirectAligned as solve } from '../../src/numerics/klu.js';
import { createStreamtubeStationOrdering } from '../../src/euler/streamtube-station-ordering.js';
const baseline = process.argv[2];
if (!baseline) throw new Error('Supply the baseline source root.');
const load = name => import(pathToFileURL(path.resolve(baseline, 'src', name)).href);
const oldPolygon = (await load('geometry/simple-polygon.js')).requirePositiveSimplePolygon;
const oldCell = (await load('euler/streamtube-linearization.js')).linearizeStreamtubeCell;
const oldBody = (await load('euler/streamtube-coupled.js')).createCoupledStreamtubeBody;
const oldSolve = (await load('numerics/klu.js')).solveSparseDirectAligned;
const parameters = { lower: [{x:0,y:0},{x:1,y:.01},{x:2,y:.03}], upper: [{x:0,y:1},{x:1,y:1.1},{x:2,y:1.15}],
  densities:[1,1.01], massFlow:.7, stagnationEnthalpy:4, gamma:1.4, geometryDomain:'positive-simple', transportSpeeds:[.65,.66] };
const linearizers = [oldCell(parameters),cell(parameters)];
const tangents = Array.from({length:256}, (_,i) => ({ lower: Array.from({length:3},(_,j)=>({x:Math.sin(i+j)*.01,y:Math.cos(i-j)*.01})),
  upper: Array.from({length:3},(_,j)=>({x:Math.sin(i-j)*.02,y:Math.cos(i+j)*.02})), densities:[i*.0001,-i*.0002],
  massFlow:.003, stagnationEnthalpy:.01, transportSpeeds:[.01,-.02] }));
for(const tangent of tangents) assert.deepEqual(linearizers[1].apply(tangent),linearizers[0].apply(tangent));
const retained = linearizers[1].apply(tangents[0]), retainedCopy = structuredClone(retained);
linearizers[1].apply(tangents[1]); assert.deepEqual(retained,retainedCopy);
const polygons=Array.from({length:256},(_,i)=>Array.from({length:i%2?4:6},(_,j)=>{
  const a=2*Math.PI*j/(i%2?4:6),r=1+.2*Math.sin(i+j);return{x:r*Math.cos(a),y:r*Math.sin(a)};
}));
for(const p of polygons)assert.equal(polygon(p),oldPolygon(p));
function timed(fn,count){let sum=0;const start=performance.now();for(let i=0;i<count;i++)sum+=fn(i);return {milliseconds:performance.now()-start,checksum:sum};}
const polygonRuns=[],derivativeRuns=[];
for(const which of [0,1,1,0]){
 const fn=[oldPolygon,polygon][which];polygonRuns.push({method:which?'optimized':'before',...timed(i=>fn(polygons[i%256]),400000)});
 derivativeRuns.push({method:which?'optimized':'before',...timed(i=>linearizers[which].apply(tangents[i%256]).streamwiseResidual,100000)});
}
const fixture=new URL('../../docs/solver-reliability/nlr64x24-subsonic/performance/late-coupled-checkpoint.json.gz',import.meta.url);
const r=JSON.parse(gunzipSync(fs.readFileSync(fixture))).restart;
const systems=[oldBody,body].map(make=>make(r.input,{...r.options,initialEuler:r.initialEuler,initialBL:r.initialBL}));
const jacobianRuns=[];let original,matrix;
for(const which of [0,1,1,0]){
 const start=performance.now(),s=systems[which],a=s.jacobian(s.initial),milliseconds=performance.now()-start;
 original??=a; assert.deepEqual(a,original);matrix=a;
 jacobianRuns.push({method:which?'optimized':'before',milliseconds});
}
const s=systems[1],rhs=s.evaluate(s.initial).residual.map(v=>-v);
const matching=createStreamtubeStationOrdering({matrix,layout:s.euler.layout,stations:s.bl.stations});
const rows=new Int32Array(matrix.n);for(let k=0;k<matrix.n;k++)rows[matching.Quser[k]]=matching.Puser[k];
const sparseRuns=[];let direction;
for(const which of [0,1,1,0]){
 const start=performance.now(),v=[oldSolve,solve][which](matrix,rhs,rows),milliseconds=performance.now()-start;
 direction??=v.x;assert.deepEqual(v.x,direction);assert(v.relativeResidual<=1e-10);
 sparseRuns.push({method:which?'optimized':'before',milliseconds,relativeResidual:v.relativeResidual});
}
console.log(JSON.stringify({baseline,polygonRuns,derivativeRuns,jacobianRuns,sparseRuns,
 unknowns:matrix.n,nonzeros:matrix.values.length,bitIdenticalJacobian:true,bitIdenticalLinearDirection:true},null,2));
