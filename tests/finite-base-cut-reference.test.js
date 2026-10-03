import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { streamtubeCutEdge, streamtubeCutReferencePoints } from '../src/euler/streamtube-cut-reference.js';
import { matchPassageDensityGuides } from '../src/euler/passage-density-guides.js';
import { createPanelPotentialBlocks } from '../src/geometry/panel-potential-blocks.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { createContourArc } from '../src/geometry/contour-arc.js';
import { naca4 } from '../src/geometry/airfoil.js';
const close=(a,b)=>assert.ok(Math.abs(a-b)<2e-13,`scalar mismatch ${a-b}`);

test('captured NLR wake request exposes lower-corner reference and midpoint reference retains exact solid banks',()=>{
  const r=JSON.parse(fs.readFileSync(new URL('../docs/nlr-finite-base/finer-potential-path/guide-prefix.json',import.meta.url)));
  const b=r.input.elements[1],curve=createSurfaceContourCurve(b.points,b),profile={curve,stag:curve.length/2};
  const center=streamtubeCutEdge(profile,b,'wake'),request=r.error.guideRequest;
  assert.ok(request.x<center.x);assert.equal(request.referencePoints[0].x,curve.evaluate(curve.length).point.x);
  const guides={lower:structuredClone(request.referencePoints),upper:structuredClone(request.referencePoints)};
  guides.upper[0]={...guides.upper[0],...curve.evaluate(0).point};
  const original=JSON.stringify(guides),points=streamtubeCutReferencePoints(b,guides,'wake',0,guides.lower.length-1);
  assert.equal(points[0].x,center.x);assert.equal(points[0].y,center.y);
  assert.equal(points[0].potential,guides.lower[0].potential);
  assert.ok(points.slice(1).every(p=>p.x>center.x));
  assert.equal(JSON.stringify(guides),original);
});

test('sharp and incoming reference arithmetic remains exactly unchanged',()=>{
  const lower=[{x:1,y:2,potential:3},{x:4,y:5,potential:6}],upper=[{x:1.1,y:2.2},{x:4,y:5}],guides={lower,upper};
  for(const [body,end] of [[{},'wake'],[{trailingEdge:{kind:'finite-base'}},'upstream']]){
    const p=streamtubeCutReferencePoints(body,guides,end,0,1);assert.deepEqual(p,lower);assert.equal(p[0],lower[0]);
  }
  const point={x:1,y:2},profile={curve:{branch:()=>({point})},stag:0};
  assert.equal(streamtubeCutEdge(profile,{},'wake'),point);
  assert.equal(streamtubeCutEdge(profile,{trailingEdge:{kind:'finite-base'}},'upstream'),point);
});

test('passage spacing for an oblique finite base stays downstream of the shared center and preserves both TE corners',()=>{
  const source=naca4('0012',40),last=source.length-1,angle=-.5,c=Math.cos(angle),sn=Math.sin(angle);
  const transform=p=>({x:c*p.x-sn*p.y,y:sn*p.x+c*p.y});
  const surface=source.map((p,i)=>transform({x:p.x,y:p.y+(i<last/2?1:-1)*.05*p.x}));
  const body={points:[...surface,{...surface[0]}],trailingEdge:{kind:'finite-base',upperIndex:0,lowerIndex:last}};
  const curve=createSurfaceContourCurve(body.points,body),arc=createContourArc(curve),stag=curve.length/2,origin=arc.at(stag);
  const phase=s=>1+Math.abs(s-stag)/stag+(s<stag?.2:0)*Math.abs(s-stag)/stag;
  const profile={curve,stag,phase,upstream:'upstream',wake:'wake'},center=streamtubeCutEdge(profile,body,'wake');
  const blocks=createPanelPotentialBlocks({bodies:[{leading:1,trailing:{upper:2.2,lower:2},inlet:0,outletIncrement:1}],outer:[{inlet:0,outlet:3.5},{inlet:0,outlet:3.5}]});
  const N=32,ranks=Array.from({length:3*N+1},(_,i)=>i/N),le=curve.branch('upper',0,stag).point;
  const input={bodies:[{...body,leadingIndex:N,trailingIndex:2*N,surfaceFractions:{upper:[],lower:[]}}]};
  const guides=[{upper:[],lower:[]}],descriptors={};
  for(const side of ['upper','lower']){
    const at=f=>(side==='upper'?-1:1)*(arc.at(curve.branch(side,f,stag).parameter)-origin);
    descriptors[side]={at,rows:[{rank:1,value:0},{rank:2,value:at(1)}],fractionAtPosition:s=>{
      let lo=0,hi=1;for(let k=0;k<56;k++){const mid=.5*(lo+hi);if(at(mid)<s)lo=mid;else hi=mid;}return .5*(lo+hi);}};
    for(let i=0;i<ranks.length;i++){
      if(i<N)guides[0][side].push({x:le.x-1+i/N,y:le.y,potential:i/N});
      else if(i>2*N)guides[0][side].push({x:center.x+(i-2*N)/N,y:center.y,potential:(side==='upper'?2.2:2)+(i-2*N)/N});
      else{const f=(i-N)/N,v=curve.branch(side,f,stag);input.bodies[0].surfaceFractions[side].push(f);guides[0][side].push({...v.point,potential:phase(v.parameter)});}
    }
  }
  const outer=[-2,2].map(y=>ranks.map(r=>({x:r-1,y,potential:r}))),before=structuredClone({guides,outer}),requests=[];
  const demand=Array.from({length:17},(_,i)=>.5*(1-Math.cos(Math.PI*i/16)));
  const fitted=matchPassageDensityGuides({input,profiles:[profile],guides,outer,blocks,ranks,surfaceMaps:[descriptors],resolvedDemands:[{upper:demand,lower:demand}],
    sampleCutX:(path,x)=>{if(path==='wake'){requests.push(x);assert.ok(x>center.x,'wake spacing must not request a point behind its center origin');}
      return{x,y:path==='wake'?center.y:le.y,potential:path==='wake'?x-center.x:x-le.x+1};}});
  assert.ok(requests.length>0);assert.deepEqual(outer,before.outer);
  for(const side of ['upper','lower'])assert.deepEqual(guides[0][side][2*N],before.guides[0][side][2*N]);
  for(let i=2*N+1;i<ranks.length;i++){assert.equal(guides[0].upper[i].x,guides[0].lower[i].x);assert.equal(guides[0].upper[i].y,guides[0].lower[i].y);close(guides[0].upper[i].potential-guides[0].lower[i].potential,.2);}
  assert.equal(fitted.physicalAcceptance,false);
});
