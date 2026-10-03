import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {quadCoupledNcritResult} from '../src/ui/quad-coupled-ncrit.js';

test('Run and legacy continuation requests both use direct final-condition assembly', async () => {
 const key='__directWorker',messages=[],calls=[];
 const modules={
  '../euler/streamtube-coupled-assembly.js':{solveCoupledStreamtubeAssembly(input,options){
   calls.push({input:structuredClone(input),direct:options.direct});
   return {sourceCase:input};
  }},
  '../ui/quad-coupled-result.js':{quadCoupledResultForDisplay:r=>r},
  '../ui/quad-mesh-progress.js':{createQuadMeshProgress:()=>({stage(){},mesh:m=>m})},
 };
 const self={postMessage:m=>messages.push(m)};
 globalThis[key]={self,modules,quadCoupledNcritResult};
 const source=fs.readFileSync(new URL('../src/worker/solver.js',import.meta.url),'utf8')
  .replace(/import \{([^}]+)\} from '[^']+';/g,(_,names)=>`const {${names.replace(/ as /g, ': ')}} = globalThis.${key};`)
  .replace(/await import\('([^']+)'\)/g,(_,path)=>`globalThis.${key}.modules[${JSON.stringify(path)}]`);
 try{
  await import('data:text/javascript;base64,'+Buffer.from(`const self=globalThis.${key}.self;\n${source}`).toString('base64'));
  const caseData={flowModel:'streamtube-grid',quadBoundaryLayers:true,mach:.74,alpha:2.68,
   gridIntervals:128,gridTubes:11,ncrit:4};
  for(const task of ['solve','continue-coupled'])
   await self.onmessage({data:{id:1,task,caseData,parentResult:{sourceCase:{mach:.2,alpha:0}}}});
  assert.equal(messages.filter(m=>m.type==='result').length,2,JSON.stringify(messages));
  assert.equal(calls.length,2);
  for(const call of calls){assert.deepEqual(call.input,caseData);assert.equal(call.direct,true);}
 }finally{delete globalThis[key];}
});
