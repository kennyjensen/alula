// SPDX-License-Identifier: GPL-2.0-or-later
// A full-state restart fixes the mesh and wake discretization. In particular,
// changing a requested spacing must never silently reuse the previous grid.
export function assertPotentialMeshSettings(mesh,requested={}){
  if(mesh.topology!=='constrained-triangular-multielement')return;
  if((requested.type??'cross-line')!=='triangular')throw new Error('Restart mesh topology differs from the requested mesh.');
  const {padding=4,surfaceScale=1,boundaryScale=surfaceScale,growth=.4,maxVertices=20000}=requested;
  const points=mesh.contours.flat(),span=Math.max(...points.map(p=>p.x))-Math.min(...points.map(p=>p.x));
  const expected={padding,surfaceScale,boundaryScale,growth,farSpacing:requested.farSpacing??padding*span/3};
  for(const [key,value]of Object.entries(expected)){
    const actual=mesh.controls?.[key];
    if(!Number.isFinite(actual)||Math.abs(actual-value)>1e-12*Math.max(Math.abs(actual),Math.abs(value)))
      throw new Error(`Restart mesh ${key}=${actual} differs from requested ${value}. Use a BL seed without initial/meshSeed to generate the requested grid.`);
  }
  if(mesh.vertices.length>maxVertices)throw new Error('Restart mesh exceeds the requested vertex budget.');
}

export function assertPotentialRestart(input,controls,saved){
  const original=saved.config?.input,state=saved.debug??saved;
  if(!original||JSON.stringify(input.elements.map(e=>e.points))!==JSON.stringify(original.elements.map(e=>e.points)))
    throw new Error('Full-state restart requires identical element coordinates and station counts.');
  const lengths=p=>p.wakeLengths??p.elements.map(()=>p.wakeLength??2);
  if((input.wakeCount??24)!==(original.wakeCount??24)||JSON.stringify(lengths(input))!==JSON.stringify(lengths(original)))
    throw new Error('Full-state restart requires identical wake station counts and lengths.');
  if(input.wakePaths&&JSON.stringify(input.wakePaths)!==JSON.stringify(saved.wakePaths))
    throw new Error('Full-state restart would overwrite the explicitly requested wake paths.');
  if((controls.farfield??'fixed')!==(saved.config.controls.farfield??'fixed'))
    throw new Error('Full-state restart requires the same farfield unknown count.');
  assertPotentialMeshSettings(state.originalMesh,controls.mesh);
}
