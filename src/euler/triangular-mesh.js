// SPDX-License-Identifier: GPL-2.0-or-later
import cdt2d from '../../third_party/cdt2d/cdt2d.js';
import { prepareContour,validateAssembly,pointInside,segmentsTouch } from '../geometry/airfoil.js';
import { buildMesh } from './mesh.js';

const key=(a,b)=>a<b?`${a}:${b}`:`${b}:${a}`;
const triangleQuality=(cell,vertices)=>{
  let squareSum=0;cell.vertices.forEach((id,i)=>{const a=vertices[id],b=vertices[cell.vertices[(i+1)%3]];squareSum+=(a.x-b.x)**2+(a.y-b.y)**2;});
  return 4*Math.sqrt(3)*cell.area/squareSum;
};
const projection=(p,a,b)=>{
  const dx=b.x-a.x,dy=b.y-a.y,l2=dx*dx+dy*dy,t=Math.max(0,Math.min(1,((p.x-a.x)*dx+(p.y-a.y)*dy)/l2));
  return{distance:Math.hypot(p.x-a.x-t*dx,p.y-a.y-t*dy),length:Math.sqrt(l2)};
};

// Constrained triangulation preserves all body and wake segments. Background
// spacing grows with distance from those curves, avoiding the leading-edge
// skew and long thin downstream cells of the reference cross-line topology.
// Each triangle is still passed through the same conservative mesh checks.
export function triangularAirfoilMesh(contours,{wakePaths,padding=4,surfaceScale=1,boundaryScale=surfaceScale,growth=.4,farSpacing,
  maxVertices=20000}={}){
  if(![padding,surfaceScale,boundaryScale,growth].every(Number.isFinite)||padding<=0||surfaceScale<=0||boundaryScale<=0||growth<=0
    ||!Number.isInteger(maxVertices)||maxVertices<100)throw new Error('Invalid triangular mesh controls.');
  const prepared=contours.map(c=>prepareContour(c));validateAssembly(prepared);
  if(wakePaths&&wakePaths.length!==prepared.length)throw new Error('Supply one triangular-mesh wake per body.');
  const bodyPoints=prepared.flat(),xmin=Math.min(...bodyPoints.map(p=>p.x)),xmax=Math.max(...bodyPoints.map(p=>p.x)),chord=xmax-xmin;
  const bounds={xmin:xmin-padding*chord,xmax:xmax+padding*chord,ymin:Math.min(...bodyPoints.map(p=>p.y))-padding*chord,ymax:Math.max(...bodyPoints.map(p=>p.y))+padding*chord};
  const coarse=farSpacing??padding*chord/3;
  if(!(coarse>0)||!Number.isFinite(coarse))throw new Error('Invalid farfield mesh spacing.');
  const vertices=[],edges=[],tags=[],curves=[],ids=new Map(),wakeVertices=[],wakeInterpolants=[];
  const add=p=>{
    const k=`${p.x.toPrecision(15)},${p.y.toPrecision(15)}`;
    if(ids.has(k))return ids.get(k);
    if(vertices.length>=maxVertices)throw new Error('Triangular mesh vertex limit.');
    const index=vertices.length;vertices.push({x:p.x,y:p.y});ids.set(k,index);return index;
  };
  const edge=(a,b,tag,spacing=true)=>{
    if(a===b)throw new Error('Collapsed mesh constraint.');
    if(spacing)curves.push([vertices[a],vertices[b]]);
    // Refine the outer discretization along the unchanged polygon, as well
    // as normal to it. Keeping a full original edge while shrinking the
    // adjacent cell height does not form a uniform mesh-refinement sequence.
    const steps=spacing?Math.max(1,Math.ceil(1/boundaryScale)):1,p=vertices[a],q=vertices[b];
    let previous=a;
    for(let j=1;j<=steps;j++){
      const fraction=j/steps,next=j===steps?b:add({x:p.x+fraction*(q.x-p.x),y:p.y+fraction*(q.y-p.y)});
      edges.push([previous,next]);tags.push(tag);previous=next;
      if(j<steps&&tag.type==='wake-cut')wakeInterpolants.push({vertex:next,element:tag.element,segment:tag.segment,fraction});
    }
  };
  prepared.forEach((points,element)=>{
    const loop=points.slice(0,-1).map(add);
    loop.forEach((a,i)=>edge(a,loop[(i+1)%loop.length],{type:'wall',element}));
  });
  const corners=[{x:bounds.xmin,y:bounds.ymin},{x:bounds.xmax,y:bounds.ymin},{x:bounds.xmax,y:bounds.ymax},{x:bounds.xmin,y:bounds.ymax}];
  const perimeter=[];
  corners.forEach((a,i)=>{
    const b=corners[(i+1)%4],steps=Math.ceil(Math.hypot(b.x-a.x,b.y-a.y)/coarse);
    for(let j=0;j<steps;j++)perimeter.push(add({x:a.x+(b.x-a.x)*j/steps,y:a.y+(b.y-a.y)*j/steps}));
  });
  perimeter.forEach((a,i)=>edge(a,perimeter[(i+1)%perimeter.length],{type:'farfield'},false));
  const allWakes=[];
  (wakePaths??[]).forEach((path,element)=>{
    if(path.length<2||path.some(p=>!Number.isFinite(p.x)||!Number.isFinite(p.y))
      ||Math.hypot(path[0].x-prepared[element][0].x,path[0].y-prepared[element][0].y)>chord*1e-10
      ||path.some((p,i)=>i>0&&p.x<=path[i-1].x))throw new Error('Triangular mesh needs downstream x-monotone wakes from each trailing edge.');
    const points=path.map(p=>({...p}));
    if(points.some(p=>p.x>=bounds.xmax||p.y<=bounds.ymin||p.y>=bounds.ymax))throw new Error('Modeled wake leaves the mesh domain.');
    // dm/ds vanishes at the modeled wake endpoint; its constant-mass
    // downstream continuation has no source. The continuous scalar potential
    // needs no artificial cut from that point to the farfield.
    const list=points.map(add);wakeVertices.push(list);allWakes.push(points);
    list.slice(1).forEach((b,i)=>edge(list[i],b,{type:'wake-cut',element,segment:i}));
  });
  // Reject crossings explicitly. The mesher must not silently merge wakes
  // or create a fluid passage through a body.
  for(let i=0;i<edges.length;i++)for(let j=i+1;j<edges.length;j++){
    const a=edges[i],b=edges[j];if(a.some(v=>b.includes(v)))continue;
    if(tags[i].type!=='wake-cut'&&tags[j].type!=='wake-cut')continue;
    if(tags[i].type==='farfield'||tags[j].type==='farfield')continue;
    if(segmentsTouch(vertices[a[0]],vertices[a[1]],vertices[b[0]],vertices[b[1]],chord*1e-12))throw new Error('Wake constraints cross a body or another wake.');
  }
  const spacing=p=>{
    let h=coarse;
    for(const [a,b] of curves){const d=projection(p,a,b);h=Math.min(h,surfaceScale*d.length+growth*d.distance);}
    return h;
  };
  const candidates=[],visit=(x0,y0,x1,y1,depth=0)=>{
    const p={x:.5*(x0+x1),y:.5*(y0+y1)},h=spacing(p);
    if(Math.hypot(x1-x0,y1-y0)>1.5*h){
      if(depth>30)throw new Error('Unresolved triangular-mesh length scale.');
      visit(x0,y0,p.x,p.y,depth+1);visit(p.x,y0,x1,p.y,depth+1);visit(x0,p.y,p.x,y1,depth+1);visit(p.x,p.y,x1,y1,depth+1);
    }else if(!prepared.some(body=>pointInside(p,body)))candidates.push({...p,h});
    if(candidates.length>4*maxVertices)throw new Error('Triangular mesh candidate limit.');
  };
  visit(bounds.xmin,bounds.ymin,bounds.xmax,bounds.ymax);
  candidates.sort((a,b)=>a.h-b.h||a.x-b.x||a.y-b.y);
  for(const p of candidates){
    if(curves.some(([a,b])=>projection(p,a,b).distance<.35*p.h))continue;
    if(vertices.some(v=>Math.hypot(p.x-v.x,p.y-v.y)<.55*p.h))continue;
    add(p);
  }
  // Split every T-junction before invoking cdt2d; it requires a clean planar
  // graph, including shared trailing-edge vertices.
  const splitEdges=[],splitTags=[];
  edges.forEach(([a,b],i)=>{
    const p=vertices[a],q=vertices[b],dx=q.x-p.x,dy=q.y-p.y,l2=dx*dx+dy*dy,chain=[{id:a,t:0},{id:b,t:1}];
    vertices.forEach((v,id)=>{
      if(id===a||id===b)return;const t=((v.x-p.x)*dx+(v.y-p.y)*dy)/l2;
      if(t>0&&t<1&&Math.abs((v.x-p.x)*dy-(v.y-p.y)*dx)<1e-12*chord*Math.sqrt(l2))chain.push({id,t});
    });
    chain.sort((a,b)=>a.t-b.t);for(let j=1;j<chain.length;j++){splitEdges.push([chain[j-1].id,chain[j].id]);splitTags.push(tags[i]);}
  });
  const triangles=cdt2d(vertices.map(p=>[p.x,p.y]),splitEdges).filter(t=>{
    const p={x:(vertices[t[0]].x+vertices[t[1]].x+vertices[t[2]].x)/3,y:(vertices[t[0]].y+vertices[t[1]].y+vertices[t[2]].y)/3};
    return !prepared.some(body=>pointInside(p,body));
  });
  const boundaries=[];splitEdges.forEach(([a,b],i)=>{if(splitTags[i].type!=='wake-cut')boundaries.push({a,b,...splitTags[i]});});
  const mesh=buildMesh(vertices,triangles,boundaries),byEdge=new Map(mesh.faces.map((f,i)=>[key(f.a,f.b),i])),cuts=[];
  splitEdges.forEach(([a,b],i)=>{
    if(splitTags[i].type!=='wake-cut')return;const face=byEdge.get(key(a,b));
    if(face===undefined||mesh.faces[face].neighbor===null)throw new Error('Triangulation lost a connected wake constraint.');
    cuts.push({a,b,face,...splitTags[i]});
  });
  return{...mesh,contours:prepared,cuts,bounds,wakeVertices,wakeInterpolants,wakePaths:allWakes,modeledWakeCounts:(wakePaths??[]).map(w=>w.length),topology:'constrained-triangular-multielement',
    deformationReferenceQuality:mesh.cells.map(c=>triangleQuality(c,mesh.vertices)),
    controls:{padding,surfaceScale,boundaryScale,growth,farSpacing:coarse,maxVertices}};
}

// Keep the triangulation and background points fixed during small wake
// updates. This makes the discrete geometry residual continuous and permits
// a potential warm start. Reject any deformation that inverts a cell; the
// caller can then generate and validate a new triangulation explicitly.
// Also reject severe quality loss before cells invert; positive area alone
// does not prevent an inaccurate velocity reconstruction in a sliver cell.
export function deformTriangularWakes(seed,wakePaths){
  if(seed.topology!=='constrained-triangular-multielement'||wakePaths.length!==seed.wakeVertices.length)return null;
  const vertices=seed.vertices.map(p=>({...p}));
  for(let e=0;e<wakePaths.length;e++){
    const path=wakePaths[e],ids=seed.wakeVertices[e],count=seed.modeledWakeCounts[e];
    if(path.length!==count)return null;
    for(let j=0;j<count;j++)vertices[ids[j]]={...path[j]};
  }
  for(const {vertex,element,segment,fraction} of seed.wakeInterpolants??[]){
    const a=wakePaths[element][segment],b=wakePaths[element][segment+1];
    vertices[vertex]={x:a.x+fraction*(b.x-a.x),y:a.y+fraction*(b.y-a.y)};
  }
  const boundaries=seed.faces.filter(f=>f.neighbor===null).map(f=>({a:f.a,b:f.b,...f.boundary}));
  let mesh;
  try{mesh=buildMesh(vertices,seed.cells.map(c=>c.vertices),boundaries);}catch{return null;}
  const reference=seed.deformationReferenceQuality??seed.cells.map(c=>triangleQuality(c,seed.vertices));
  if(mesh.cells.some((c,i)=>triangleQuality(c,vertices)<.25*reference[i]))return null;
  return{...seed,...mesh,deformationReferenceQuality:reference,wakePaths:seed.wakeVertices.map(ids=>ids.map(i=>({...vertices[i]})))};
}
