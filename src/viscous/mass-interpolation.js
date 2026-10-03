// SPDX-License-Identifier: GPL-2.0-or-later
// A C1 cubic interpolant of signed displacement mass. Derivatives reproduce
// quadratic mass profiles on nonuniform arc coordinates. The two TE limits
// use one-sided stencils; the upper and lower sides are never joined.
export function nodalMassSlopes(s,start=0){
  if(s.length<3||s.some((v,i)=>!Number.isFinite(v)||(i&&v<=s[i-1])))throw new Error('Mass interpolation needs ordered arc coordinates.');
  return s.map((value,i)=>{
    const first=Math.max(0,Math.min(i-1,s.length-3)),ids=[first,first+1,first+2];
    return new Map(ids.map(j=>{
      const [k,l]=ids.filter(n=>n!==j);
      return[start+j,((value-s[k])+(value-s[l]))/((s[j]-s[k])*(s[j]-s[l]))];
    }));
  });
}

export function addMassWeights(target,source,scale=1){
  for(const [i,v]of source)target.set(i,(target.get(i)??0)+scale*v);
  return target;
}

// Exact integral of dm/ds on a fractional segment [lo,hi]. Taking a
// difference of the Hermite mass itself enforces the endpoint mass budget,
// including after arbitrary subdivision of a fixed body panel.
export function integrateMassDerivative(node,length,left,right,lo,hi){
  const basis=t=>[2*t**3-3*t*t+1,-2*t**3+3*t*t,length*(t**3-2*t*t+t),length*(t**3-t*t)];
  const a=basis(lo),b=basis(hi),w=new Map([[node,b[0]-a[0]],[node+1,b[1]-a[1]]]);
  addMassWeights(w,left,b[2]-a[2]);addMassWeights(w,right,b[3]-a[3]);return w;
}
