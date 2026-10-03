// SPDX-License-Identifier: GPL-2.0-or-later
import {evaluateAffine} from '../../src/potential/reconstruction.js';
import {isentropicState} from '../../src/potential/isentropic.js';

// Independent control-volume audit of the reconstructed outer field. The
// full-potential solver enforces mass, not discrete momentum. Therefore this
// defect must be measured under refinement, not included in its Newton norm.
// In freestream dynamic-pressure units the momentum-defect flux is
// 2 rho (V - V_inf) (V.n) + Cp n. Subtracting V_inf times continuity reduces
// cancellation and does not change the continuum closed-volume balance.
export function potentialMomentumAudit(system,x,{alpha=system.alpha??0}={}){
  const {mesh,reconstruction,conditions}=system,a=alpha*Math.PI/180,freestream=[Math.cos(a),Math.sin(a)];
  const groups={},rule=[[-Math.sqrt(.6),5/18],[0,4/9],[Math.sqrt(.6),5/18]];
  const reconstructed=(f,p)=>{
    const cells=f.neighbor===null?[f.owner]:[f.owner,f.neighbor],v=[0,0];
    for(const cell of cells)reconstruction.velocity(cell,p).forEach((form,k)=>v[k]+=evaluateAffine(form,x)/cells.length);
    return v;
  };
  let facePressureDifference=0;const evaluated=system.evaluate(x);
  for(const [i,f]of mesh.faces.entries()){
    if(f.neighbor!==null)continue;
    const kind=f.boundary.type,group=groups[kind]??={x:0,y:0,mass:0,pressureX:0,pressureY:0};
    const center=reconstructed(f,f),actual=system.faces[i].map(form=>evaluateAffine(form,x));
    const p=mesh.vertices[f.a],q=mesh.vertices[f.b];let pressure=0;
    for(const [s,weight]of rule){
      const point={x:f.x+s*(q.x-p.x)/2,y:f.y+s*(q.y-p.y)/2},value=reconstructed(f,point);
      // Direct normal coupling is constant along a face. Recover its exact
      // correction from the public midpoint form, independently of the
      // integrated flux implementation. Walls retain constant normal speed.
      const delta=value.map((v,k)=>v-center[k]);
      if(kind==='wall'){const normal=delta[0]*f.nx+delta[1]*f.ny;delta[0]-=normal*f.nx;delta[1]-=normal*f.ny;}
      const [u,v]=actual.map((n,k)=>n+delta[k]),state=isentropicState(u,v,conditions);
      const length=weight*f.length,mass=state.rho*(u*f.nx+v*f.ny)*length,cp=state.cp*length;
      group.x+=2*(u-freestream[0])*mass+cp*f.nx;group.y+=2*(v-freestream[1])*mass+cp*f.ny;
      group.mass+=mass;group.pressureX+=cp*f.nx;group.pressureY+=cp*f.ny;pressure+=cp;
    }
    facePressureDifference=Math.max(facePressureDifference,Math.abs(pressure-evaluated.pressureIntegrals[i]));
  }
  const total={x:0,y:0,mass:0};for(const g of Object.values(groups))for(const k of Object.keys(total))total[k]+=g[k];
  return{scope:'Outer-flow momentum defect including every wall and both sides of each wake. This is a discretization audit, not total viscous drag or an experimental error bound.',
    groups,total,momentumDefect:Math.hypot(total.x,total.y),facePressureDifference,
    dragComponent:total.x*Math.cos(a)+total.y*Math.sin(a),liftComponent:total.y*Math.cos(a)-total.x*Math.sin(a)};
}
