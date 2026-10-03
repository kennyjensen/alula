import test from 'node:test';import assert from 'node:assert/strict';
import {quadCoupledCoefficients} from '../src/ui/quad-coupled-coefficients.js';
const near=(a,b,tol=1e-12)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
const fixture=()=>({alpha:0,mach:0,referenceChord:1,momentReference:{x:.25,y:0},stagnation:[{x:0,y:0}],
 surfaces:[{element:0,side:'upper',stations:[{x:.5,y:.1,cp:1},{x:1,y:0,cp:1}]},
 {element:0,side:'lower',stations:[{x:.5,y:-.1,cp:1},{x:1,y:0,cp:1}]}],
 wakes:[{element:0,stations:[{index:0,theta:.01,deltaStar:.02,ue:1}]}]});
test('uniform pressure gives no force/moment; recovered wake drag is twice momentum thickness',()=>{
 const f=fixture(),r=quadCoupledCoefficients(f);
 for(const key of ['cx','cy','cl','cm','pressureIntegralDrag'])near(r[key],0);
 near(r.cd,.02);assert.equal(r.physicalAcceptance,false);
});
test('known symmetric triangular upper pressure load has correct lift, moment and reference scaling',()=>{
 const f=fixture();f.surfaces[0].stations[0].cp=2;
 const r=quadCoupledCoefficients(f);near(r.cx,0);near(r.cl,-.5);near(r.cm,.125);near(r.cd,.02);
 const shifted=quadCoupledCoefficients({...f,momentReference:{x:.35,y:0}});near(shifted.cm,r.cm+.1*r.cy);
 const ref=quadCoupledCoefficients({...f,referenceChord:2});near(ref.cl,r.cl/2);near(ref.cm,r.cm/4);near(ref.cd,r.cd/2);
 const rotated=quadCoupledCoefficients({...f,alpha:90});near(rotated.cl,0);near(rotated.pressureIntegralDrag,-.5);near(rotated.cd,.02);
});
test('separate wake losses sum once; invalid wake leaves pressure forces available but total drag unavailable',()=>{
 const f=fixture();f.stagnation.push({x:2,y:0});
 f.surfaces.push(...f.surfaces.map(s=>({...s,element:1,stations:s.stations.map(p=>({...p,x:p.x+2}))})));
 f.wakes.push({element:1,stations:[{index:1,theta:.005,deltaStar:.01,ue:1}]});
 let r=quadCoupledCoefficients(f);near(r.cd,.03);assert.equal(r.wakes.length,2);
 f.wakes[0].stations[0].theta=NaN;r=quadCoupledCoefficients(f);
 assert.equal(r.cd,null);assert.ok(Number.isFinite(r.cl));assert.ok(Number.isFinite(r.cm));assert.equal(r.warnings.length,1);
});
