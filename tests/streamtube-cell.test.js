import test from 'node:test';
import assert from 'node:assert/strict';
import {streamtubeCellGeometry,streamtubeSection,evaluateStreamtubeCell} from '../src/euler/streamtube-cell.js';

const close=(a,b,tol=2e-12)=>assert.ok(Math.abs(a-b)<tol*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
const dot=(a,b)=>a.x*b.x+a.y*b.y;
const lower=[{x:0,y:0},{x:.4,y:0},{x:1,y:0}],upper=lower.map(p=>({x:p.x+.13,y:.2}));

test('intrinsic streamtube geometry keeps signed areas and rejects folded streamlines',()=>{
  const g=streamtubeCellGeometry(lower,upper);
  g.normalAreas.forEach(a=>close(a,.2));close(g.pressureCurvature,0);
  assert.throws(()=>streamtubeCellGeometry(upper,lower),/Folded/);
  assert.throws(()=>streamtubeCellGeometry(lower,lower),/Folded/);
  assert.throws(()=>streamtubeCellGeometry(lower.slice(1),upper),/three finite/);
  assert.throws(()=>streamtubeCellGeometry(lower,[upper[0],{x:-2,y:.2},upper[2]]),/Folded/);
});

test('streamtube elimination preserves mass and energy and differentiates all section inputs',()=>{
  const parameters={density:1.1,massFlow:.19,normalArea:.23,stagnationEnthalpy:9,gamma:1.4};
  const s=streamtubeSection(parameters);
  close(s.rho*s.q*parameters.normalArea,parameters.massFlow);
  close(parameters.gamma/(parameters.gamma-1)*s.p/s.rho+.5*s.q*s.q,parameters.stagnationEnthalpy);
  for(const key of ['density','massFlow','normalArea','stagnationEnthalpy']){
    const h=1e-5*parameters[key],p=streamtubeSection({...parameters,[key]:parameters[key]+h}),m=streamtubeSection({...parameters,[key]:parameters[key]-h});
    for(const variable of ['q','p'])close(s.derivatives[variable][key],(p[variable]-m[variable])/(2*h),2e-9);
  }
  assert.throws(()=>streamtubeSection({...parameters,density:0}),/Invalid/);
  assert.throws(()=>streamtubeSection({...parameters,stagnationEnthalpy:.01}),/Nonpositive/);
});

test('physical thermal rejection retains the eliminated energy state without classifying invalid inputs as recoverable',()=>{
  const parameters={density:1.02,massFlow:11.5,normalArea:1,stagnationEnthalpy:63,gamma:1.4};
  assert.throws(()=>streamtubeSection(parameters),error=>{
    assert.equal(error.code,'streamtube-static-enthalpy');
    const d=error.diagnostics;
    assert.equal(d.rho,parameters.density);
    for(const key of ['massFlow','normalArea','stagnationEnthalpy','gamma'])assert.equal(d[key],parameters[key]);
    close(d.rho*d.q*d.normalArea,d.massFlow);
    close(d.enthalpy+.5*d.q*d.q,d.stagnationEnthalpy);
    assert.ok(d.enthalpy<0);assert.ok(d.q>d.maximumPhysicalSpeed);
    return true;
  });
  for(const invalid of [{density:0},{normalArea:-1},{stagnationEnthalpy:NaN}])
    assert.throws(()=>streamtubeSection({...parameters,...invalid}),error=>{
      assert.equal(error.code,undefined);assert.match(error.message,/Invalid/);return true;
    });
});

test('uniform streamtubes are invariant under rotation, translation and geometric scale',()=>{
  const rho=1.2,q=.9,p=3,gamma=1.4,h0=gamma/(gamma-1)*p/rho+.5*q*q;
  for(const angle of [0,.43,-2.1])for(const length of [.1,1,4]){
    const map=({x,y})=>({x:3+length*(x*Math.cos(angle)-y*Math.sin(angle)),y:-2+length*(x*Math.sin(angle)+y*Math.cos(angle))});
    const c=evaluateStreamtubeCell({lower:lower.map(map),upper:upper.map(map),densities:[rho,rho],massFlow:rho*q*.2*length,stagnationEnthalpy:h0,gamma});
    close(c.streamwiseResidual,0);close(c.interfacePressure.lower,p);close(c.interfacePressure.upper,p);
    c.states.forEach(s=>{close(s.q,q);close(s.p,p);});
  }
});

test('eliminated side pressures recover the independent conservative vector momentum balance on a curved skew cell',()=>{
  const l=[{x:0,y:0},{x:.45,y:.025},{x:1.1,y:.08}],u=[{x:.1,y:.25},{x:.58,y:.3},{x:1.18,y:.4}];
  for(const factor of [0,.1,.2]){
    const c=evaluateStreamtubeCell({lower:l,upper:u,densities:[1.05,1.02],massFlow:.27,stagnationEnthalpy:8,pressureCorrectionFactor:factor});
    const {geometry:g,states:[a,b],interfacePressure:pi}=c;
    // Integrate outward momentum directly. Upper/lower bent faces each
    // have constant interface pressure; their normal measures telescope.
    const flux={x:0,y:0};
    const pressure=(p,v,sign)=>{flux.x+=sign*p*v.y;flux.y-=sign*p*v.x;};
    pressure(a.p,g.sections[0],-1);pressure(b.p,g.sections[1],1);
    pressure(pi.lower,g.sides.lower,1);pressure(pi.upper,g.sides.upper,-1);
    for(const key of ['x','y'])flux[key]+=.27*(b.q*g.directions[1][key]-a.q*g.directions[0][key]);
    close(dot(flux,g.transverse)/g.area,0);
    close(dot(flux,g.streamwise)/g.area,c.streamwiseResidual);
    close(pi.lower+pi.upper,a.p+b.p+2*c.pressureCorrection);
  }
});

test('the undissipated streamtube cell preserves exact normal-shock mass, momentum and enthalpy jumps',()=>{
  const gamma=1.4,rho=1,p=1,mach=2,q=mach*Math.sqrt(gamma*p/rho);
  const ratio=(gamma+1)*mach*mach/((gamma-1)*mach*mach+2),pressureRatio=1+2*gamma/(gamma+1)*(mach*mach-1);
  const straightUpper=lower.map(v=>({...v,y:.2}));
  const c=evaluateStreamtubeCell({lower,upper:straightUpper,densities:[rho,rho*ratio],massFlow:rho*q*.2,stagnationEnthalpy:gamma/(gamma-1)*p/rho+.5*q*q,gamma});
  close(c.states[0].p,p);close(c.states[1].p,p*pressureRatio);close(c.states[1].q,q/ratio);
  close(c.streamwiseResidual,0);
  assert.ok(Math.log(pressureRatio)-gamma*Math.log(ratio)>0);
  // This local jump identity does not establish stable shock capture.
});

test('intrinsic side-pressure recovery converges to the independent compressible irrotational vortex',()=>{
  const gamma=1.4,mach=.3,pInf=1/(gamma*mach*mach),h0=1/((gamma-1)*mach*mach)+.5;
  const exact=r=>{const q=1/r,t=1+.5*(gamma-1)*mach*mach*(1-q*q),rho=t**(1/(gamma-1));return{q,rho,p:pInf*t**(gamma/(gamma-1))};};
  const errors=[];
  for(const h of [.2,.1,.05,.025]){
    const ri=1-h/2,ro=1+h/2,point=(r,t)=>({x:r*Math.cos(t),y:r*Math.sin(t)});
    const l=[-h,0,h].map(t=>point(ro,t)),u=[-h,0,h].map(t=>point(ri,t));
    // Independent Simpson integration of rho*q through a radial section.
    const n=1000,dr=(ro-ri)/n;let integral=0;
    for(let i=0;i<=n;i++){const s=exact(ri+i*dr);integral+=(i===0||i===n?1:i%2?4:2)*s.rho*s.q;}
    const massFlow=integral*dr/3,rho=exact(Math.cos(h/2)).rho;
    const c=evaluateStreamtubeCell({lower:l,upper:u,densities:[rho,rho],massFlow,stagnationEnthalpy:h0,gamma});
    close(c.streamwiseResidual,0);
    errors.push(Math.max(Math.abs(c.interfacePressure.lower-exact(ro).p),Math.abs(c.interfacePressure.upper-exact(ri).p)));
  }
  for(let i=1;i<errors.length;i++)assert.ok(errors[i]<.3*errors[i-1],String(errors));
  assert.ok(errors.at(-1)<1e-3,String(errors));
});
