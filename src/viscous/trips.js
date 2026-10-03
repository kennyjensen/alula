// SPDX-License-Identifier: GPL-2.0-or-later
import { splind,sinvrt } from './xfoil/spline.js';

// XIFSET's geometric inversion before its upstream-of-stagnation fallback.
// Keep the actual material trip location when stagnation moves during the
// assembly solve; an initially ineffective trip can become downstream later.
export function geometricTripArc(ctx,side){
  const fraction=ctx.XSTRIP[side],end=ctx.S[side===1?1:ctx.N];
  if(fraction>=1)return end;
  const chx=ctx.XTE-ctx.XLE,chy=ctx.YTE-ctx.YLE,chord2=chx*chx+chy*chy;
  const s=ctx.S.subarray(1,ctx.N+1),x=Float64Array.from(s,(_,i)=>((ctx.X[i+1]-ctx.XLE)*chx+(ctx.Y[i+1]-ctx.YLE)*chy)/chord2),derivative=new Float64Array(ctx.N);
  splind(x,derivative,s,ctx.N,-999,-999);
  return sinvrt(ctx.SLE+(end-ctx.SLE)*fraction,fraction,x,derivative,s,ctx.N);
}

export function tripDistance(surface,stagnation){
  const direction=surface.side===1?-1:1;
  const end=surface.side===1?surface.body.s[0]:surface.body.s.at(-1);
  const distance=direction*(surface.tripArc-stagnation),te=direction*(end-stagnation);
  return distance>=0?Math.min(distance,te):te;
}

export function surfaceTransition(system,surface,states,stagnation){
  const id=surface.transitionId,tripS=tripDistance(surface,stagnation);
  const input=surface.transition===0?system.leadingTransitionInput(states[id],tripS)
    :{upstream:states[surface.ids[surface.transition-1]],downstream:states[id],regime:'transition',tripS};
  const transition=system.kernel.interval(input).transition;
  const arc=stagnation+(surface.side===1?-1:1)*transition.s,{body}=surface;
  let i=0;while(i<body.s.length-2&&body.s[i+1]<arc)i++;
  const f=(arc-body.s[i])/(body.s[i+1]-body.s[i]),a=body.points[i],b=body.points[i+1];
  return{...transition,point:{x:a.x+f*(b.x-a.x),y:a.y+f*(b.y-a.y)}};
}
