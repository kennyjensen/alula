// SPDX-License-Identifier: GPL-2.0-or-later
// AssemblyScript port of paired boundary SLOR. Float64 throughout; no fast math.
// The JS reference defines the equations and tests this kernel independently.
let nx:i32=0, nt:i32=0, ng:i32=0, capacity:i32=0;
let mode:i32=1;
let physical:bool=false, grape:bool=false, lengthScale:f64=1;
let nodes=new Float64Array(0), trial=new Float64Array(0), xi=new Float64Array(0), eta=new Float64Array(0);
let background=new Float64Array(0), weights=new Float64Array(0), active=new Float64Array(0), corners=new Float64Array(0);
let dx=new Float64Array(0), de=new Float64Array(0), groups=new Float64Array(0);
let controls=new Float64Array(0), met=new Float64Array(0), cellSines=new Float64Array(0);
let lower=new Float64Array(0), diagonal=new Float64Array(0), upper=new Float64Array(0), rhs=new Float64Array(0), scales=new Float64Array(0);
let tu=new Float64Array(0), tr=new Float64Array(0), solution=new Float64Array(0);
let lu=new Float64Array(0), rowScale=new Float64Array(0), pivots=new Int32Array(0), value=new Float64Array(0);
let events=new Float64Array(0), summaries=new Float64Array(0);
export let eventCount:i32=0, errorCode:i32=0, errorStation:i32=0, errorSide:i32=0;
export let residual:f64=0, merit:f64=0, minSine:f64=0, maxUpdate:f64=0;
export function allocate(x:i32,t:i32,g:i32,maxRows:i32,controlMode:i32,upwind:bool,scale:f64):void {
  nx=x;nt=t;ng=g;capacity=2*maxRows;mode=controlMode;physical=mode==2||mode==4;grape=upwind;lengthScale=scale;
  const points=(nx+1)*(nt+1), n=nx-1, b=capacity*capacity;
  nodes=new Float64Array(2*points);trial=new Float64Array(2*points);
  xi=new Float64Array(nx+1);eta=new Float64Array(nt+1);background=new Float64Array(points);
  weights=new Float64Array(2*(nt+1));active=new Float64Array(2*(nx+1));corners=new Float64Array(2*(nx+1));
  dx=new Float64Array(6*(nx+1));de=new Float64Array(6*(nt+1));groups=new Float64Array(2*ng);
  controls=new Float64Array(10*(nx+1));met=new Float64Array(5*points);cellSines=new Float64Array(nx*nt);
  lower=new Float64Array(n*b);diagonal=new Float64Array(n*b);upper=new Float64Array(n*b);
  rhs=new Float64Array(n*capacity);scales=new Float64Array(n*capacity);tu=new Float64Array(n*b);tr=new Float64Array(n*capacity);solution=new Float64Array(n*capacity);
  lu=new Float64Array(b);rowScale=new Float64Array(capacity);pivots=new Int32Array(capacity);value=new Float64Array(capacity);
  events=new Float64Array(ng*21*5);summaries=new Float64Array(ng*6);
}
// Stable ABI: coordinates, xi, eta, background, boundary weights, activity,
// corner flags, xi/eta derivative stencils, row groups, trial log, row log.
export function pointer(which:i32):usize {
  switch(which){case 0:return nodes.dataStart;case 1:return xi.dataStart;case 2:return eta.dataStart;case 3:return background.dataStart;
    case 4:return weights.dataStart;case 5:return active.dataStart;case 6:return corners.dataStart;case 7:return dx.dataStart;
    case 8:return de.dataStart;case 9:return groups.dataStart;case 10:return events.dataStart;case 11:return summaries.dataStart;}
  return 0;
}
function pos(i:i32,j:i32):i32{return i*(nt+1)+j;}
function point(a:Float64Array,i:i32,j:i32,k:i32):f64{return a[2*pos(i,j)+k];}
function first(a:Float64Array,i:i32,j:i32,k:i32,streamwise:bool):f64 {
  const p=point(a,i,j,k);
  return streamwise ? dx[6*i]*(point(a,i-1,j,k)-p)+dx[6*i+2]*(point(a,i+1,j,k)-p)
    : de[6*j]*(point(a,i,j-1,k)-p)+de[6*j+2]*(point(a,i,j+1,k)-p);
}
function second(a:Float64Array,i:i32,j:i32,k:i32,streamwise:bool):f64 {
  const p=point(a,i,j,k);
  return streamwise ? dx[6*i+3]*(point(a,i-1,j,k)-p)+dx[6*i+5]*(point(a,i+1,j,k)-p)
    : de[6*j+3]*(point(a,i,j-1,k)-p)+de[6*j+5]*(point(a,i,j+1,k)-p);
}
function mixed(a:Float64Array,i:i32,j:i32,k:i32):f64 {
  const h=first(a,i,j,k,false);
  return dx[6*i]*(first(a,i-1,j,k,false)-h)+dx[6*i+2]*(first(a,i+1,j,k,false)-h);
}
function sf(i:i32,di:i32,d:f64):f64 {
  if(!grape)return dx[6*i+di+1];
  if(d<0)return di==1?0:(di==0?1:-1)/(xi[i]-xi[i-1]);
  return di==-1?0:(di==0?-1:1)/(xi[i+1]-xi[i]);
}
function sourceFirst(a:Float64Array,i:i32,j:i32,k:i32,d:f64):f64 {
  const p=point(a,i,j,k);
  return sf(i,-1,d)*(point(a,i-1,j,k)-p)+sf(i,1,d)*(point(a,i+1,j,k)-p);
}
function sample(a:Float64Array,side:i32,i:i32):bool {
  const j=side==0?0:nt, sign=side==0?1:-1;
  const l=xi[i]-xi[i-1],r=xi[i+1]-xi[i], width=l+r;
  const px=point(a,i,j,0),py=point(a,i,j,1);
  const lx=px-point(a,i-1,j,0),rx=point(a,i+1,j,0)-px,ly=py-point(a,i-1,j,1),ry=point(a,i+1,j,1)-py;
  const tx=(lx+rx)/width,ty=(ly+ry)/width,sx=2*(rx/r-lx/l)/width,sy=2*(ry/r-ly/l)/width;
  const gamma=tx*tx+ty*ty, speed=Math.sqrt(gamma), normalX=-sign*ty/speed,normalY=sign*tx/speed;
  const d1=sign*(eta[j+sign]-eta[j]),d2=sign*(eta[j+2*sign]-eta[j]);
  const A=2*d2/(d1*d1*(d2-d1)),B=-2*d1/(d2*d2*(d2-d1)),m=2*(1/d1+1/d2);
  const kx=A*(point(a,i,j+sign,0)-px)+B*(point(a,i,j+2*sign,0)-px),ky=A*(point(a,i,j+sign,1)-py)+B*(point(a,i,j+2*sign,1)-py);
  const curvature=(sx*normalX+sy*normalY)/gamma,kn=kx*normalX+ky*normalY,disc=1-4*curvature*(kn/m)/m;
  if(!(kn>0)||!(disc>64*f64.EPSILON)||!isFinite(disc)||!isFinite(kn)||!isFinite(m)){errorCode=1;errorStation=i;errorSide=side;return false;}
  const speedN=(2*kn/m)/(1+Math.sqrt(disc)),alpha=speedN*speedN,tangentK=kx*tx+ky*ty;
  // Keep reference operation order (dot divided by alpha, not reciprocal multiplication).
  let stretch=-(sx*tx+sy*ty)/gamma-tangentK/alpha;
  if(!(speedN>0)||!isFinite(speedN)||!isFinite(alpha)||!isFinite(stretch)
    ||!isFinite(kx-m*speedN*normalX)||!isFinite(ky-m*speedN*normalY)){errorCode=2;return false;}
  const slope=m*Math.sqrt(disc),gx=-tx/alpha+2*(tangentK/alpha)*normalX/speedN/slope,gy=-ty/alpha+2*(tangentK/alpha)*normalY/speedN/slope;
  const c=5*(side*(nx+1)+i);
  if(physical)stretch/=gamma;
  controls[c]=stretch;
  controls[c+1]=physical?A*gx/gamma:A*gx;controls[c+2]=physical?A*gy/gamma:A*gy;
  controls[c+3]=physical?B*gx/gamma:B*gx;controls[c+4]=physical?B*gy/gamma:B*gy;
  for(let k=0;k<5;k++)if(!isFinite(controls[c+k])){errorCode=2;return false;}
  return true;
}
function boundaryControls(a:Float64Array):bool {
  for(let side=0;side<2;side++)for(let i=1;i<nx;i++){
    const ix=side*(nx+1)+i;
    const needed=active[ix]!=0&&corners[ix]==0 || i>1&&active[ix-1]!=0&&corners[ix-1]!=0 || i<nx-1&&active[ix+1]!=0&&corners[ix+1]!=0;
    if(needed&&!sample(a,side,i))return false;
  }
  return true;
}
function feedback(i:i32,j:i32):f64 {
  let result=background[pos(i,j)];
  for(let side=0;side<2;side++){
    const ix=side*(nx+1)+i;
    if(active[ix]==0)continue;
    const control=corners[ix]!=0?.5*(controls[5*(ix-1)]+controls[5*(ix+1)]):controls[5*ix];
    result+=weights[2*j+side]*(control-background[pos(i,side==0?0:nt)]);
  }
  return result;
}
function controlDerivative(i:i32,j:i32,k:i32,di:i32,v:i32):f64 {
  let result:f64=0;
  for(let side=0;side<2;side++){
    const ix=side*(nx+1)+i,distance=side==0?k:nt-k;
    if(active[ix]==0||distance<1||distance>2)continue;
    if(corners[ix]!=0){if(di!=0)result+=.5*weights[2*j+side]*controls[5*(ix+di)+1+2*(distance-1)+v];}
    else if(di==0)result+=weights[2*j+side]*controls[5*ix+1+2*(distance-1)+v];
  }
  return result;
}
function metrics(a:Float64Array,firstRow:i32=1,lastRow:i32=nt-1,refreshControls:bool=true):bool {
  if(refreshControls&&(mode==1||mode==2)&&!boundaryControls(a))return false;
  for(let i=1;i<nx;i++)for(let j=firstRow;j<=lastRow;j++){
    const gx=first(a,i,j,0,true),gy=first(a,i,j,1,true),hx=first(a,i,j,0,false),hy=first(a,i,j,1,false);
    const alpha=hx*hx+hy*hy,beta=gx*hx+gy*hy,gamma=gx*gx+gy*gy,J=gx*hy-gy*hx;
    const drift=(physical?J*J:alpha)*(mode==0?0:mode>=3?background[pos(i,j)]:feedback(i,j)),p=5*pos(i,j);
    if(!(alpha>0&&gamma>0)||J==0||!isFinite(J)||!isFinite(alpha)||!isFinite(beta)||!isFinite(gamma)||!isFinite(drift)){errorCode=2;errorStation=i;return false;}
    met[p]=alpha;met[p+1]=beta;met[p+2]=gamma;met[p+3]=J;met[p+4]=drift;
  }
  return true;
}
function operator(a:Float64Array,i:i32,j:i32,k:i32):f64 {
  const p=5*pos(i,j);
  return met[p]*second(a,i,j,k,true)-2*met[p+1]*mixed(a,i,j,k)+met[p+2]*second(a,i,j,k,false)+met[p+4]*sourceFirst(a,i,j,k,met[p+4]);
}
function quality(a:Float64Array,firstCell:i32=0,lastCell:i32=nt-1):bool {
  minSine=Infinity;
  for(let i=0;i<nx;i++)for(let j=0;j<nt;j++){
    let cellMinimum=cellSines[i*nt+j];
    if(j>=firstCell&&j<=lastCell){
      cellMinimum=Infinity;
      for(let k=0;k<4;k++){
        const p=k==0?pos(i,j):k==1?pos(i+1,j):k==2?pos(i+1,j+1):pos(i,j+1);
        const q=k==0?pos(i+1,j):k==1?pos(i+1,j+1):k==2?pos(i,j+1):pos(i,j);
        const r=k==0?pos(i+1,j+1):k==1?pos(i,j+1):k==2?pos(i,j):pos(i+1,j);
        const ax=a[2*q]-a[2*p],ay=a[2*q+1]-a[2*p+1],bx=a[2*r]-a[2*q],by=a[2*r+1]-a[2*q+1];
        const sine=(ax*by-ay*bx)/(Math.hypot(ax,ay)*Math.hypot(bx,by));
        cellMinimum=Math.min(cellMinimum,isFinite(sine)?sine:-Infinity);
      }
      cellSines[i*nt+j]=cellMinimum;
    }
    minSine=Math.min(minSine,cellMinimum);
  }
  return minSine>1e-12;
}
export function evaluate():bool {
  errorCode=0;
  if(!quality(nodes)){errorCode=3;return false;}
  if(!metrics(nodes))return false;
  residual=0;let sum:f64=0;
  for(let i=1;i<nx;i++)for(let j=1;j<nt;j++){
    const scale=(met[5*pos(i,j)]+met[5*pos(i,j)+2])*lengthScale;
    const x=operator(nodes,i,j,0)/scale,y=operator(nodes,i,j,1)/scale;
    residual=Math.max(residual,Math.max(Math.abs(x),Math.abs(y)));sum=sum+x*x+y*y;
  }
  merit=.5*sum;
  if(!isFinite(residual)||!isFinite(merit)){errorCode=2;return false;}
  return true;
}
function assemble(start:i32,count:i32):void {
  const size=2*count,b=size*size;
  for(let i=1;i<nx;i++)for(let jr=0;jr<count;jr++){
    const j=start+jr,p=5*pos(i,j),alpha=met[p],beta=met[p+1],gamma=met[p+2],J=met[p+3],drift=met[p+4];
    for(let r=0;r<2;r++){rhs[(i-1)*size+2*jr+r]=-operator(nodes,i,j,r);scales[(i-1)*size+2*jr+r]=(alpha+gamma)*lengthScale;}
    for(let kr=0;kr<count;kr++){
      const k=start+kr,delta=k-j,w1=Math.abs(delta)<=1?de[6*j+delta+1]:0,w2=Math.abs(delta)<=1?de[6*j+delta+4]:0;
      for(let di=-1;di<=1;di++)for(let v=0;v<2;v++){
        const dg=k==j?dx[6*i+di+1]:0,dh=di==0?w1:0;
        const g=first(nodes,i,j,v,true),h=first(nodes,i,j,v,false),da=2*h*dh,db=h*dg+g*dh,dc=2*g*dg;
        const dJ=v==0?first(nodes,i,j,1,false)*dg-first(nodes,i,j,1,true)*dh:-first(nodes,i,j,0,false)*dg+first(nodes,i,j,0,true)*dh;
        const dDrift=(physical?2*drift/J*dJ:drift/alpha*da)+(physical?J*J:alpha)*controlDerivative(i,j,k,di,v);
        const scalar=alpha*dx[6*i+di+4]*(k==j?1:0)-2*beta*dx[6*i+di+1]*w1+gamma*w2*(di==0?1:0)+drift*sf(i,di,drift)*(k==j?1:0);
        for(let r=0;r<2;r++){
          const val=(r==v?scalar:0)+da*second(nodes,i,j,r,true)-2*db*mixed(nodes,i,j,r)+dc*second(nodes,i,j,r,false)+dDrift*sourceFirst(nodes,i,j,r,drift);
          const index=(i-1)*b+(2*jr+r)*size+2*kr+v;
          if(di==-1)lower[index]=val;else if(di==0)diagonal[index]=val;else upper[index]=val;
        }
      }
    }
  }
}
function factor(size:i32):bool {
  for(let r=0;r<size;r++){
    let scale:f64=0;for(let c=0;c<size;c++)scale=Math.max(scale,Math.abs(lu[r*size+c]));
    if(!(scale>0)||!isFinite(scale))return false;rowScale[r]=scale;
    for(let c=0;c<size;c++)lu[r*size+c]/=scale;
  }
  for(let k=0;k<size;k++){
    let pivot=k;for(let i=k+1;i<size;i++)if(Math.abs(lu[i*size+k])>Math.abs(lu[pivot*size+k]))pivot=i;
    if(Math.abs(lu[pivot*size+k])<64*f64.EPSILON)return false;pivots[k]=pivot;
    if(pivot!=k)for(let j=0;j<size;j++){const tmp=lu[k*size+j];lu[k*size+j]=lu[pivot*size+j];lu[pivot*size+j]=tmp;}
    for(let i=k+1;i<size;i++){const f=lu[i*size+k]/lu[k*size+k];lu[i*size+k]=f;if(f!=0)for(let j=k+1;j<size;j++)lu[i*size+j]-=f*lu[k*size+j];}
  }
  return true;
}
function solve(size:i32):bool {
  for(let i=0;i<size;i++)value[i]/=rowScale[i];
  for(let k=0;k<size;k++)if(pivots[k]!=k){const tmp=value[k];value[k]=value[pivots[k]];value[pivots[k]]=tmp;}
  for(let i=0;i<size;i++)for(let j=0;j<i;j++)value[i]-=lu[i*size+j]*value[j];
  for(let i=size-1;i>=0;i--){for(let j=i+1;j<size;j++)value[i]-=lu[i*size+j]*value[j];value[i]/=lu[i*size+i];if(!isFinite(value[i]))return false;}
  return true;
}
function blockSolve(size:i32):bool {
  const n=nx-1,b=size*size;
  for(let i=0;i<n;i++){
    for(let k=0;k<b;k++)lu[k]=diagonal[i*b+k];
    for(let r=0;r<size;r++)tr[i*size+r]=rhs[i*size+r];
    if(i>0)for(let r=0;r<size;r++)for(let k=0;k<size;k++){
      const a=lower[i*b+r*size+k];tr[i*size+r]-=a*tr[(i-1)*size+k];
      for(let c=0;c<size;c++)lu[r*size+c]-=a*tu[(i-1)*b+k*size+c];
    }
    if(!factor(size))return false;
    for(let r=0;r<size;r++)value[r]=tr[i*size+r];if(!solve(size))return false;
    for(let r=0;r<size;r++)tr[i*size+r]=value[r];
    if(i+1<n)for(let c=0;c<size;c++){
      for(let r=0;r<size;r++)value[r]=upper[i*b+r*size+c];if(!solve(size))return false;
      for(let r=0;r<size;r++)tu[i*b+r*size+c]=value[r];
    }
  }
  for(let i=n-1;i>=0;i--)for(let r=0;r<size;r++){
    let v=tr[i*size+r];if(i+1<n)for(let c=0;c<size;c++)v-=tu[i*b+r*size+c]*solution[(i+1)*size+c];
    if(!isFinite(v))return false;solution[i*size+r]=v;
  }
  return true;
}
// One call per complete sweep; no JS calls in the numerical inner loops.
export function sweep(omega:f64,tolerance:f64):bool {
  errorCode=0;eventCount=0;maxUpdate=0;
  if(!metrics(nodes))return false;quality(nodes);
  for(let g=0;g<ng;g++){
    const start=<i32>groups[2*g],count=<i32>groups[2*g+1],size=2*count;
    assemble(start,count);
    const globalControl=start<=2||start+count-1>=nt-2;
    let lineResidual:f64=0,base:f64=0;
    for(let i=0;i<nx-1;i++){let sum:f64=0;for(let k=0;k<size;k++){const v=rhs[i*size+k]/scales[i*size+k];lineResidual=Math.max(lineResidual,Math.abs(v));sum+=v*v;}base+=sum;}
    base*=.5;summaries[6*g]=base;summaries[6*g+5]=lineResidual;
    if(lineResidual<=tolerance){summaries[6*g+4]=1;continue;}summaries[6*g+4]=0;
    if(!blockSolve(size)){errorCode=4;return false;}
    let accepted=false;
    for(let halving=0;halving<=20;halving++){
      const fraction=omega*Math.pow(2,-halving);trial.set(nodes);let movement:f64=0;
      for(let i=1;i<nx;i++)for(let r=0;r<count;r++)for(let k=0;k<2;k++){
        const delta=fraction*solution[(i-1)*size+2*r+k];trial[2*pos(i,start+r)+k]+=delta;movement=Math.max(movement,Math.abs(delta)/lengthScale);
      }
      let rejection=0,trialMerit:f64=Infinity;
      if(!quality(trial,start-1,start+count-1))rejection=3;
      else if(!metrics(trial,globalControl?1:max<i32>(1,start-1),globalControl?nt-1:min<i32>(nt-1,start+count),globalControl))rejection=errorCode;
      else {trialMerit=0;for(let i=1;i<nx;i++)for(let r=0;r<count;r++)for(let k=0;k<2;k++){
        const v=operator(trial,i,start+r,k)/scales[(i-1)*size+2*r+k];trialMerit+=.5*(v*v);
      }}
      const sufficient=rejection==0&&trialMerit<=(1-2e-4*fraction)*base;
      const e=5*eventCount++;events[e]=g;events[e+1]=fraction;events[e+2]=trialMerit;events[e+3]=rejection;events[e+4]=sufficient?1:0;
      if(sufficient){nodes.set(trial);maxUpdate=Math.max(maxUpdate,movement);summaries[6*g+1]=trialMerit;summaries[6*g+2]=fraction;summaries[6*g+3]=minSine;accepted=true;break;}
    }
    if(!accepted){errorCode=5;return false;}
  }
  errorCode=0;return true;
}

// Fixed-boundary Giles scalar SLOR, used by the harmonic/spacing fallback.
export function scalarSweep(omega:f64):bool {
  errorCode=0;maxUpdate=0;const n=nx-1;
  if(!metrics(nodes))return false;
  for(let j=1;j<nt;j++){
    if(!metrics(nodes,j,j))return false;
    for(let i=1;i<nx;i++){
      const k=i-1,p=5*pos(i,j),a=met[p],b=met[p+1],c=met[p+2],d=met[p+4],b1=de[6*j+1],b2=de[6*j+4];
      lower[k]=a*dx[6*i+3]-2*b*dx[6*i]*b1+d*sf(i,-1,d);
      upper[k]=a*dx[6*i+5]-2*b*dx[6*i+2]*b1+d*sf(i,1,d);
      diagonal[k]=a*dx[6*i+4]-2*b*dx[6*i+1]*b1+c*b2+d*sf(i,0,d);
      rhs[2*k]=-operator(nodes,i,j,0);rhs[2*k+1]=-operator(nodes,i,j,1);
    }
    for(let k=0;k<n;k++){
      const scale=Math.abs(lower[k])+Math.abs(diagonal[k])+Math.abs(upper[k]);
      if(!(Math.abs(diagonal[k])>32*f64.EPSILON*scale)){errorCode=4;return false;}
      if(k+1<n){const f=lower[k+1]/diagonal[k];diagonal[k+1]-=f*upper[k];for(let v=0;v<2;v++)rhs[2*(k+1)+v]-=f*rhs[2*k+v];}
    }
    for(let k=n-1;k>=0;k--)for(let v=0;v<2;v++){
      rhs[2*k+v]=(rhs[2*k+v]-(k==n-1?0:upper[k]*rhs[2*(k+1)+v]))/diagonal[k];
      const delta=omega*rhs[2*k+v];nodes[2*pos(k+1,j)+v]+=delta;maxUpdate=Math.max(maxUpdate,Math.abs(delta)/lengthScale);
    }
  }
  return true;
}
