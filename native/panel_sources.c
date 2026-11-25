// SPDX-License-Identifier: GPL-2.0-or-later
#include "panel_sources.h"
#include <float.h>

/* Integrate log(z-xi) for the source streamfunction or 1/(z-xi) for
 * velocity. Outside two panel lengths the series avoids cancellation in
 * the recurrence for the higher moments. */
static double wrap(double angle,double center)
{
    return angle+2*ALULA_PI*floor(.5+(center-angle)/(2*ALULA_PI));
}

void alula_source_psi(const Panel *p,AlulaPoint target,double heading,double chart,int degree,double psi[3])
{
    AlulaPoint d=subtract(target,p->a),t=p->tangent;
    double x=dot(d,t)/p->length,y=cross(t,d)/p->length;
    double r0=hypot(x,y);
    double t0=wrap(atan2(y,x),chart);
    if (r0>2 && degree==0) {
        /* Exact constant-source integral. Form the angle and radius
         * differences directly so a short, distant panel keeps its sign. */
        double rr=x*x+y*y;
        double angle=atan2(-y,rr-x),log_ratio=-.5*log1p((1-2*x)/rr);
        psi[0]=t0+(x-1)*angle+y*log_ratio;
    } else if (r0>2) {
        double rr=x*x+y*y,zr=x/rr,zi=-y/rr,pr=1,pi=0;
        for (int k=0;k<=degree;++k) psi[k]=t0/(k+1);
        for (int n=1;n<=60;++n) {
            double next=pr*zr-pi*zi;
            pi=pr*zi+pi*zr; pr=next;
            for (int k=0;k<=degree;++k) psi[k]-=pi/(n*(k+n+1));
            if ((fabs(pr)+fabs(pi))/(n*(n+1))<DBL_EPSILON*.01) break;
        }
    } else if (r0<1e-14 || hypot(x-1,y)<1e-14) {
        double angle=r0<1e-14?wrap(ALULA_PI,chart):t0;
        for (int k=0;k<=degree;++k) psi[k]=angle/(k+1);
    } else {
        double t1=wrap(atan2(y,x-1),t0);
        double re=log(r0/hypot(x-1,y)), im=t0-t1;
        for (int k=0;k<=degree;++k) {
            double next=x*re-y*im-1.0/(k+1);
            im=x*im+y*re; re=next;
            psi[k]=(t1+im)/(k+1);
        }
    }
    for (int k=0;k<=degree;++k) psi[k]=p->length*(psi[k]+heading/(k+1))/(2*ALULA_PI);
}

void alula_source_velocity(const Panel *p,AlulaPoint target,int degree,int principal,
    int midpoint,AlulaPoint velocity[3])
{
    AlulaPoint d=subtract(target,p->a),t=p->tangent;
    double x=midpoint?.5:dot(d,t)/p->length,y=midpoint?0:cross(t,d)/p->length;
    double r2=x*x+y*y, re[3]={0}, im[3]={0};
    if (degree && r2>4) {
        double zr=x/r2,zi=-y/r2,pr=1,pi=0;
        for (int n=0;n<60;++n) {
            double next=pr*zr-pi*zi;
            pi=pr*zi+pi*zr; pr=next;
            for (int k=0;k<=degree;++k) { re[k]+=pr/(k+n+1); im[k]+=pi/(k+n+1); }
            if ((fabs(pr)+fabs(pi))/(n+1)<DBL_EPSILON*.001*(fabs(re[0])+fabs(im[0]))) break;
        }
    } else {
        double b2=(x-1)*(x-1)+y*y;
        if (r2<1e-28) re[0]=-log(p->length);
        else if (b2<1e-28) re[0]=log(p->length);
        else {
            re[0]=.5*log(r2/b2);
            im[0]=principal && fabs(y)<1e-12?0:-atan2(y,x*(x-1)+y*y);
        }
        for (int k=1;k<=degree;++k) {
            re[k]=x*re[k-1]-y*im[k-1]-1.0/k;
            im[k]=x*im[k-1]+y*re[k-1];
        }
    }
    for (int k=0;k<=degree;++k)
        velocity[k]=(AlulaPoint){(re[k]*t.x+im[k]*t.y)/(2*ALULA_PI),
                                  (re[k]*t.y-im[k]*t.x)/(2*ALULA_PI)};
}
