// SPDX-License-Identifier: GPL-2.0-or-later
#include "contour.h"
#include "panel_internal.h"

#include <stdlib.h>
#include <string.h>

static AlulaPoint scale(AlulaPoint p, double s) { return (AlulaPoint){s*p.x,s*p.y}; }
static AlulaPoint add(AlulaPoint a, AlulaPoint b) { return (AlulaPoint){a.x+b.x,a.y+b.y}; }

void alula_contour_spline_free(AlulaContourSpline *curve)
{
    free(curve->arc);
    free(curve->points);
    free(curve->second);
    *curve = (AlulaContourSpline){0};
}

static int create_spline(AlulaContourSpline *curve,const AlulaPoint *points,int n,int parabolic)
{
    if (!curve || !points || n < 4 || n > ALULA_PANEL_MAX_TOTAL+1) return 0;
    *curve = (AlulaContourSpline){.count=n};
    curve->arc = calloc((size_t)n,sizeof *curve->arc);
    curve->points = malloc((size_t)n*sizeof *curve->points);
    curve->second = calloc((size_t)n,sizeof *curve->second);
    double *diag = malloc((size_t)n*sizeof *diag);
    if (!curve->arc || !curve->points || !curve->second || !diag) goto fail;
    memcpy(curve->points,points,(size_t)n*sizeof *points);
    for (int i=1; i<n; ++i) {
        AlulaPoint delta = subtract(points[i],points[i-1]);
        curve->arc[i] = curve->arc[i-1]+hypot(delta.x,delta.y);
        if (!isfinite(curve->arc[i]) || !(curve->arc[i]>curve->arc[i-1])) goto fail;
    }
    for (int i=1; i<n-1; ++i) {
        double l=curve->arc[i]-curve->arc[i-1], r=curve->arc[i+1]-curve->arc[i];
        diag[i]=2*(l+r);
        if (parabolic && i==1) diag[i]+=l;
        if (parabolic && i==n-2) diag[i]+=r;
        AlulaPoint f=scale(subtract(scale(subtract(points[i+1],points[i]),1/r),
            scale(subtract(points[i],points[i-1]),1/l)),6);
        if (i>1) {
            double a=l/diag[i-1];
            diag[i]-=a*l;
            f=subtract(f,scale(curve->second[i-1],a));
        }
        curve->second[i]=f;
    }
    for (int i=n-2; i>0; --i)
        curve->second[i]=scale(subtract(curve->second[i],
            scale(curve->second[i+1],curve->arc[i+1]-curve->arc[i])),1/diag[i]);
    if (parabolic) {
        curve->second[0]=curve->second[1];
        curve->second[n-1]=curve->second[n-2];
    }
    free(diag);
    return 1;
fail:
    free(diag);
    alula_contour_spline_free(curve);
    return 0;
}

int alula_contour_spline_create(AlulaContourSpline *curve,const AlulaPoint *points,int n)
{
    return create_spline(curve,points,n,0);
}

int alula_contour_spline_parabolic(AlulaContourSpline *curve,const AlulaPoint *points,int n)
{
    return create_spline(curve,points,n,1);
}

static int interval(const AlulaContourSpline *s,double parameter)
{
    int lo=0,hi=s->count-1;
    while(hi-lo>1) {
        int mid=lo+(hi-lo)/2;
        if(s->arc[mid]>parameter) hi=mid; else lo=mid;
    }
    return lo;
}

static AlulaPoint tangent(AlulaPoint slope,AlulaPoint d0,AlulaPoint d1,double h,double a,double b)
{
    return add(slope,scale(add(scale(d0,1-3*a*a),scale(d1,3*b*b-1)),h/6));
}

void alula_contour_spline_value(const AlulaContourSpline *s, double arc,
    AlulaPoint *p, AlulaPoint *t, AlulaPoint *dd)
{
    int lo=interval(s,arc),hi=lo+1;
    double h=s->arc[hi]-s->arc[lo], b=(arc-s->arc[lo])/h,a=1-b;
    AlulaPoint d0=s->second[lo],d1=s->second[hi];
    AlulaPoint shape=add(scale(d0,a*a*a-a),scale(d1,b*b*b-b));
    if (p) *p=add(add(scale(s->points[lo],a),scale(s->points[hi],b)),scale(shape,h*h/6));
    if (t) *t=tangent(scale(subtract(s->points[hi],s->points[lo]),1/h),d0,d1,h,a,b);
    if (dd) *dd=add(scale(d0,a),scale(d1,b));
}

static double integrate_arc(const AlulaContourSpline *s,int segment,double end)
{
    static const double nodes[]={.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363};
    static const double weights[]={.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763};
    double start=s->arc[segment],h=s->arc[segment+1]-start;
    double sum=0,mid=.5*(start+end),half=.5*(end-start);
    AlulaPoint slope=scale(subtract(s->points[segment+1],s->points[segment]),1/h);
    AlulaPoint d0=s->second[segment],d1=s->second[segment+1];
    for (int i=0;i<4;++i) {
        double l=((mid-half*nodes[i])-start)/h,r=((mid+half*nodes[i])-start)/h;
        AlulaPoint lo=tangent(slope,d0,d1,h,1-l,l),hi=tangent(slope,d0,d1,h,1-r,r);
        sum+=weights[i]*(hypot(lo.x,lo.y)+hypot(hi.x,hi.y));
    }
    return half*sum;
}

void alula_contour_arc_build(const AlulaContourSpline *s,double *a)
{
    a[0]=0;
    for (int i=1;i<s->count;++i) a[i]=a[i-1]+integrate_arc(s,i-1,s->arc[i]);
}

double alula_contour_arc(const AlulaContourSpline *s,const double *a,double parameter)
{
    if (!(parameter>=0 && parameter<=s->arc[s->count-1])) return NAN;
    int lo=interval(s,parameter);
    return a[lo]+integrate_arc(s,lo,parameter);
}
