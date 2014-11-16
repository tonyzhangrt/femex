/*

! This file is part of AGMG software package,
! Release 3.0 built on "Aug 18 2011"
!
!    AGMG is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AGMG is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AGMG.  If not, see <http://www.gnu.org/licenses/>.
!
! Up-to-date copies of the AGMG package can be obtained
! from the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG>.
!
! You can acknowledge, citing references [1] [2], and [3], the contribution
! of this package in any scientific publication dependent upon it use.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, Report GANMN 10-03, Universite Libre de Bruxelles,
!    Brussels, Belgium, 2010.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, Report GANMN 11-01, Universite Libre de Bruxelles, Brussels,
!    Belgium, 2011.
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! AGMG Copyright (C) 2011 Yvan NOTAY
!
*/
#include "mex.h"
#include "matrix.h"
/* computational subroutine */
extern void dagmg_(int *,double *,int *,int *,double *,double *,int *,int *,
                   int *,int *,double *);

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    mwIndex *irs,*jcs;
    double  *as, *x, *f, *fc, *tolp, *iprintp, *iterp, *nrestp, *x0, *ijobp;
    double  tol;
    mwSize  nn, ncol, nsiz;
    int     n, i, k, nrest, iter, iprint, ijob, ijb, ijbe;
    static int    *ja, *ia, nz, np=0;
    static double *a; 
  
    /* get pointers to input*/
    f       = mxGetPr(prhs[1]);
    iprintp = mxGetPr(prhs[2]);
    nrestp  = mxGetPr(prhs[3]);
    iterp   = mxGetPr(prhs[4]);
    tolp    = mxGetPr(prhs[5]);
    ijobp   = mxGetPr(prhs[6]);
    iprint=(int)*iprintp;
    nrest =(int)*nrestp;
    iter  =(int)*iterp;
    tol   =*tolp;

    ijb  =(int)*ijobp;
    if (ijb >=100) {ijob=ijb-100;ijbe=ijob;}
    else {ijob=ijb; 
	if (ijb >=0) {ijbe=ijb+100;}
        else {ijbe=ijb;}
    }
    if (ijob <= 1 && np>0) {
        np=0;
    }
    if (ijob >= 0) {
  
    if (ijob < 3) {
    nn      = mxGetN(prhs[0]);
    as      = mxGetPr(prhs[0]);
    irs     = mxGetIr(prhs[0]);
    jcs     = mxGetJc(prhs[0]);
    /* create workspace, get pointers to it*/
    n   = nn;
    nz  = jcs[n];
    nsiz= nz*sizeof(int);
    ja  = mxMalloc(nsiz); 
    nsiz= (n+1)*sizeof(int);
    ia  = mxMalloc(nsiz); 
        for (i=0 ; i<=n ; i++) {
	    ia[i] = jcs[i]+1;}
        for (i=0 ; i<jcs[n] ; i++) {
            ja[i] = irs[i]+1;}
    }
    else {n=np;nn=np;}

    if (ijob != 1) {
    /* create workspace, get pointers to it*/
    nsiz= n*sizeof(double);
    fc  = mxMalloc(nsiz); 
    /* copy r.h.s. to avoid overwriting */
    for(i=0 ; i<nn ; i++){fc[i]=f[i];}
    /* create a new array and set the output pointer to it */
    ncol = 1;
    plhs[0] = mxCreateDoubleMatrix(nn, ncol, mxREAL);
    x = mxGetPr(plhs[0]);
    /* process initial guess, if any */
       if (nrhs > 7 && ijob < 3) {
        x0 = mxGetPr(prhs[7]);
	    for(i=0 ; i<nn ; i++){
		x[i]=x0[i];}
	    ijbe=ijbe+10;
       }
    }
    }

	dagmg_(&n,as,ja,ia,fc,x,&ijbe,&iprint,&nrest,&iter,&tol);
    
    if (ijob == 0 || ijob >1) {
       if (ijob == 3) {iter=0;}
       plhs[1] = mxCreateDoubleMatrix(ncol, ncol, mxREAL);
       iterp   = mxGetPr(plhs[1]);
       *iterp  = iter;
       if (iter < 0) {iter=-iter;}
       nn = iter+1;
       plhs[2] = mxCreateDoubleMatrix(nn, ncol, mxREAL);
       f = mxGetPr(plhs[2]);
       for (i=0 ; i <= iter; i++) {f[i]=fc[i];}
    }
    /* keep n in local memory if ijob == 1 */
    else if (ijob == 1) {
        np=n;
        }
}







