// Geodesic-loxodrome tensor interpolation
//
// Syntax: Y = tenInterpTeem(D,NUM_STEPS,w,interp_type)
//
// Inputs:
//    D - tensor data (6x8): 6 tensor coefficients at 8 cubic points
//    NUM_STEPS - number of interpolated steps between two interpolants
//    w - interpolated (x,y,z) locations (1x3) in [0 NUM_STEPS]
//    interp_type - tensor interpolation method
//      2: Ki
//      3: Ri
//
// Output:
//    Y - interpolated tensor coefficients (1x6)
//          (D11,D12,D13,D22,D23,D33)
//
// Compile:
//    mex tenInterpTeem.c -I/Users/gahm/Documents/MATLAB/teem/include -L/Users/gahm/Documents/MATLAB/teem/lib -lteem
//
// Teem Library:
//    http://teem.sourceforge.net
//
// Written by JK Gahm, UCLA. 02/06/2013.

#include <math.h> /* Needed for the ceil() prototype. */
#include "mex.h"
#include <teem/ten.h>

//#define NUM_STEPS	500		// number of steps for geolox
#define MAX_ITER	500     // max number of iterations of geolox

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  
  /* Declare variables. */
  int i,j,k,l,m,n,nzmax,*irs,*jcs,isfull;
  double *D, *w, *sr;
  double percent_sparse;
  double alpha=6E-5, alpha1;
  double T1[7], T2[7], d;
  tenInterpParm *tip;
  int num_steps, ptype, interp_type;
  double ten[8][7];
  int wx, wy, wz;
  Nrrd *npath1, *npath2, *npath3, *npath4;
  double (*lup1)(const void *, size_t), (*lup2)(const void *, size_t), (*lup3)(const void *, size_t), (*lup4)(const void *, size_t);

  /* Check for proper number of input and output arguments. */    
  if (nrhs != 4) {
    mexErrMsgTxt("Four input arguments required.\ntenInterpTeem(D,NUM_STEPS,weight,interp_type)");
  } 
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.\ntenInterpTeem(D,NUM_STEPS,weight)");
  }
  
  // input
  D = mxGetPr(prhs[0]);
  num_steps = (int)(*mxGetPr(prhs[1]));
  w = mxGetPr(prhs[2]);
  wx = (int)w[0]; wy = (int)w[1]; wz = (int)w[2];
  
  for (i=0; i < 8; i++) 
    for (j=0; j < 7; j++) {
      if (j==0)
        ten[i][j] = 12;
      else
        ten[i][j] = D[i*6+j-1];
    }
  
  interp_type = (int)(*mxGetPr(prhs[3]));
  if (interp_type == 0)
    ptype = tenInterpTypeLinear;
  else if (interp_type == 1)
    ptype = tenInterpTypeLogLinear;
  else if (interp_type == 2)
    ptype = tenInterpTypeGeoLoxK;
  else if (interp_type == 3)
    ptype = tenInterpTypeGeoLoxR;
  
  // output
  plhs[0] = mxCreateDoubleMatrix(1,6,mxREAL);
  sr  = mxGetPr(plhs[0]);
  
  // setup for interp
  tip = tenInterpParmNew();
  if ((ptype == tenInterpTypeLinear) || (ptype == tenInterpTypeLogLinear))
    tip->numSteps = num_steps+1;
  else {
    tip->numSteps = num_steps;
    tip->maxIter = MAX_ITER;
  }
  npath1 = nrrdNew();	npath2 = nrrdNew();	npath3 = nrrdNew(); npath4 = nrrdNew();
    
  // interpolation along x-axis
  tenInterpTwoDiscrete_d(npath1, ten[0], ten[1], ptype, tip->numSteps, tip);
  tenInterpTwoDiscrete_d(npath2, ten[2], ten[3], ptype, tip->numSteps, tip);
  tenInterpTwoDiscrete_d(npath3, ten[4], ten[5], ptype, tip->numSteps, tip);
  tenInterpTwoDiscrete_d(npath4, ten[6], ten[7], ptype, tip->numSteps, tip);
  lup1 = nrrdDLookup[npath1->type];
  lup2 = nrrdDLookup[npath2->type];
  lup3 = nrrdDLookup[npath3->type];
  lup4 = nrrdDLookup[npath4->type];

  for (l = 0; l < 7; l++) {
    i = 7*wx + l;
    ten[0][l] = lup1(npath1->data, i);
    ten[1][l] = lup2(npath2->data, i);
    ten[2][l] = lup3(npath3->data, i);
    ten[3][l] = lup4(npath4->data, i);
  }

  // interpolation along y-axis
  tenInterpTwoDiscrete_d(npath1, ten[0], ten[1], ptype, tip->numSteps, tip);
  tenInterpTwoDiscrete_d(npath2, ten[2], ten[3], ptype, tip->numSteps, tip);
  lup1 = nrrdDLookup[npath1->type];
  lup2 = nrrdDLookup[npath2->type];

  for (l = 0; l < 7; l++) {
    i = 7*wy + l;
    ten[0][l] = lup1(npath1->data, i);
    ten[1][l] = lup2(npath2->data, i);
  }

  // interpolate along z-axis
  tenInterpTwoDiscrete_d(npath1, ten[0], ten[1], ptype, tip->numSteps, tip);
  lup1 = nrrdDLookup[npath1->type];
  
  // final output
  for (l = 1; l < 7; l++) {
    i = 7*wz + l;
    sr[l-1] = lup1(npath1->data, i);
  }
  
  nrrdNuke(npath1);	nrrdNuke(npath2);	nrrdNuke(npath3);	nrrdNuke(npath4);
}
