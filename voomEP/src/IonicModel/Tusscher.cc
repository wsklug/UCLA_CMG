#include "Tusscher.h"

namespace voom {
  void Tusscher::ComputeRates() {
    _Algebraic[7] = 1.00000/(1.00000+(exp(((_State[0]+20.0000)/7.00000))));
    _Algebraic[20] =  1102.50*(exp((- (pow((_State[0]+27.0000), 2.00000))/225.000)))+200.000/(1.00000+(exp(((13.0000 - _State[0])/10.0000))))+180.000/(1.00000+(exp(((_State[0]+30.0000)/10.0000))))+20.0000;
    _Rates[12] = (_Algebraic[7] - _State[12])/_Algebraic[20];
    _Algebraic[8] = 0.670000/(1.00000+(exp(((_State[0]+35.0000)/7.00000))))+0.330000;
    _Algebraic[21] =  562.000*(exp((- (pow((_State[0]+27.0000), 2.00000))/240.000)))+31.0000/(1.00000+(exp(((25.0000 - _State[0])/10.0000))))+80.0000/(1.00000+(exp(((_State[0]+30.0000)/10.0000))));
    _Rates[13] = (_Algebraic[8] - _State[13])/_Algebraic[21];
    _Algebraic[9] = 0.600000/(1.00000+(pow((_State[10]/0.0500000), 2.00000)))+0.400000;
    _Algebraic[22] = 80.0000/(1.00000+(pow((_State[10]/0.0500000), 2.00000)))+2.00000;
    _Rates[14] = (_Algebraic[9] - _State[14])/_Algebraic[22];
    _Algebraic[10] = 1.00000/(1.00000+(exp(((_State[0]+20.0000)/5.00000))));
    _Algebraic[23] =  85.0000*(exp((- (pow((_State[0]+45.0000), 2.00000))/320.000)))+5.00000/(1.00000+(exp(((_State[0] - 20.0000)/5.00000))))+3.00000;
    _Rates[15] = (_Algebraic[10] - _State[15])/_Algebraic[23];
    _Algebraic[11] = 1.00000/(1.00000+(exp(((20.0000 - _State[0])/6.00000))));
    _Algebraic[24] =  9.50000*(exp((- (pow((_State[0]+40.0000), 2.00000))/1800.00)))+0.800000;
    _Rates[16] = (_Algebraic[11] - _State[16])/_Algebraic[24];
    _Algebraic[0] = 1.00000/(1.00000+(exp(((- 26.0000 - _State[0])/7.00000))));
    _Algebraic[13] = 450.000/(1.00000+(exp(((- 45.0000 - _State[0])/10.0000))));
    _Algebraic[26] = 6.00000/(1.00000+(exp(((_State[0]+30.0000)/11.5000))));
    _Algebraic[34] =  1.00000*_Algebraic[13]*_Algebraic[26];
    _Rates[4] = (_Algebraic[0] - _State[4])/_Algebraic[34];
    _Algebraic[1] = 1.00000/(1.00000+(exp(((_State[0]+88.0000)/24.0000))));
    _Algebraic[14] = 3.00000/(1.00000+(exp(((- 60.0000 - _State[0])/20.0000))));
    _Algebraic[27] = 1.12000/(1.00000+(exp(((_State[0] - 60.0000)/20.0000))));
    _Algebraic[35] =  1.00000*_Algebraic[14]*_Algebraic[27];
    _Rates[5] = (_Algebraic[1] - _State[5])/_Algebraic[35];
    _Algebraic[2] = 1.00000/(1.00000+(exp(((- 5.00000 - _State[0])/14.0000))));
    _Algebraic[15] = 1400.00/ pow((1.00000+(exp(((5.00000 - _State[0])/6.00000)))), 1.0 / 2);
    _Algebraic[28] = 1.00000/(1.00000+(exp(((_State[0] - 35.0000)/15.0000))));
    _Algebraic[36] =  1.00000*_Algebraic[15]*_Algebraic[28]+80.0000;
    _Rates[6] = (_Algebraic[2] - _State[6])/_Algebraic[36];
    _Algebraic[3] = 1.00000/(pow((1.00000+(exp(((- 56.8600 - _State[0])/9.03000)))), 2.00000));
    _Algebraic[16] = 1.00000/(1.00000+(exp(((- 60.0000 - _State[0])/5.00000))));
    _Algebraic[29] = 0.100000/(1.00000+(exp(((_State[0]+35.0000)/5.00000))))+0.100000/(1.00000+(exp(((_State[0] - 50.0000)/200.000))));
    _Algebraic[37] =  1.00000*_Algebraic[16]*_Algebraic[29];
    _Rates[7] = (_Algebraic[3] - _State[7])/_Algebraic[37];
    _Algebraic[4] = 1.00000/(pow((1.00000+(exp(((_State[0]+71.5500)/7.43000)))), 2.00000));
    _Algebraic[17] = (_State[0]<- 40.0000 ?  0.0570000*(exp((- (_State[0]+80.0000)/6.80000))) : 0.00000);
    _Algebraic[30] = (_State[0]<- 40.0000 ?  2.70000*(exp(( 0.0790000*_State[0])))+ 310000.*(exp(( 0.348500*_State[0]))) : 0.770000/( 0.130000*(1.00000+(exp(((_State[0]+10.6600)/- 11.1000))))));
    _Algebraic[38] = 1.00000/(_Algebraic[17]+_Algebraic[30]);
    _Rates[8] = (_Algebraic[4] - _State[8])/_Algebraic[38];
    _Algebraic[5] = 1.00000/(pow((1.00000+(exp(((_State[0]+71.5500)/7.43000)))), 2.00000));
    _Algebraic[18] = (_State[0]<- 40.0000 ? (( ( - 25428.0*(exp(( 0.244400*_State[0]))) -  6.94800e-06*(exp(( - 0.0439100*_State[0]))))*(_State[0]+37.7800))/1.00000)/(1.00000+(exp(( 0.311000*(_State[0]+79.2300))))) : 0.00000);
    _Algebraic[31] = (_State[0]<- 40.0000 ? ( 0.0242400*(exp(( - 0.0105200*_State[0]))))/(1.00000+(exp(( - 0.137800*(_State[0]+40.1400))))) : ( 0.600000*(exp(( 0.0570000*_State[0]))))/(1.00000+(exp(( - 0.100000*(_State[0]+32.0000))))));
    _Algebraic[39] = 1.00000/(_Algebraic[18]+_Algebraic[31]);
    _Rates[9] = (_Algebraic[5] - _State[9])/_Algebraic[39];
    _Algebraic[6] = 1.00000/(1.00000+(exp(((- 8.00000 - _State[0])/7.50000))));
    _Algebraic[19] = 1.40000/(1.00000+(exp(((- 35.0000 - _State[0])/13.0000))))+0.250000;
    _Algebraic[32] = 1.40000/(1.00000+(exp(((_State[0]+5.00000)/5.00000))));
    _Algebraic[40] = 1.00000/(1.00000+(exp(((50.0000 - _State[0])/20.0000))));
    _Algebraic[42] =  1.00000*_Algebraic[19]*_Algebraic[32]+_Algebraic[40];
    _Rates[11] = (_Algebraic[6] - _State[11])/_Algebraic[42];
    _Algebraic[55] = (( (( _Constants[21]*_Constants[10])/(_Constants[10]+_Constants[22]))*_State[2])/(_State[2]+_Constants[23]))/(1.00000+ 0.124500*(exp((( - 0.100000*_State[0]*_Constants[2])/( _Constants[0]*_Constants[1]))))+ 0.0353000*(exp((( - _State[0]*_Constants[2])/( _Constants[0]*_Constants[1])))));
    _Algebraic[25] =  (( _Constants[0]*_Constants[1])/_Constants[2])*(log((_Constants[11]/_State[2])));

    _Algebraic[50] =  _Constants[16]*(pow(_State[7], 3.00000))*_State[8]*_State[9]*(_State[0] - _Algebraic[25]);
    _Algebraic[51] =  _Constants[17]*(_State[0] - _Algebraic[25]);
    _Algebraic[56] = ( _Constants[24]*( (exp((( _Constants[27]*_State[0]*_Constants[2])/( _Constants[0]*_Constants[1]))))*(pow(_State[2], 3.00000))*_Constants[12] -  (exp((( (_Constants[27] - 1.00000)*_State[0]*_Constants[2])/( _Constants[0]*_Constants[1]))))*(pow(_Constants[11], 3.00000))*_State[3]*_Constants[26]))/( ((pow(_Constants[29], 3.00000))+(pow(_Constants[11], 3.00000)))*(_Constants[28]+_Constants[12])*(1.00000+ _Constants[25]*(exp((( (_Constants[27] - 1.00000)*_State[0]*_Constants[2])/( _Constants[0]*_Constants[1]))))));
    _Rates[2] =  (( - 1.00000*(_Algebraic[50]+_Algebraic[51]+ 3.00000*_Algebraic[55]+ 3.00000*_Algebraic[56]))/( 1.00000*_Constants[4]*_Constants[2]))*_Constants[3];
    _Algebraic[33] =  (( _Constants[0]*_Constants[1])/_Constants[2])*(log((_Constants[10]/_State[1])));
    _Algebraic[44] = 0.100000/(1.00000+(exp(( 0.0600000*((_State[0] - _Algebraic[33]) - 200.000)))));
    _Algebraic[45] = ( 3.00000*(exp(( 0.000200000*((_State[0] - _Algebraic[33])+100.000))))+(exp(( 0.100000*((_State[0] - _Algebraic[33]) - 10.0000)))))/(1.00000+(exp(( - 0.500000*(_State[0] - _Algebraic[33])))));
    _Algebraic[46] = _Algebraic[44]/(_Algebraic[44]+_Algebraic[45]);
    _Algebraic[47] =  _Constants[13]*_Algebraic[46]* pow((_Constants[10]/5.40000), 1.0 / 2)*(_State[0] - _Algebraic[33]);
    _Algebraic[54] =  _Constants[20]*_State[16]*_State[15]*(_State[0] - _Algebraic[33]);
    _Algebraic[48] =  _Constants[14]* pow((_Constants[10]/5.40000), 1.0 / 2)*_State[4]*_State[5]*(_State[0] - _Algebraic[33]);
    _Algebraic[41] =  (( _Constants[0]*_Constants[1])/_Constants[2])*(log(((_Constants[10]+ _Constants[9]*_Constants[11])/(_State[1]+ _Constants[9]*_State[2]))));
    _Algebraic[49] =  _Constants[15]*(pow(_State[6], 2.00000))*(_State[0] - _Algebraic[41]);
    _Algebraic[52] = ( (( _Constants[18]*_State[11]*_State[12]*_State[13]*_State[14]*4.00000*(_State[0] - 15.0000)*(pow(_Constants[2], 2.00000)))/( _Constants[0]*_Constants[1]))*( 0.250000*_State[10]*(exp((( 2.00000*(_State[0] - 15.0000)*_Constants[2])/( _Constants[0]*_Constants[1])))) - _Constants[12]))/((exp((( 2.00000*(_State[0] - 15.0000)*_Constants[2])/( _Constants[0]*_Constants[1])))) - 1.00000);
    _Algebraic[43] =  (( 0.500000*_Constants[0]*_Constants[1])/_Constants[2])*(log((_Constants[12]/_State[3])));
    _Algebraic[53] =  _Constants[19]*(_State[0] - _Algebraic[43]);
    _Algebraic[58] = ( _Constants[32]*(_State[0] - _Algebraic[33]))/(1.00000+(exp(((25.0000 - _State[0])/5.98000))));
    _Algebraic[57] = ( _Constants[30]*_State[3])/(_State[3]+_Constants[31]);
    // dV/dt
    _Rates[0] = - (_Algebraic[47]+_Algebraic[54]+_Algebraic[48]+_Algebraic[49]+_Algebraic[52]+_Algebraic[55]+_Algebraic[50]+_Algebraic[51]+_Algebraic[56]+_Algebraic[53]+_Algebraic[58]+_Algebraic[57]+_Algebraic[12]);

    _Rates[1] =  (( - 1.00000*((_Algebraic[47]+_Algebraic[54]+_Algebraic[48]+_Algebraic[49]+_Algebraic[58]+_Algebraic[12]) -  2.00000*_Algebraic[55]))/( 1.00000*_Constants[4]*_Constants[2]))*_Constants[3];
    _Algebraic[59] = _Constants[44]/(1.00000+(pow(_Constants[42], 2.00000))/(pow(_State[3], 2.00000)));
    _Algebraic[60] =  _Constants[43]*(_State[17] - _State[3]);
    _Algebraic[61] =  _Constants[41]*(_State[10] - _State[3]);
    _Algebraic[63] = 1.00000/(1.00000+( _Constants[45]*_Constants[46])/(pow((_State[3]+_Constants[46]), 2.00000)));
    _Rates[3] =  _Algebraic[63]*((( (_Algebraic[60] - _Algebraic[59])*_Constants[51])/_Constants[4]+_Algebraic[61]) - ( 1.00000*((_Algebraic[53]+_Algebraic[57]) -  2.00000*_Algebraic[56])*_Constants[3])/( 2.00000*1.00000*_Constants[4]*_Constants[2]));
    _Algebraic[62] = _Constants[38] - (_Constants[38] - _Constants[39])/(1.00000+(pow((_Constants[37]/_State[17]), 2.00000)));
    _Algebraic[65] =  _Constants[34]*_Algebraic[62];
    _Rates[18] =  - _Algebraic[65]*_State[10]*_State[18]+ _Constants[36]*(1.00000 - _State[18]);
    _Algebraic[64] = _Constants[33]/_Algebraic[62];
    _Algebraic[66] = ( _Algebraic[64]*(pow(_State[10], 2.00000))*_State[18])/(_Constants[35]+ _Algebraic[64]*(pow(_State[10], 2.00000)));
    _Algebraic[67] =  _Constants[40]*_Algebraic[66]*(_State[17] - _State[10]);
    _Algebraic[68] = 1.00000/(1.00000+( _Constants[47]*_Constants[48])/(pow((_State[17]+_Constants[48]), 2.00000)));
    _Rates[17] =  _Algebraic[68]*(_Algebraic[59] - (_Algebraic[67]+_Algebraic[60]));
    _Algebraic[69] = 1.00000/(1.00000+( _Constants[49]*_Constants[50])/(pow((_State[10]+_Constants[50]), 2.00000)));
    _Rates[10] =  _Algebraic[69]*((( - 1.00000*_Algebraic[52]*_Constants[3])/( 2.00000*1.00000*_Constants[52]*_Constants[2])+( _Algebraic[67]*_Constants[51])/_Constants[52]) - ( _Algebraic[61]*_Constants[4])/_Constants[52]);
  }
  
  /*!
    Update as suggested by Whiteley. We have equations of the form
    \f[
    \frac{du_i}{dt} = a_i(V^n) + b_i(V^m)u_i
    \f]
    which will be solved explicitly using backward difference. Other nonlinear
    ODE's will be solved as stated below
    \f[
    \frac{d}{dt} State[i] = f( State[0],.. State[k])
    \f]
    We rewrite this as
    \f[
    State^{n+1}[i] - State^n[i] - \Delta t * f() = 0
    \f]
    We will use Newton Raphson on these ODE's. Since the ODE's are very 
    nonlinear we will use numerical derivatives to construct the Hessian.

  */
  void Tusscher::UpdateStateVariables(const Real dt) {
    /* First update gating equations of the form du_i/dt = a_i + b_i u_i
       We have 19 gating variables. We can neglect 0 which is Voltage. Hence
       we totally have 18 variables left. We have States 4,5,6,7,8,9,11,12,13,
       14, 15 and 16 of this form
    */
    int ind[]={4,5,6,7,8,9,11,12,13,14,15,16}; // 12 values
    const int nLinear = 12; // Variables needing Linear Solution
    const int nVar = 6; // Variables needing Newton Raphson
    int dind[]={34,35,36,37,38,39,42,20,21,22,23,24};// Index for denominator
    int NRind[]={1,2,3,10,17,18}; // 6 variable
    const Real TOLER = 1e-6;
    Real resid = 100., fxph[nVar], fxmh[nVar], oldRates[nVar], StateN[nVar];;
    Real Jacobian[nVar*nVar], B[nVar];
    int N = nVar, NRHS = 1, LDA = nVar, LDB=nVar;
    int IPIV[nVar], INFO;
    
    for(int i = 0; i < nVar; i++)  
      StateN[i] = _State[ NRind[i] ];

    ComputeRates();
    for(int i = 0; i < nLinear; i++) {
      Real b = -1./_Algebraic[ dind[i] ];
      Real a = _Algebraic[i]/ _Algebraic[ dind[i] ];
      // Forward Euler Update
      _State[ ind[i] ] = (_State[ ind[i] ] + dt *a)/
	(1 - dt* b);
    }

    // For rest of the state variables
    ComputeRates();
    for(int i = 0; i < nVar; i++)  
      oldRates[i] = _Rates[ NRind[i] ];  

    while ( resid > TOLER ) {
      /* We will use numerical derivatives to get the Hessian. Vary each
	 state variable and see how each rate changes
	 We have _Rates[i] = f(_State[1], _State[2] .._State[6], u^(n+1) )
	 _State*[i] - f() *dt - _State[i] = g() = 0
	 Use NR on g()
      */ 
      for(int i = 0; i < nVar; i++) {
	const Real eps = 1e-3*_State[ NRind[i] ];
	const Real twoeps = 2. * eps;
	_State[ NRind[i] ] += eps;
	ComputeRates();
	for(int j = 0; j < nVar; j++) 
	  fxph[j] =_State[ NRind[j] ] - StateN[j] - dt*_Rates[ NRind[j] ];
	
	_State[ NRind[i] ] -= 2*eps;
	ComputeRates();
	for(int j = 0; j < nVar; j++) 
	  fxmh[j] =_State[ NRind[j] ] - StateN[j] - dt*_Rates[ NRind[j] ];
	
	_State[ NRind[i] ] += eps;
	// Fill in terms for the Jacobian
	for(int j = 0; j < nVar; j++) {
	  Jacobian[j*nVar + i] = (fxph[j] - fxmh[j])/twoeps;
	}
      } // End of Jacobian
      
      for(int i = 0; i < nVar; i++) 
	B[i] = _State[ NRind[i] ] - StateN[i] - dt* oldRates[i];
      
      // Solve J^-1 f(x_n)
      dgesv_(&N, &NRHS, Jacobian, &LDA, IPIV, B, &LDB, &INFO);
      if (INFO != 0 ){
	std::cout << "DGESV_: Info = " << INFO << std::endl;
      }
      // Update Solution
      for(int i = 0; i < nVar; i++) _State[ NRind[i] ] -= B[i];
      ComputeRates();
      // Store Old Rates
      for(int i = 0; i < nVar; i++) oldRates[i] = _Rates[ NRind[i] ];
      // Compute Residual
      resid = 0;
      for(int i = 0; i < nVar; i++) resid += B[i] * B[i];
      resid = sqrt( resid );
    } // while loop
  }

  /*!
    Tusscher model computes current as pA/pF. We will need to convert it to
    mA/cc and return to body class
  */
  Real Tusscher::Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt,
			     Real volt, Real istim) {
    // Set voltage value in State Array
    _State[0] = volt;
    // iStim is in uA/cc. We need to send pA/pF
    if(userXi) _Algebraic[12] = -istim/(C_m*Xi); 
    else _Algebraic[12] = -istim/(C_m*_Xi); 

    // Update of State Variables
    UpdateStateVariables(dt);

    return _Rates[0];
  }

  //! Get gamma
  Real Tusscher::getGamma() {
    // return 0.9094*(1. - 0.1*tanh( 6.5*log10(_State[3]) + 22.) );
    // return 0.913*(1. - 0.095*tanh( 5.*log10(_State[3]) + 16.5) );
    return 0.9135*(1. - 0.095*tanh(4.*log10(_State[3]) + 13.5) );
  }
}
