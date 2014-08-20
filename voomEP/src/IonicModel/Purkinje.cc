#include "Purkinje.h"

namespace voom {
  //! Compute the Rates
  void Purkinje::ComputeRates() {
    _Algebraic[1] = 1.00000/(1.00000+(exp(((_State[0]+47.8000)/-5.50000))));
    _Rates[6] = (_Algebraic[1] - _State[6])/_Constants[16];
    _Algebraic[3] = 1.00000/(1.00000+(exp(((_State[0]+14.6000)/-5.50000))));
    _Rates[8] = (_Algebraic[3] - _State[8])/_Constants[18];
    _Algebraic[6] = 1.00000/(1.00000+(exp(((_State[0]+-7.00000)/-9.00000))));
    _Rates[11] = (_Algebraic[6] - _State[11])/_Constants[20];
    _Algebraic[7] = 1.00000/(1.00000+(exp(((_State[0]+27.5000)/8.00000))));
    _Rates[12] = (_Algebraic[7] - _State[12])/_Constants[21];
    _Algebraic[8] = 1.00000/(1.00000+(exp(((_State[0]+25.0000)/-5.00000))));
    _Rates[13] = (_Algebraic[8] - _State[13])/_Constants[24];
    _Algebraic[9] = 1.00000/(1.00000+(exp(((_State[0]+69.0000)/3.96000))));
    _Rates[14] = (_Algebraic[9] - _State[14])/_Constants[25];
    _Algebraic[10] = 1.00000/(1.00000+(exp(((_State[0]+30.0000)/-5.00000))));
    _Rates[15] = (_Algebraic[10] - _State[15])/_Constants[27];
    _Algebraic[2] = 1.00000/(1.00000+(exp(((_State[0]+67.9000)/3.87000))));
    _Algebraic[17] =  1.42271*(exp(( -0.0511900*_State[0])));
    _Rates[7] = (_Algebraic[2] - _State[7])/_Algebraic[17];
    _Algebraic[4] = 1.00000/(1.00000+(exp(((_State[0]+31.0000)/5.54000))));
    _Algebraic[18] = 25.1000/(0.0400000+ 0.700000*(exp(( -1.00000*(pow(( 0.0280000*(_State[0]+14.5000)), 2.00000))))));
    _Rates[9] = (_Algebraic[4] - _State[9])/_Algebraic[18];
    _Algebraic[5] = 0.400000+0.600000/(1.00000+(pow((_State[1]/0.000100000), 2.00000)));
    _Algebraic[19] = 2.00000+80.0000/(1.00000+(pow((_State[1]/0.000100000), 2.00000)));
    _Rates[10] = (_Algebraic[5] - _State[10])/_Algebraic[19];
    _Algebraic[11] = 0.100000+0.900000/(1.00000+(exp(((_State[0]+75.6000)/6.30000))));
    _Algebraic[20] = 120.000+ 1.00000*(exp(((_State[0]+100.000)/25.0000)));
    _Rates[16] = (_Algebraic[11] - _State[16])/_Algebraic[20];
    _Algebraic[13] = 1.00000/(1.00000+(exp(((_State[0] - 1.50000)/-16.7000))));
    _Algebraic[22] = ((fabs((_State[0]+30.0000)))<0.0145000 ? 417.946 : 1.00000/(( 7.19000e-05*(_State[0]+30.0000))/(1.00000 - (exp(( -0.148000*(_State[0]+30.0000)))))+( 0.000131000*(_State[0]+30.0000))/((exp(( 0.0687000*(_State[0]+30.0000)))) - 1.00000)));
    _Rates[18] = (_Algebraic[13] - _State[18])/_Algebraic[22];
    _Algebraic[14] = 1.00000/(1.00000+(exp(((_State[0]+109.000)/10.0000))));
    _Algebraic[23] = 6000.00/((exp(( -1.00000*(2.90000+ 0.0400000*_State[0]))))+(exp(( 1.00000*(3.60000+ 0.110000*_State[0])))));
    _Rates[20] = (_Algebraic[14] - _State[20])/_Algebraic[23];
    _Algebraic[15] = 1.00000/(1.00000+(exp(((_State[0]+109.000)/10.0000))));
    _Algebraic[24] = 6000.00/((exp(( -1.00000*(2.90000+ 0.0400000*_State[0]))))+(exp(( 1.00000*(3.60000+ 0.110000*_State[0])))));
    _Rates[21] = (_Algebraic[15] - _State[21])/_Algebraic[24];
    _Algebraic[12] = 1.00000/(1.00000+(exp(((_State[0]+50.0000)/-7.50000))));
    _Algebraic[21] = ((fabs((_State[0]+7.00000)))>0.00100000 ? ( 0.00138000*1.00000*(_State[0]+7.00000))/(1.00000 - (exp(( -0.123000*(_State[0]+7.00000))))) : 0.00138000/0.123000);
    _Algebraic[26] = ((fabs((_State[0]+10.0000)))>0.00100000 ? ( 6.10000e-05*1.00000*(_State[0]+10.0000))/((exp(( 0.145000*(_State[0]+10.0000)))) - 1.00000) : 0.000610000/0.145000);
    _Algebraic[29] = 1.00000/(_Algebraic[21]+_Algebraic[26]);
    _Rates[17] = (_Algebraic[12] - _State[17])/_Algebraic[29];
    _Algebraic[27] = _Algebraic[13];
    _Algebraic[30] =  4.00000*_Algebraic[22];
    _Rates[19] = (_Algebraic[27] - _State[19])/_Algebraic[30];
    _Algebraic[41] = 1.00000/(1.00000+(exp(((92.0000+_State[0])/10.0000))));
    _Algebraic[42] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[43] =  _Constants[29]*(pow((_Constants[5]/5.40000), 0.800000))*_Algebraic[41]*(_State[0] - _Algebraic[42]);
    _Algebraic[32] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[33] =  _Constants[22]*_State[11]*_State[12]*(_State[0] - _Algebraic[32]);
    _Algebraic[34] = 1.00000/(1.00000+(exp(((5.00000 - _State[0])/17.0000))));
    _Algebraic[35] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[36] =  _Constants[23]*_Algebraic[34]*(_State[0] - _Algebraic[35]);
    _Algebraic[44] = 1.00000/(1.00000+(exp(((33.0000+_State[0])/22.4000))));
    _Algebraic[45] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[46] =  _Constants[30]*(pow((_Constants[5]/5.40000), 1.00000))*_Algebraic[44]*_State[17]*(_State[0] - _Algebraic[45]);
    _Algebraic[47] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[48] =  _Constants[31]*_State[18]*_State[19]*(_State[0] - _Algebraic[47]);
    _Algebraic[53] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[54] =  _Constants[37]*(_State[0] - _Algebraic[53]);
    _Algebraic[50] = 1.00000/(1.00000+(exp(((_State[0]+80.0000)/-45.0000))));
    _Algebraic[51] = 1.00000/(1.00000+(exp(((_State[0]+0.00000)/125.000))));
    _Algebraic[52] =  _Constants[36]*_Algebraic[50]*_Algebraic[51]*(1.00000/(1.00000+(pow((1.90000/_Constants[5]), 1.45000))))*(1.00000/(1.00000+(pow((31.9800/_State[4]), 1.00000))));
    _Algebraic[60] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[62] =  _Constants[43]*_State[20]*(_State[0] - _Algebraic[60]);

    _Rates[5] = ( ( -1.00000*_Algebraic[33]+ -1.00000*_Algebraic[36]+ -1.00000*_Algebraic[46]+ -1.00000*_Algebraic[48]+ -1.00000*_Algebraic[43]+ -1.00000*_Algebraic[54]+ -1.00000*_Algebraic[62]+ -1.00000*_Algebraic[0]+ 2.00000*_Algebraic[52])*1000.00)/( _Constants[1]*_Constants[61]);
    _Algebraic[28] =  0.500000*_Constants[58]*(log((_Constants[3]/_State[1])));
    _Algebraic[31] =  _Constants[19]*_State[8]*_State[9]*_State[10]*(_State[0] - _Algebraic[28]);
    _Algebraic[16] =  0.500000*_Constants[58]*(log((_Constants[3]/_State[1])));
    _Algebraic[25] =  _Constants[17]*_State[6]*_State[7]*(_State[0] - _Algebraic[16]);
    _Algebraic[37] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[38] =  _Constants[26]*_State[13]*_State[14]*(_State[0] - _Algebraic[37]);
    _Algebraic[39] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[40] =  _Constants[28]*_State[15]*_State[16]*(_State[0] - _Algebraic[39]);
    _Algebraic[59] =  _Constants[40]*(1.00000/(1.00000+(pow((_Constants[41]/_State[1]), _Constants[42]))));
    _Algebraic[57] =  0.500000*_Constants[58]*(log((_Constants[3]/_State[1])));
    _Algebraic[58] =  _Constants[39]*(_State[0] - _Algebraic[57]);
    _Algebraic[64] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[66] =  _Constants[44]*_State[21]*(_State[0] - _Algebraic[64]);
    _Algebraic[49] = ( 512.009*_Constants[33]*( (exp((( _Constants[35]*(_Constants[32] - 2.00000)*_State[0])/_Constants[58])))*(pow(_State[4], _Constants[32]))*_Constants[3] -  _State[1]*(exp((( (_Constants[35] - 1.00000)*(_Constants[32] - 2.00000)*_State[0])/_Constants[58])))*(pow(_Constants[4], _Constants[32]))))/( (1.00000+ _Constants[34]*( _State[1]*(pow(_Constants[4], _Constants[32]))+ _Constants[3]*(pow(_State[4], _Constants[32]))))*(1.00000+_State[1]/0.00690000));
    _Algebraic[55] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[56] =  _Constants[38]*(_State[0] - _Algebraic[55]);
    // This is dV/dt
    _Rates[0] =  (( - 1.00000*1.00000)/_Constants[6])*
      (_Algebraic[43]+ _Algebraic[33]+_Algebraic[36]+_Algebraic[46]+
       _Algebraic[48]+_Algebraic[54]+_Algebraic[52]+_Algebraic[31]+
       _Algebraic[25]+_Algebraic[38]+_Algebraic[40]+_Algebraic[59]+
       _Algebraic[58]+_Algebraic[66]+_Algebraic[62]+_Algebraic[49]+
       _Algebraic[56]+_Algebraic[0]);
    _Algebraic[61] =  _Constants[45]*(_State[3]/(1.00000+(pow((_Constants[46]/_State[1]), 2.00000))));
    _Algebraic[63] =  _Constants[47]*(((pow((_State[1]/_Constants[48]), _Constants[50])) - (pow((_State[3]/_Constants[49]), _Constants[50])))/(1.00000+(pow((_State[1]/_Constants[48]), _Constants[50]))+(pow((_State[3]/_Constants[49]), _Constants[50]))));
    _Algebraic[65] =  _Constants[51]*(((pow((_State[2]/_Constants[52]), _Constants[54])) - (pow((_State[3]/_Constants[53]), _Constants[54])))/(1.00000+(pow((_State[2]/_Constants[52]), _Constants[54]))+(pow((_State[3]/_Constants[53]), _Constants[54]))));
    _Algebraic[67] =  _Constants[55]*(_State[3] - _State[2]);
    _Rates[3] = ( ( -1.00000*_Algebraic[61]+ -1.00000*_Algebraic[67]+_Algebraic[65]+_Algebraic[63])*1000.00)/( 2.00000*_Constants[1]*_Constants[62]);
    _Rates[4] = ( ( -1.00000*_Algebraic[38]+ -1.00000*_Algebraic[40]+ -3.00000*_Algebraic[52]+ -3.00000*_Algebraic[49]+ -1.00000*_Algebraic[66]+ -1.00000*_Algebraic[56])*1000.00)/( _Constants[1]*_Constants[61]);
    _Algebraic[68] =  _Constants[56]*(_State[1] - _State[2]);

    _Rates[1] = ( ( -1.00000*_Algebraic[31]+ -1.00000*_Algebraic[25]+ -1.00000*_Algebraic[58]+ 1.00000*_Algebraic[61]+ -1.00000*_Algebraic[63]+ -1.00000*_Algebraic[68]+ 2.00000*_Algebraic[49]+ -1.00000*_Algebraic[59])*1000.00)/( 2.00000*_Constants[1]*_Constants[59]);
    _Rates[2] = ( ( -1.00000*_Algebraic[65]+_Algebraic[68]+_Algebraic[67])*1000.00)/( 2.00000*_Constants[1]*_Constants[60]);
  }
  
  void Purkinje::ComputeVariables() {
    _Algebraic[1] = 1.00000/(1.00000+(exp(((_State[0]+47.8000)/-5.50000))));
    _Algebraic[3] = 1.00000/(1.00000+(exp(((_State[0]+14.6000)/-5.50000))));
    _Algebraic[6] = 1.00000/(1.00000+(exp(((_State[0]+-7.00000)/-9.00000))));
    _Algebraic[7] = 1.00000/(1.00000+(exp(((_State[0]+27.5000)/8.00000))));
    _Algebraic[8] = 1.00000/(1.00000+(exp(((_State[0]+25.0000)/-5.00000))));
    _Algebraic[9] = 1.00000/(1.00000+(exp(((_State[0]+69.0000)/3.96000))));
    _Algebraic[10] = 1.00000/(1.00000+(exp(((_State[0]+30.0000)/-5.00000))));
    _Algebraic[2] = 1.00000/(1.00000+(exp(((_State[0]+67.9000)/3.87000))));
    _Algebraic[17] =  1.42271*(exp(( -0.0511900*_State[0])));
    _Algebraic[4] = 1.00000/(1.00000+(exp(((_State[0]+31.0000)/5.54000))));
    _Algebraic[18] = 25.1000/(0.0400000+ 0.700000*(exp(( -1.00000*(pow(( 0.0280000*(_State[0]+14.5000)), 2.00000))))));
    _Algebraic[5] = 0.400000+0.600000/(1.00000+(pow((_State[1]/0.000100000), 2.00000)));
    _Algebraic[19] = 2.00000+80.0000/(1.00000+(pow((_State[1]/0.000100000), 2.00000)));
    _Algebraic[11] = 0.100000+0.900000/(1.00000+(exp(((_State[0]+75.6000)/6.30000))));
    _Algebraic[20] = 120.000+ 1.00000*(exp(((_State[0]+100.000)/25.0000)));
    _Algebraic[13] = 1.00000/(1.00000+(exp(((_State[0] - 1.50000)/-16.7000))));
    _Algebraic[22] = ((fabs((_State[0]+30.0000)))<0.0145000 ? 417.946 : 1.00000/(( 7.19000e-05*(_State[0]+30.0000))/(1.00000 - (exp(( -0.148000*(_State[0]+30.0000)))))+( 0.000131000*(_State[0]+30.0000))/((exp(( 0.0687000*(_State[0]+30.0000)))) - 1.00000)));
    _Algebraic[14] = 1.00000/(1.00000+(exp(((_State[0]+109.000)/10.0000))));
    _Algebraic[23] = 6000.00/((exp(( -1.00000*(2.90000+ 0.0400000*_State[0]))))+(exp(( 1.00000*(3.60000+ 0.110000*_State[0])))));
    _Algebraic[15] = 1.00000/(1.00000+(exp(((_State[0]+109.000)/10.0000))));
    _Algebraic[24] = 6000.00/((exp(( -1.00000*(2.90000+ 0.0400000*_State[0]))))+(exp(( 1.00000*(3.60000+ 0.110000*_State[0])))));
    _Algebraic[12] = 1.00000/(1.00000+(exp(((_State[0]+50.0000)/-7.50000))));
    _Algebraic[21] = ((fabs((_State[0]+7.00000)))>0.00100000 ? ( 0.00138000*1.00000*(_State[0]+7.00000))/(1.00000 - (exp(( -0.123000*(_State[0]+7.00000))))) : 0.00138000/0.123000);
    _Algebraic[26] = ((fabs((_State[0]+10.0000)))>0.00100000 ? ( 6.10000e-05*1.00000*(_State[0]+10.0000))/((exp(( 0.145000*(_State[0]+10.0000)))) - 1.00000) : 0.000610000/0.145000);
    _Algebraic[29] = 1.00000/(_Algebraic[21]+_Algebraic[26]);
    _Algebraic[27] = _Algebraic[13];
    _Algebraic[30] =  4.00000*_Algebraic[22];
    _Algebraic[41] = 1.00000/(1.00000+(exp(((92.0000+_State[0])/10.0000))));
    _Algebraic[42] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[43] =  _Constants[29]*(pow((_Constants[5]/5.40000), 0.800000))*_Algebraic[41]*(_State[0] - _Algebraic[42]);
    _Algebraic[32] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[33] =  _Constants[22]*_State[11]*_State[12]*(_State[0] - _Algebraic[32]);
    _Algebraic[34] = 1.00000/(1.00000+(exp(((5.00000 - _State[0])/17.0000))));
    _Algebraic[35] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[36] =  _Constants[23]*_Algebraic[34]*(_State[0] - _Algebraic[35]);
    _Algebraic[44] = 1.00000/(1.00000+(exp(((33.0000+_State[0])/22.4000))));
    _Algebraic[45] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[46] =  _Constants[30]*(pow((_Constants[5]/5.40000), 1.00000))*_Algebraic[44]*_State[17]*(_State[0] - _Algebraic[45]);
    _Algebraic[47] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[48] =  _Constants[31]*_State[18]*_State[19]*(_State[0] - _Algebraic[47]);
    _Algebraic[53] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[54] =  _Constants[37]*(_State[0] - _Algebraic[53]);
    _Algebraic[50] = 1.00000/(1.00000+(exp(((_State[0]+80.0000)/-45.0000))));
    _Algebraic[51] = 1.00000/(1.00000+(exp(((_State[0]+0.00000)/125.000))));
    _Algebraic[52] =  _Constants[36]*_Algebraic[50]*_Algebraic[51]*(1.00000/(1.00000+(pow((1.90000/_Constants[5]), 1.45000))))*(1.00000/(1.00000+(pow((31.9800/_State[4]), 1.00000))));
    _Algebraic[60] =  _Constants[58]*(log((_Constants[5]/_State[5])));
    _Algebraic[62] =  _Constants[43]*_State[20]*(_State[0] - _Algebraic[60]);
    _Algebraic[28] =  0.500000*_Constants[58]*(log((_Constants[3]/_State[1])));
    _Algebraic[31] =  _Constants[19]*_State[8]*_State[9]*_State[10]*(_State[0] - _Algebraic[28]);
    _Algebraic[16] =  0.500000*_Constants[58]*(log((_Constants[3]/_State[1])));
    _Algebraic[25] =  _Constants[17]*_State[6]*_State[7]*(_State[0] - _Algebraic[16]);
    _Algebraic[37] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[38] =  _Constants[26]*_State[13]*_State[14]*(_State[0] - _Algebraic[37]);
    _Algebraic[39] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[40] =  _Constants[28]*_State[15]*_State[16]*(_State[0] - _Algebraic[39]);
    _Algebraic[59] =  _Constants[40]*(1.00000/(1.00000+(pow((_Constants[41]/_State[1]), _Constants[42]))));
    _Algebraic[57] =  0.500000*_Constants[58]*(log((_Constants[3]/_State[1])));
    _Algebraic[58] =  _Constants[39]*(_State[0] - _Algebraic[57]);
    _Algebraic[64] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[66] =  _Constants[44]*_State[21]*(_State[0] - _Algebraic[64]);
    _Algebraic[49] = ( 512.009*_Constants[33]*( (exp((( _Constants[35]*(_Constants[32] - 2.00000)*_State[0])/_Constants[58])))*(pow(_State[4], _Constants[32]))*_Constants[3] -  _State[1]*(exp((( (_Constants[35] - 1.00000)*(_Constants[32] - 2.00000)*_State[0])/_Constants[58])))*(pow(_Constants[4], _Constants[32]))))/( (1.00000+ _Constants[34]*( _State[1]*(pow(_Constants[4], _Constants[32]))+ _Constants[3]*(pow(_State[4], _Constants[32]))))*(1.00000+_State[1]/0.00690000));
    _Algebraic[55] =  _Constants[58]*(log((_Constants[4]/_State[4])));
    _Algebraic[56] =  _Constants[38]*(_State[0] - _Algebraic[55]);
    _Algebraic[61] =  _Constants[45]*(_State[3]/(1.00000+(pow((_Constants[46]/_State[1]), 2.00000))));
    _Algebraic[63] =  _Constants[47]*(((pow((_State[1]/_Constants[48]), _Constants[50])) - (pow((_State[3]/_Constants[49]), _Constants[50])))/(1.00000+(pow((_State[1]/_Constants[48]), _Constants[50]))+(pow((_State[3]/_Constants[49]), _Constants[50]))));
    _Algebraic[65] =  _Constants[51]*(((pow((_State[2]/_Constants[52]), _Constants[54])) - (pow((_State[3]/_Constants[53]), _Constants[54])))/(1.00000+(pow((_State[2]/_Constants[52]), _Constants[54]))+(pow((_State[3]/_Constants[53]), _Constants[54]))));
    _Algebraic[67] =  _Constants[55]*(_State[3] - _State[2]);
    _Algebraic[68] =  _Constants[56]*(_State[1] - _State[2]);
  }

  // Forward Euler Adaptive Update of state variables
  // Newton Raphson update of State Variables
  void Purkinje::UpdateStateVariables(const Real dt) {
    // TIme adaptive Forward Euler
    const Real MAX_DV_DT = 1.5;
    const Real EULER_ADAPT_FCT = 10;
    ComputeRates();
    if ( fabs(_Rates[0]) > MAX_DV_DT ){
      for(int i = 0; i < EULER_ADAPT_FCT; i++) {
	ComputeRates();
	for(int j = 0; j < 22; j++) _State[j] += _Rates[j]*dt/EULER_ADAPT_FCT;
      }
    }else {
      for(int j = 0; j < 22; j++) _State[j] += _Rates[j]*dt;
    }
  }

  /*!
    Purkinje model computes current as uA/uF. 
  */
  Real Purkinje::Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, 
			     Real volt, Real istim){

    // Set voltage value in State array
    _State[0]     = volt;
    const Real conversion = -1E6; //uA to pA
    // iStim should go as pA
    if (userXi) _Algebraic[0] = conversion * istim * _SurfaceArea/Xi;
    else _Algebraic[0] = conversion * istim * _SurfaceArea/_Xi;

    // Update State Variables
    UpdateStateVariables(dt);
    return _Rates[0];
  }

  //! Get Gamma active contraction values
  Real Purkinje::getGamma() {
//    return 0.2127*(4. - 0.7*tanh(2.*log(_State[1]) + 3.));
    return 1.;
  }
  
}
