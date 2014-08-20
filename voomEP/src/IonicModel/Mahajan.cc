/* UCLA Cell Model */
#include "Mahajan.h"
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

namespace voom {
  Mahajan::Mahajan(Real *Constants, int pos) {
    // Define the Initial State Variables
    _nVar = 26;
    _State.resize( 26 );


/*
_State[0] =-87.1504;_State[1] =0.001079;_State[2] =0.99065;
_State[3] =0.99386;_State[4] =0.007104;_State[5] =0.049736;
_State[6] =0.1061;_State[7] =1.5135;_State[8] =0.22881;
_State[9] =0.26008;_State[10] =104.6294;_State[11] =97.4697;
_State[12] =0.006656;_State[13] =1.8e-005;_State[14] =0.97933;
_State[15] =0.001213;_State[16] =3.4e-005;_State[17] =0.004137;
_State[18] =0.015269;_State[19] =11.4519;_State[20] =0.003648;
_State[21] =0.99308;_State[22] =0.003653;_State[23] =0.17323;
_State[24] =22.3956;_State[25] =20.0081;
*/

switch(pos){
case(0):
{
_State[0] =-81.4359;_State[1] =0.002789;_State[2] =0.96679;
_State[3] =0.97859;_State[4] =0.014994;_State[5] =0.044487;
_State[6] =0.091992;_State[7] =2.1933;_State[8] =0.23277;
_State[9] =0.24868;_State[10] =101.6418;_State[11] =95.4278;
_State[12] =0.009495;_State[13] =3.7e-005;_State[14] =0.96094;
_State[15] =0.003718;_State[16] =6.7e-005;_State[17] =0.010571;
_State[18] =0.024659;_State[19] =11.5909;_State[20] =0.00533;
_State[21] =0.99047;_State[22] =0.005333;_State[23] =0.17488;
_State[24] =21.5955;_State[25] =20.0206;
break;
}
case(1):
{
_State[0] =-81.3846;_State[1] =0.002813;_State[2] =0.96642;
_State[3] =0.97829;_State[4] =0.015109;_State[5] =0.04988;
_State[6] =0.10153;_State[7] =2.3513;_State[8] =0.23616;
_State[9] =0.2529;_State[10] =102.5325;_State[11] =96.0961;
_State[12] =0.010394;_State[13] =3.7e-005;_State[14] =0.95815;
_State[15] =0.004448;_State[16] =6.7e-005;_State[17] =0.01261;
_State[18] =0.024683;_State[19] =11.442;_State[20] =0.005349;
_State[21] =0.99012;_State[22] =0.005352;_State[23] =0.17323;
_State[24] =21.8749;_State[25] =20.2219;
break;
}
case(2):
{
_State[0] =-81.3238;_State[1] =0.002841;_State[2] =0.96596;
_State[3] =0.97791;_State[4] =0.015254;_State[5] =0.056799;
_State[6] =0.11328;_State[7] =2.5207;_State[8] =0.23954;
_State[9] =0.25736;_State[10] =103.4612;_State[11] =96.8014;
_State[12] =0.011385;_State[13] =3.7e-005;_State[14] =0.95449;
_State[15] =0.005413;_State[16] =6.7e-005;_State[17] =0.015284;
_State[18] =0.024706;_State[19] =11.2566;_State[20] =0.00537;
_State[21] =0.98962;_State[22] =0.005373;_State[23] =0.17158;
_State[24] =22.1689;_State[25] =20.4268;
break;
}
case(3):
{
_State[0] =-81.2657;_State[1] =0.002868;_State[2] =0.96553;
_State[3] =0.97752;_State[4] =0.015406;_State[5] =0.064019;
_State[6] =0.12504;_State[7] =2.6684;_State[8] =0.24235;
_State[9] =0.26131;_State[10] =104.2505;_State[11] =97.4054;
_State[12] =0.01227;_State[13] =3.7e-005;_State[14] =0.95055;
_State[15] =0.006458;_State[16] =6.8e-005;_State[17] =0.018158;
_State[18] =0.024723;_State[19] =11.0791;_State[20] =0.005391;
_State[21] =0.98902;_State[22] =0.005394;_State[23] =0.17018;
_State[24] =22.4276;_State[25] =20.6034;
break;
}
case(4):
{
_State[0] =-81.1929;_State[1] =0.002903;_State[2] =0.96497;
_State[3] =0.97694;_State[4] =0.015625;_State[5] =0.073994;
_State[6] =0.14059;_State[7] =2.8407;_State[8] =0.2456;
_State[9] =0.26605;_State[10] =105.1459;_State[11] =98.097;
_State[12] =0.013328;_State[13] =3.7e-005;_State[14] =0.94486;
_State[15] =0.00798;_State[16] =6.8e-005;_State[17] =0.022313;
_State[18] =0.024742;_State[19] =10.8644;_State[20] =0.005417;
_State[21] =0.98806;_State[22] =0.005421;_State[23] =0.16859;
_State[24] =22.7359;_State[25] =20.8156;

break;
}
case(5):
{
_State[0] =-81.1078;_State[1] =0.002944;_State[2] =0.96431;
_State[3] =0.97609;_State[4] =0.01595;_State[5] =0.087401;
_State[6] =0.1603;_State[7] =3.0189;_State[8] =0.24905;
_State[9] =0.27153;_State[10] =106.0478;_State[11] =98.7911;
_State[12] =0.014455;_State[13] =3.8e-005;_State[14] =0.93663;
_State[15] =0.010192;_State[16] =6.8e-005;_State[17] =0.028302;
_State[18] =0.024766;_State[19] =10.6216;_State[20] =0.005448;
_State[21] =0.98648;_State[22] =0.005452;_State[23] =0.16682;
_State[24] =23.0895;_State[25] =21.0607;

break;
}
case(6):
{
_State[0] =-81.2626;_State[1] =0.00287;_State[2] =0.9655;
_State[3] =0.97756;_State[4] =0.015387;_State[5] =0.061634;
_State[6] =0.12278;_State[7] =2.7185;_State[8] =0.24183;
_State[9] =0.2588;_State[10] =104.2029;_State[11] =97.6338;
_State[12] =0.01267;_State[13] =3.7e-005;_State[14] =0.95198;
_State[15] =0.006047;_State[16] =6.8e-005;_State[17] =0.017057;
_State[18] =0.024804;_State[19] =11.0783;_State[20] =0.005392;
_State[21] =0.98941;_State[22] =0.005395;_State[23] =0.16991;
_State[24] =22.2542;_State[25] =20.5438;
break;
}
case(7):
{
_State[0] =-81.1873;_State[1] =0.002905;_State[2] =0.96493;
_State[3] =0.97701;_State[4] =0.015593;_State[5] =0.071141;
_State[6] =0.13799;_State[7] =2.9031;_State[8] =0.24508;
_State[9] =0.26334;_State[10] =105.1225;_State[11] =98.3761;
_State[12] =0.013821;_State[13] =3.7e-005;_State[14] =0.94672;
_State[15] =0.007447;_State[16] =6.8e-005;_State[17] =0.020885;
_State[18] =0.024838;_State[19] =10.8496;_State[20] =0.005419;
_State[21] =0.9886;_State[22] =0.005422;_State[23] =0.16858;
_State[24] =22.5504;_State[25] =20.7522;
break;
}
case(8):
{
_State[0] =-81.1085;_State[1] =0.002943;_State[2] =0.96432;
_State[3] =0.97631;_State[4] =0.015852;_State[5] =0.082505;
_State[6] =0.15518;_State[7] =3.0765;_State[8] =0.2482;
_State[9] =0.26807;_State[10] =105.9653;_State[11] =99.0544;
_State[12] =0.01493;_State[13] =3.8e-005;_State[14] =0.94008;
_State[15] =0.009226;_State[16] =6.8e-005;_State[17] =0.025709;
_State[18] =0.024873;_State[19] =10.6197;_State[20] =0.005448;
_State[21] =0.98743;_State[22] =0.005451;_State[23] =0.16728;
_State[24] =22.8568;_State[25] =20.9674;
break;
}
default:
{
_State[0] =-81.3238;_State[1] =0.002841;_State[2] =0.96596;
_State[3] =0.97791;_State[4] =0.015254;_State[5] =0.056799;
_State[6] =0.11328;_State[7] =2.5207;_State[8] =0.23954;
_State[9] =0.25736;_State[10] =103.4612;_State[11] =96.8014;
_State[12] =0.011385;_State[13] =3.7e-005;_State[14] =0.95449;
_State[15] =0.005413;_State[16] =6.7e-005;_State[17] =0.015284;
_State[18] =0.024706;_State[19] =11.2566;_State[20] =0.00537;
_State[21] =0.98962;_State[22] =0.005373;_State[23] =0.17158;
_State[24] =22.1689;_State[25] =20.4268;

}

}

_Constants = Constants;
  
    _Xi = 1.55E-4/2.58E-8;
  }
  //! Compute the Right hand side of the function
  void Mahajan::ucla_rhsfun(const Real h, const Real istim) {
    Real xik1, xito, xinak, xinacaq, xdif, xiup, xileak;
    Real po, rxa, xicaq, Qr, xina, xikr;
    Real xiks, xinaca, xica;
    
    const Real xkon=0.0327f;
    const Real xkoff=0.0196f;
    const Real btrop=70.0f;
    //conversion factor between micro molar/ms to micro amps/ micro farads
    const Real wca=8.0313f; 
    
    // note: sodium is in m molar so need to divide by 
    const Real xrr=(1.0f/wca)/1000.0f;   
    Real v_old = _State[0]; // For AP onset detection.
    
    // -------- Compute ion channel currents.
    // Sodium current
    xina = comp_ina(h);

    // Inward rectifier K+ current IK1.
    xik1 = comp_ik1();
    // Inward K+ current Ito
    xito = comp_ito(h);
    // Delayed rectifier K+ currents (rapid and slow).
    xikr = comp_ikr(h);
    xiks = comp_iks(h);

    // s
    // --------- Compute pump currents.
    // Na-K pump.
    xinak = comp_inak();
    // Na-Ca exchange ion flux.
    xinacaq = comp_inaca();
    xinaca = wca*xinacaq; // Convert flux to current.
    
    // --------- Compute Ca cycling.
    xdif = (_State[8]-_State[9])/_Constants[14]; //diffusion from submembrane to myoplasm
    // Nonlinear buffering (Troponin).
    _Rates[24] = xkon*_State[9]*(btrop - _State[24]) - xkoff*_State[24];
    _Rates[25] = xkon*_State[8]*(btrop - _State[25]) - xkoff*_State[25];
    
    xiup = comp_iup();
    xileak = comp_ileak();
    
    po = comp_icalpo(h);
    rxa = comp_rxa(); // iCa in the paper.
    
    xicaq = _Constants[4]*po*rxa;// Ca current in micro M/ms
    
    // Derivatives of Ca gating variables.
    _Rates[8] = comp_inst_buffer(_State[8])*		\
      (50.0f*(_State[12]-xdif-xicaq+xinacaq)-_Rates[25]);
    _Rates[9] = comp_inst_buffer(_State[9])*(xdif - xiup + xileak - _Rates[24]);
    _Rates[10] = -_State[12] + xiup - xileak; // SR load dynamics
    _Rates[11] = (_State[10] - _State[11])/_Constants[17]; // NSR-JSR relaxation dynamics
    Qr = comp_Q();
    _Rates[12] = comp_dir(po, Qr, rxa, _Rates[10]); // dJrel/dt
    _Rates[7] = comp_dcp(po, Qr, rxa);

    // Total Ca current.
    xica = 2.0f*wca*xicaq;
  
    // sodium dynamics
    _Rates[19] = -xrr*(xina + 3.0f*xinak + 3.0f*xinaca);
  
   Real ena = (1.0/_Constants[24])*log(_Constants[0]/_State[19]);

 
    // Compute dv, which is simply the sum of the ion currents, since these
    // are given in uA/uF. (C = 1 uF/cm2)
    _Rates[0] = -(xina+xik1+xikr+xiks+xito+xinaca+xinak+xica) + istim;
    //    std::cout << "dVdT: " << _Rates[0] << "\n";
  }  
  
  /* -------------------------------------------------------------
   * Model currents
   */
  //-----------	sodium current following Hund-Rudy -------------------
  Real Mahajan::comp_ina(const Real dt){
    Real ena, a, am, bm, ah, bh, aj, bj;
    Real taum, tauj, tauh;
    Real xina;
    ena = (1.0/_Constants[24])*log(_Constants[0]/_State[19]);
    
    // Newer formulation (Leonid Livschitz).
   

// regular Ina SSA and SSI curves

   a = 1.0 - 1.0/(1.0+exp(-(_State[0]+40.0)/0.24));
    ah = a*0.135*exp((80.0+_State[0])/(-6.8));
    bh = (1.0-a)/(0.13*(1+exp((_State[0]+10.66)/(-11.1)))) +            \
      a*(3.56*exp(0.079*_State[0])+3.1*1.0e5*exp(0.35*_State[0]));
    aj =  a*(-1.2714e5*exp(0.2444*_State[0])-3.474e-5*exp(-0.04391*_State[0]))*(_State[0]+37.78)/(1.0+exp(0.311*(_State[0]+79.23)));
    bj = (1.0-a)*(0.3*exp(-2.535e-7*_State[0])/(1+exp(-0.1*(_State[0]+32)))) + \
      a*(0.1212*exp(-0.01052*_State[0])/(1+exp(-0.1378*(_State[0]+40.14))));

if (fabsf(_State[0]+47.13) < 0.01)
      am = 3.2;
    else
      am = 0.32*(_State[0]+47.13)/(1.0 - exp(-0.1*(_State[0]+47.13)));
    bm = 0.08*exp(-_State[0]/11.0);

    xina= _Constants[15]*_State[2]*_State[3]*_State[1]*_State[1]*_State[1]*(_State[0] - ena);


// SSA and SSI curves shifted left (more negative) by 10 mV to make cell more excitable
/* 
   Real vm,vh,mshift,hshift;
  
   mshift = 0.0;
   hshift = 0.0;


//    mshift = 5.0;
//    hshift = 10.0;



    vm = _State[0] + mshift;
    vh = _State[0] + hshift;
       
    a = 1.0 - 1.0/(1.0+exp(-(vh+40.0)/0.24));
    ah = a*0.135*exp((80.0+vh)/(-6.8));
    bh = (1.0-a)/(0.13*(1+exp((vh+10.66)/(-11.1)))) +		\
      a*(3.56*exp(0.079*vh)+3.1*1.0e5*exp(0.35*vh));
    aj =  a*(-1.2714e5*exp(0.2444*_State[0])-3.474e-5*exp(-0.04391*_State[0]))*(_State[0]+37.78)/(1.0+exp(0.311*(_State[0]+79.23)));
    bj = (1.0-a)*(0.3*exp(-2.535e-7*_State[0])/(1+exp(-0.1*(_State[0]+32)))) + \
      a*(0.1212*exp(-0.01052*_State[0])/(1+exp(-0.1378*(_State[0]+40.14))));
    
    if (fabsf(vm+47.13) < 0.01)
      am = 3.2;
    else
      am = 0.32*(vm+47.13)/(1.0 - exp(-0.1*(vm+47.13)));
    bm = 0.08*exp(-vm/11.0);    
    
    xina= _Constants[15]*_State[2]*_State[3]*_State[1]*_State[1]*_State[1]*(_State[0] - ena);
*/  

  
    // Rush-Larsen method.
    tauh = 1.0/(ah+bh);
    tauj = 1.0/(aj+bj);
    taum = 1.0/(am+bm); 
    _Rates[2] = (ah*tauh - (ah*tauh - _State[2])*exp(-dt/tauh) - _State[2])/dt;
    _Rates[3] = (aj*tauj - (aj*tauj - _State[3])*exp(-dt/tauj) - _State[3])/dt;
    _Rates[1] = (am*taum - (am*taum - _State[1])*exp(-dt/taum) - _State[1])/dt; 
    
    return xina;
  }
  
  //-------------- Ikr following Shannon------------------
  Real Mahajan::comp_ikr(const Real dt)
  {
    const Real ek = (1.0/_Constants[24])*log (_Constants[2]/_Constants[1]);// K reversal potential
  const Real gss = sqrt(_Constants[2]/5.4);

  Real xkrv1, xkrv2, taukr, xkrinf, rg, xikr;
    
    if (fabsf(_State[0] + 7.0) < 0.001/0.123)
      xkrv1=0.00138/0.123;
    else
      xkrv1=0.00138*(_State[0] + 7.0)/(1.0 - exp(-0.123*(_State[0] + 7.0)));
    
    if (fabsf(_State[0] + 10.0) < 0.001/0.145 )
      xkrv2=0.00061/0.145;
    else
      xkrv2=0.00061*(_State[0]+10.0)/(exp(0.145*(_State[0]+10.0)) - 1.0);
    
    taukr = 1.0/(xkrv1 + xkrv2);
    xkrinf = 1.0/(1.0 + exp(-(_State[0] + 50.0)/7.5));
    rg = 1.0/(1.0 + exp((_State[0] + 33.0)/22.4));
    xikr = _Constants[9]*gss*_State[4]*rg*(_State[0] - ek);
    // Rush-Larsen method.
    _Rates[4] = (xkrinf - (xkrinf - _State[4])*exp(-dt/taukr) - _State[4])/dt;
    return xikr;
  }
  
  // ----- Iks modified from Shannon, with new Ca dependence------------
  Real Mahajan::comp_iks(Real dt) {
    const Real prnak=0.018330;
    Real eks, xs1ss, xs2ss, tauxs1, tauxs2, gksx, xiks;
    
    eks = (1.0/_Constants[24])*log((_Constants[2]+prnak*_Constants[0])/(_Constants[1]+prnak*_State[19]));
    xs1ss = 1.0/(1.0 + exp(-(_State[0]-1.50)/16.70));
    xs2ss = xs1ss;
    
    // Adapted from original.
    if (fabsf(_State[0] + 30.0) < 0.001/0.0687)
      tauxs1=1.0/(0.0000719/0.148+0.000131/0.0687);
    else
      tauxs1 = 1.0/(0.0000719*(_State[0] + 30.0)/(1.0 - exp(-0.148*(_State[0] + 30.0))) + 
		    0.000131*(_State[0] + 30.0)/(exp(0.0687*(_State[0] + 30.0)) - 1.0));
    tauxs2 = 4.0*tauxs1;
    gksx = 0.433*(1.0 + 0.8/(1.0 + pow((Real)(0.5/_State[9]),3))); // (q_Ks)
    
    xiks = _Constants[10]*gksx*_State[5]*_State[6]*(_State[0]-eks);
    
    // Rush-Larsen method.
    _Rates[5] = (xs1ss - (xs1ss - _State[5])*exp(-dt/tauxs1) - _State[5])/dt;
    _Rates[6] = (xs2ss - (xs2ss - _State[6])*exp(-dt/tauxs2) - _State[6])/dt;
    
    return xiks;
  }
  
  //------Ik1 following Luo-Rudy formulation (from Shannon model) ------
  Real Mahajan::comp_ik1 ( ){
    Real ek, aki, bki, xkin, xik1;
   // const Real gki = sqrt(_Constants[2]/5.4);
    Real gk1_factor; 
  const Real gki = sqrt(_Constants[2]/5.4);
  
//   ek = (1.0/_Constants[24])*log(_Constants[2]/_Constants[1]); // K reversal potential

 Real P_K = 100.0;
 Real P_Na = 1.0;

// Real P_Na = 0.0;


    ek = (1.0/_Constants[24])*log( (P_K*_Constants[2]+P_Na*_Constants[0])/(P_K*_Constants[1]+P_Na*_State[19])); // K reversal potential

    aki =1.02/ ( 1.0+exp ( 0.2385* ( _State[0]-ek-59.215 ) ) );
    bki = (0.49124*exp(0.08032*(_State[0] - ek + 5.476)) +
	   exp(0.061750*(_State[0] - ek - 594.31)))/(1.0 + exp(-0.5143*(_State[0] - ek + 4.753)));
    xkin = aki/(aki + bki);
//  gk1_factor = 0.55  + 1.25*_State[9];
// xik1 = _Constants[11]*gki*gk1_factor*xkin* (_State[0] - ek);
 

 xik1 = _Constants[11]*gki*xkin* (_State[0] - ek);


    return xik1;
  }
  
  //------- Ito fast following Shannon et. al. 2005 -----------
  //------- Ito slow following Shannon et. al. 2005 -----------
  Real Mahajan::comp_ito(Real dt) {
    // -- Ito,s
    Real ek = 1.0/_Constants[24]*log(_Constants[2]/_Constants[1]);// K reversal potential
    Real rt1 = -(_State[0] + 3.0)/15.0;
    Real rt2 = (_State[0] + 33.5)/10.0;
    Real rt3 = (_State[0] + 60.0)/10.0;
    Real xtos_inf = 1.0/(1.0 + exp(rt1));
    Real ytos_inf = 1.0/(1.0 + exp(rt2));
    Real rs_inf = 1.0/(1.0 + exp(rt2));
    Real txs = 9.0/(1.0+exp(-rt1)) + 0.5;
    Real tys = 3000.0/(1.0+exp(rt3)) + 30.0;
    Real xitos = _Constants[6]*_State[22]*(_State[23] + 0.5*rs_inf)*(_State[0] - ek); // ito slow (original)
    
    // -- Ito,f
    Real xtof_inf = xtos_inf;
    Real ytof_inf = ytos_inf;
    Real txf = 3.5*exp(-(_State[0]/30.0)*(_State[0]/30.0)) + 1.5;
    Real tyf = 20.0/(1.0+exp(rt2)) + 20.0;
    Real xitof = _Constants[7]*_State[20]*_State[21]*(_State[0] - ek); // ito fast
    
    // Rush-Larsen method.
    _Rates[22] = (xtos_inf - (xtos_inf-_State[22])*exp(-dt/txs) -_State[22])/dt;
    _Rates[23] = (ytos_inf - (ytos_inf-_State[23])*exp(-dt/tys) -_State[23])/dt;
    _Rates[20] = (xtof_inf - (xtof_inf-_State[20])*exp(-dt/txf) -_State[20])/dt;
    _Rates[21] = (ytof_inf - (ytof_inf-_State[21])*exp(-dt/tyf) -_State[21])/dt;
    
    return xitos + xitof;
  }
  
  // -------Inak (sodium-potassium exchanger) following Shannon --------------
  Real Mahajan::comp_inak ( ){
    Real fnak, xinak;
    const Real xkmko=1.5;	//these are Inak constants adjusted to fit
    //the experimentally measured dynamic restitution curve
    const Real xkmnai=12.0;
    const Real sigma = ( exp ( _Constants[0]/67.3 )-1.0 ) /7.0;
    
    fnak = 1.0/ ( 1+0.1245*exp ( -0.1*_State[0]*_Constants[24] ) +0.0365*sigma*exp ( -_State[0]*_Constants[24] ) );
    xinak = _Constants[12]*fnak* ( 1./ ( 1.+ ( xkmnai/_State[19] ) ) ) *_Constants[2]/ ( _Constants[2]+xkmko );
    
    return xinak;
  }
  
  // --- Inaca (sodium-calcium exchange) following Shannon and Hund-Rudy------
  //	Note: all concentrations are in mM
  Real Mahajan::comp_inaca ( ) {
    Real zw3, zw4, aloss, yz1, yz2, yz3, yz4, zw8, xinacaq, csm;
    const Real xkdna=0.3; // (cnaca)
    const Real xmcao=1.3;
    const Real xmnao=87.5;
    const Real xmnai=12.3;
    const Real xmcai=0.0036;
    csm = _State[8]/1000.0; // cs is in uMol, but these eqns are in mMol
    
    zw3 = pow(_State[19],3)*_Constants[3]*exp(_State[0]*0.35*_Constants[24]) - pow(_Constants[0],3)*csm*exp(_State[0]*(0.35-1.)*_Constants[24]);
    zw4 = 1.0+0.2*exp ( _State[0]* ( 0.35-1.0 ) *_Constants[24] );
    aloss = 1.0/(1.0+pow((xkdna/_State[8]),3));
    yz1 = xmcao*pow ( _State[19],3 ) +pow ( xmnao,3 ) *csm;
    yz2 = pow ( xmnai,3 ) *_Constants[3]* ( 1.0+csm/xmcai );
    yz3 = xmcai*pow ( _Constants[0],3 ) * ( 1.0+pow ( ( _State[19]/xmnai ),3 ));
    yz4 = pow ( _State[19],3 ) *_Constants[3]+pow ( _Constants[0],3 ) *csm;
    zw8 = yz1+yz2+yz3+yz4;
    
    xinacaq = _Constants[8]*aloss*zw3/ ( zw4*zw8 );
    
    return xinacaq;
  }
  
  
  //	compute driving force (iCa)
  Real Mahajan::comp_rxa() {
    Real za, factor, rxa, csm;
    Real factor1;
    const Real pca = 0.00054;
    csm = _State[8]/1000.0;
    
    za = _State[0]*2.0*_Constants[24];
    factor1 = 4.0*pca*_Constants[23]*_Constants[23]/(_Constants[22]*_Constants[20]); // Temperature is a factor here!
    factor = _State[0]*factor1;
    if (fabsf(za) < 0.001)
      rxa = factor1*(csm*exp(za) - 0.341*(_Constants[3]))/(2.0*_Constants[24]);
    else
      rxa = factor*(csm*exp(za) - 0.341*(_Constants[3]))/(exp(za) - 1.0);
    return rxa;
  }
  
  // ------ Markovian Ca current --------------------------------
  //	Markov model:All parameters have been fitted directly to
  //	experimental current traces using a multidimensional current fitting
  //	routine.
  Real Mahajan::comp_icalpo(Real dt) {
    Real P_o, P_o_inf, P_r, P_s, P_tmp;
    Real alpha, beta, f_c_p, R_V, T_Ca, tau_Ca, tau_Ba;
    Real s1, s2, s2_p, k1;
    Real k3, k3_p, k6, k5, k6_p, k5_p, k4, k4_p;
    
    // Original parameter values from table 4
    const Real cp_tilde = 3.0;   // << Does not appear in Table 4;
    const Real cp_bar = 6.09365; // << Change necessary to reproduce figures (also in CellML)
    const Real tau_po=1.0;
    const Real r1 = 0.3;
    const Real r2 = 3.0;
    const Real s1_p=0.00195;
    const Real k1_p=0.00413;
    const Real k2 = 1.03615e-4;
    const Real k2_p = 0.00224;
    const Real T_Ba = 450.0;
    
    P_o_inf = 1.0/(1.0 + exp(-_State[0]/8.0));
    P_r = 1.0-1.0/(1.0+exp(-(_State[0]+40.0)/4.0)); // Erratum in Mahajan et al.
    P_s = 1.0/(1.0+exp(-(_State[0]+40.0)/11.32));
    P_tmp = exp(-(_State[0]+40.0)/3.0);    // For readability (not in original publication)
    
    alpha = P_o_inf/tau_po;
    // MOD EDL2010MAY03: Avoid potential division by zero below by taking BETA
    // larger than 0.0f  
    beta = max(1.0e-6,(1.0 - P_o_inf)/tau_po);
    
    f_c_p = 1.0/(1.0 + pow(((Real)(cp_tilde/_State[7])),3)); // f(c_p)
    
    R_V = 10.0 + 4954.0*exp(_State[0]/15.6);
    T_Ca = (78.0329 + 0.1*(1.0 + pow(((Real)(_State[7]/cp_bar)),4)))/(1.0 + pow(((Real)(_State[7]/cp_bar)),4));
    
    tau_Ca = (R_V - T_Ca)*P_r + T_Ca;
    tau_Ba = (R_V - T_Ba)*P_r + T_Ba;
    
    s1 = 0.0182688*f_c_p;  // << Change necessary to reproduce figures (also in CellML)
    k1 = 0.024168*f_c_p;   // << idem
    s2 = s1*(k2/k1)*(r1/r2);
    s2_p =  s1_p*(k2_p/k1_p)*(r1/r2);
    
    k3 = P_tmp/(3.0*(1.0+P_tmp));
    k3_p = k3;
    k5 = (1.0f - P_s)/tau_Ca;
    k5_p = (1.0f - P_s)/tau_Ba;
    k6 = f_c_p*P_s/tau_Ca;
    k6_p = P_s/tau_Ba;
    k4 = k3*(alpha/beta)*(k1/k2)*(k5/k6);
    k4_p = k3_p*(alpha/beta)*(k1_p/k2_p)*(k5_p/k6_p);
    
    P_o = 1.0-_State[15]-_State[17]-_State[16]-_State[18]-_State[13]-_State[14];
    
    // Derivatives of Markovian state variables
    _Rates[14] = beta*_State[13] + k5*_State[17] + k5_p*_State[18] - (k6+k6_p+alpha)*_State[14];
    _Rates[13] = alpha*_State[14] + k2*_State[15] + k2_p*_State[16] + r2*P_o - (beta+r1+k1_p+k1)*_State[13];
    _Rates[15] = k1*_State[13] + k4*_State[17] + s1*P_o - (k3+k2+s2)*_State[15];
    _Rates[17] = k3*_State[15] + k6*_State[14] - (k5+k4)*_State[17];
    _Rates[16] = k1_p*_State[13] + k4_p*_State[18] + s1_p*P_o - (k3_p+k2_p+s2_p)*_State[16];
    _Rates[18] = k3_p*_State[16] + k6_p*_State[14] - (k5_p+k4_p)*_State[18];
    
    return P_o;
  }
  
  
  //----- SERCA2a uptake current ------------------------------------
  Real Mahajan::comp_iup() {
    const Real xup = 0.5;// uptake threshold
    
    return _Constants[13]*_State[9]*_State[9]/(_State[9]*_State[9] + xup*xup );
  }
  
  // ---------leak from the SR--------------------------
  Real Mahajan::comp_ileak(){
    const Real gleak = 0.00002069;
    const Real kj = 50.0;
    return gleak * _State[10]*_State[10]/(_State[10]*_State[10]+kj*kj) * (_State[10]*16.667-_State[9]); //vsr/vcell=0.06
  }
  
  // ---------- buffer dynamics in the myoplasm -----------------------
  //buffering to calmodulin and SR are instantaneous, while buffering to
  //Troponin C is time dependent.These are important to have reasonable
  //Ca transient.Note: we have buffering in the suB_membrane space and
  //the myoplasm.
  Real Mahajan::comp_inst_buffer(Real c) {
    Real SR_buf, Cal_buf, Mem_buf, Sarco_buf, ATP_buf;
    const Real B_Cd = 24.0;
    const Real K_Cd = 7.0;
    const Real B_SR = 47.0;
    const Real K_SR = 0.6;
    const Real B_mem = 15.0;
    const Real K_mem = 0.3;
    const Real B_sar = 42.0;
    const Real K_sar = 13.0;
    
    SR_buf = B_SR*K_SR/((K_SR+c)*(K_SR+c)); // SR binding sites buffering.
    Cal_buf = B_Cd*K_Cd/((K_Cd+c)*(K_Cd+c));  // Calmodulin Ca buffering.
    Mem_buf = B_mem*K_mem/((K_mem+c)*(K_mem+c));  // Membrane Ca buffering.
    Sarco_buf = B_sar*K_sar/((K_sar+c)*(K_sar+c));  // SR Ca buffering
    
    return 1.0/(1.0 + SR_buf + Cal_buf + Mem_buf + Sarco_buf);
  }
  
  // --------- release-load functional dependence ----------------
  Real Mahajan::comp_Q() {
    Real bv;
    
    bv = (_Constants[19]-50.0) - _Constants[18]*_Constants[19];
    if (_State[11] < 50)
      return 0.0;
    else if ((_State[11] > 50.0) && (_State[11] < _Constants[19]))
      return _State[11] - 50.0;
    else
      return _Constants[18]*_State[11] + bv;  // bv is called s in the original publication.
  }
  
  Real Mahajan::comp_dir(Real po, Real Qr, Real rxa, Real dcj) {
    Real spark_rate, gRyRV;
    const Real ay = 0.05;
    const Real gRyR = 2.58079f;
    
    gRyRV = gRyR*exp(-ay*(_State[0]+30.0))/(1.0+exp(-ay*(_State[0]+30.0)));
    spark_rate=_Constants[5]*po*fabsf(rxa)*gRyRV; // minus and without ABS?
    
    return spark_rate*_State[10]*Qr/_Constants[19] - _State[12]*(1.0-_Constants[16]*_Rates[10]/_State[10])/_Constants[16];
  }
  
  // ----------- dyadic junction dynamics ------------------------
  Real Mahajan::comp_dcp(Real po, Real Qr, Real rxa) {
    Real gSRV, JSRtld, JCatld;
    const Real gSRbar=26841.8f; // m mol/(cm C)
    const Real ax = 0.3576f; 
    const Real gCabar = 9000.0f;; // m mol/(cm C)
    const Real taups = 0.5;
    
    gSRV  = gSRbar*exp(-ax*(_State[0]+30.0))/(1.0+exp(-ax*(_State[0]+30.0)));
    JSRtld = gSRV*Qr*po*fabsf(rxa); // minus and without fabsf?
    JCatld = _Constants[5]*gCabar*po*fabsf(rxa); // minus and without fabsf?
    
    return JSRtld + JCatld - (_State[7]-_State[8])/taups;
  }
  
  
  // Runge Kutta Update
  void Mahajan::UpdateStateVariables(const Real dt, const Real i_stim) {
    const Real MAX_DV_DT = 20.0;
    const int EULER_ADAPT_FCT = 10;
    ucla_rhsfun(dt, i_stim);
    if ( abs(_Rates[0]) > MAX_DV_DT) {
      for(int i = 0; i < EULER_ADAPT_FCT; i++) {
	ucla_rhsfun( dt/EULER_ADAPT_FCT, i_stim);
	for(int j = 0; j < _nVar; j++)
	  _State[j] += _Rates[j] * dt/EULER_ADAPT_FCT;
      }
    } else {
      for(int j = 0; j < _nVar; j++)
	_State[j] += _Rates[j] * dt;
    }
  }
  
  // Computing Ionic Current
  Real Mahajan::Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt, 
				  Real Volt, Real istim) {
      // Set Voltage Value
      _State[0] = Volt;







      // Stimulus current needs to be in uA/uF
      Real i_stim;
      i_stim = (userXi) ? istim/(C_m *Xi) : istim/(C_m*_Xi);

      // Compute Data for one time step
      UpdateStateVariables(dt, i_stim);

      //Return dV/dT
      return _Rates[0];
    }
}
