/* UCLA Cell Model */
#include "Mahajan_fail.h"
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

namespace voom {
  Mahajan_fail::Mahajan_fail(Real *Constants, int pos) {
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
void Mahajan_fail::ucla_rhsfun(double h, double istim) {
    double xik1, xito, xinak, jnaca, jd, jup, jleak;
    double po, rxa, jca, Qr, xina, xikr;
    double xiks, xinaca, xica;

    const double xkon=0.0327;
    const double xkoff=0.0196;
    const double btrop=70.0;

    xik1=comp_ik1();
    xito=comp_ito(h);//itos and itof
    xinak=comp_inak();

    jnaca = comp_inaca();
//----------- Equations for Ca cycling -------------------------
    jd = _Constants[28] * (_State[8]-_State[9])/_Constants[14];//diffusion from submembrane to myoplasm
    // Troponin kinetics

    _Rates[24] = _Constants[29] * (xkon*_State[9]*(btrop-_State[24])-xkoff*_State[24]);
    _Rates[25] = xkon*_State[8]*(btrop-_State[25])-xkoff*_State[25];

    jup = _Constants[30] * comp_iup();
    jleak = _Constants[27] * comp_ileak();

    po = comp_icalpo(h);
    rxa = comp_rxa();
    jca = _Constants[4]*po*rxa;// Ca current in micro M/ms

    _Rates[8] = comp_inst_buffer(_State[8])*(50.0*(_State[12]-jd-jca+jnaca)-_Rates[25]);
    _Rates[9] = comp_inst_buffer(_State[9])*(jd-jup+jleak-_Rates[24]);
    _Rates[10] = -_State[12]+jup-jleak;// SR load dynamics
    _Rates[11]=(_State[10]-_State[11])/_Constants[17];// NSR-JSR relaxation dynamics
    Qr=comp_Q();
    _Rates[12]=comp_dir(po, Qr, rxa, _Rates[10]);
    _Rates[7]=comp_dcp(po, Qr, rxa);

    xina=comp_ina(h);
    xikr=comp_ikr(h);
    xiks=comp_iks(h);

    // std::cout << _Rates[7] << " " << _Rates[8] << " " << _Rates[9] << " " << _Rates[10] << " " << _Rates[11] << " " << Qr << " " << _Rates[12] << endl;

//-------convert ion flow to current---------
    const double wca=8.0313;//conversion factor between micro molar/ms to micro amps/ micro farads
    xinaca=wca*jnaca;
    xica=2.0*wca*jca;
//--------sodium dynamics -------------------------
    const double xrr=(1.0/wca)/1000.0;// note: sodium is in m molar so need to divide by 1000
    _Rates[19] = (-xrr*(xina+3.0*xinak+3.0*xinaca));
// --------	dV/dt ------------------------------------
    _Rates[0] = (-(xina+xik1+xikr+xiks+xito+xinaca+xica+xinak)+ istim);
} // End of compute rhs function



/* -------------------------------------------------------------
 * Model currents
 */
//-----------	sodium current following Hund-Rudy -------------------
double Mahajan_fail::comp_ina(double dt){
    double ena, a, am, bm, ah, bh, aj, bj;
    double taum, tauj, tauh;
    double xina;

    ena = (1.0/_Constants[24])*log(_Constants[0]/_State[19]);
    am=3.2;

    if (fabs(_State[0]+47.13)>0.001) am = 0.32*(_State[0]+47.13)/(1.0-exp(-0.1*(_State[0]+47.13)));
    bm = 0.08*exp(-_State[0]/11.0);


    // Newer formulation (Leonid Livschitz).
    // regular Ina SSA and SSI curves
    a = 1.0 - 1.0/(1.0+exp(-(_State[0]+40.0)/0.24));
    ah = a*0.135*exp((80.0+_State[0])/(-6.8));
    bh = (1.0-a)/(0.13*(1+exp((_State[0]+10.66)/(-11.1)))) +		\
      a*(3.56*exp(0.079*_State[0])+3.1*1.0e5*exp(0.35*_State[0]));
    aj =  a*(-1.2714e5*exp(0.2444*_State[0])-3.474e-5*exp(-0.04391*_State[0]))*(_State[0]+37.78)/(1.0+exp(0.311*(_State[0]+79.23)));
    bj = (1.0-a)*(0.3*exp(-2.535e-7*_State[0])/(1+exp(-0.1*(_State[0]+32)))) + \
    a*(0.1212*exp(-0.01052*_State[0])/(1+exp(-0.1378*(_State[0]+40.14))));

    // Zhilin's Code:
    /*
    if(_State[0]<(-40.0))
    {
             ah=0.135*exp((80.0+_State[0])/(-6.8));
             bh=3.56*exp(0.079*_State[0])+310000.0*exp(0.35*_State[0]);
             aj=((-127140.0*exp(0.2444*_State[0])-0.00003474*exp(-0.04391*_State[0]))*(_State[0]+37.78))/(1.0+exp(0.311*(_State[0]+79.23)));
             bj=(0.1212*exp(-0.01052*_State[0]))/(1.0+exp(-0.1378*(_State[0]+40.14)));
    }
    else
    {
             ah=0.0;
             bh=1.0/(0.13*(1.0+exp((_State[0]+10.66)/(-11.1))));
             aj=0.0;
             bj=(0.3*exp(-0.0000002535*_State[0]))/(1.0+exp(-0.1*(_State[0]+32.0)));
    }
    */

    tauh = 1.0/(ah+bh);
    tauj = _Constants[34] * 1.0/(aj+bj);
    taum = 1.0/(am+bm);
    xina= _Constants[15]*(_State[2]+_Constants[32])/(1 + _Constants[32])*(_State[3]+_Constants[33])/(1 + _Constants[33])*(_State[1]+_Constants[31])/(1+_Constants[31])*(_State[1]+_Constants[31])/(1+_Constants[31])*(_State[1]+_Constants[31])/(1+_Constants[31])*(_State[0]-ena);

    _Rates[2] = (ah*tauh - (ah*tauh - _State[2])*exp(-dt/tauh) - _State[2])/dt;
    _Rates[3] = (aj*tauj - (aj*tauj - _State[3])*exp(-dt/tauj) - _State[3])/dt;
    _Rates[1] = (am*taum - (am*taum - _State[1])*exp(-dt/taum) - _State[1])/dt;

    return xina;

} // comp ina



//-------------- Ikr following Shannon------------------
double Mahajan_fail::comp_ikr(double dt)
{
    const double ek = (1.0/_Constants[24])*log(_Constants[2]/_Constants[1]);// K reversal potential
    const double gss = sqrt(_Constants[2]/5.4);

    double xkrv1 = 0.00138/0.123;
    if (fabs(_State[0]+7.0)>0.001)
        xkrv1=0.00138*(_State[0]+7.0)/( 1.-exp(-0.123*(_State[0]+7.0)));
    double xkrv2=0.00061/0.145;
    if (fabs(_State[0]+10.0)>0.001)
        xkrv2 = 0.00061*(_State[0]+10.0)/(exp( 0.145*(_State[0]+10.0))-1.0);

    double taukr = 1.0/(xkrv1+xkrv2);
    double xkrinf = 1.0/(1.0+exp(-(_State[0]+50.0)/7.5));
    double rg = 1.0/(1.0+exp((_State[0]+33.0)/22.4));
    double xikr=_Constants[9]*gss*_State[4]*rg*(_State[0]-ek);

    _Rates[4] = (xkrinf - (xkrinf - _State[4])*exp(-dt/taukr) - _State[4])/dt;
    return xikr;
} //  comp ikr



// ----- Iks modified from Shannon, with new Ca dependence------------
double Mahajan_fail::comp_iks(double dt) {
    const double prnak=0.018330;

    double eks = (1.0/_Constants[24])*log((_Constants[2]+prnak*_Constants[0])/(_Constants[1]+prnak*_State[19]));
    double xs1ss = 1.0/(1.0+exp(-(_State[0]-1.50)/16.70));
    double xs2ss = xs1ss;

    double tauxs1;
    if (fabs(_State[0]+30.0)<0.001/0.0687)
        tauxs1 = 1/(0.0000719/0.148+0.000131/0.0687);
    else
        tauxs1 = 1.0/(0.0000719*(_State[0]+30.0)/(1.0-exp(-0.148*(_State[0]+30.0)))+0.000131*(_State[0]+30.0)/(exp(0.0687*(_State[0]+30.0))-1.0));
    double tauxs2 = 4*tauxs1;
    double gksx = (1+0.8/(1+pow((0.5/_State[9]),3)));

    double xiks = _Constants[10]*gksx*_State[5]*_State[6]*(_State[0]-eks);

    _Rates[5] = (xs1ss - (xs1ss - _State[5])*exp(-dt/tauxs1) - _State[5])/dt;
    _Rates[6] = (xs2ss - (xs2ss - _State[6])*exp(-dt/tauxs2) - _State[6])/dt;

    return xiks;
} // iks



//------Ik1 following Luo-Rudy formulation (from Shannon model) ------
double Mahajan_fail::comp_ik1 ( ){
    // Zhilin's ek:
    // double ek = (1.0/_Constants[24])*log(_Constants[2]/_Constants[1]);// K reversal potential

    // Our ek:
    double P_K = 100.0;
    double P_Na = 1.0;
    double ek = (1.0/_Constants[24])*log( (P_K*_Constants[2]+P_Na*_Constants[0])/(P_K*_Constants[1]+P_Na*_State[19])); // K reversal potential

    const double gki = (sqrt(_Constants[2]/5.4));
    double aki=1.02/(1.0+exp(0.2385*(_State[0]-ek-59.215)));
    double bki=(0.49124*exp(0.08032*(_State[0]-ek+5.476))+exp(0.061750*(_State[0]-ek-594.31)))/(1.0+exp(-0.5143*(_State[0]-ek+4.753)));
    double xkin=aki/(aki+bki);
    double xik1=_Constants[11]*gki*xkin*(_State[0]-ek);
    return xik1;
} // comp ik1



//------- Ito fast following Shannon et. al. 2005 -----------
//------- Ito slow following Shannon et. al. 2005 -----------
double Mahajan_fail::comp_ito(double dt) {
    double ek = 1.0/_Constants[24]*log(_Constants[2]/_Constants[1]);// K reversal potential

    double rt1=-(_State[0]+3.0)/15.0;
    double rt2=(_State[0]+33.5)/10.0;
    double rt3=(_State[0]+60.0)/10.0;
    double xtos_inf=1.0/(1.0+exp(rt1));
    double ytos_inf=1.0/(1.0+exp(rt2));
    double rs_inf=1.0/(1.0+exp(rt2));
    double txs=9.0/(1.0+exp(-rt1)) + 0.5;
    double tys=3000.0/(1.0+exp(rt3)) + 30.0;

    double xitos=_Constants[6]*_State[22]*(_State[23]+0.5*rs_inf)*(_State[0]-ek);// ito slow

    _Rates[22] = (xtos_inf - (xtos_inf-_State[22])*exp(-dt/txs) -_State[22])/dt;
    _Rates[23] = (ytos_inf - (ytos_inf-_State[23])*exp(-dt/tys) -_State[23])/dt;

    double xtof_inf=xtos_inf;
    double ytof_inf=ytos_inf;
    double rt4=-(_State[0]/30.0)*(_State[0]/30.0);
    double rt5=(_State[0]+33.5)/10.0;
    double txf=3.5*exp(rt4)+1.5;
    double tyf=20.0/(1.0+exp(rt5))+20.0;

    double xitof=_Constants[7]*_State[20]*_State[21]*(_State[0]-ek);// ito fast

    _Rates[20] = (xtof_inf - (xtof_inf-_State[20])*exp(-dt/txf) -_State[20])/dt;
    _Rates[21] = (ytof_inf - (ytof_inf-_State[21])*exp(-dt/tyf) -_State[21])/dt;

    return xitos+xitof;
} // comp ito



// -------Inak (sodium-potassium exchanger) following Shannon --------------
double Mahajan_fail::comp_inak ( ){
    const double xkmko=1.5;	//these are Inak constants adjusted to fit
                              //the experimentally measured dynamic restitution curve
    const double xkmnai=12.0;
    const double sigma = (exp(_Constants[0]/67.3)-1.0)/7.0;
    double fnak = 1.0/(1+0.1245*exp(-0.1*_State[0]*_Constants[24])+0.0365*sigma*exp(-_State[0]*_Constants[24]));
    double xinak = _Constants[12]*fnak*(1./(1.+(xkmnai/_State[19])))*_Constants[2]/(_Constants[2]+xkmko);
    return xinak;
} // inak



// --- Inaca (sodium-calcium exchange) following Shannon and Hund-Rudy------
//     Note: all concentrations are in mM
double Mahajan_fail::comp_inaca () {
    double csm = _State[8]/1000.0;// convert micro M to mM
    double zw3=pow(_State[19],3)*_Constants[3]*exp(_State[0]*0.35*_Constants[24])-pow(_Constants[0],3)*csm*exp(_State[0]*(0.35-1.)*_Constants[24]);
    double zw4=1.0+0.2*exp(_State[0]*(0.35-1.0)*_Constants[24]);
    const double xkdna=0.3;// micro M
    double aloss=1.0/(1.0+pow((xkdna/_State[8]),3));

    const double xmcao=1.3;
    const double xmnao=87.5;
    const double xmnai=12.3;
    const double xmcai=0.0036;

    double yz1=xmcao*pow(_State[19],3)+pow(xmnao,3)*csm;
    double yz2=pow(xmnai,3)*_Constants[3]*(1.0+csm/xmcai);
    double yz3=xmcai*pow(_Constants[0],3)*(1.0+pow((_State[19]/xmnai),3));
    double yz4=pow(_State[19],3)*_Constants[3]+pow(_Constants[0],3)*csm;
    double zw8=yz1+yz2+yz3+yz4;
    double xinacaq=_Constants[8]*aloss*zw3/(zw4*zw8);

    return xinacaq;
} // comp inaca



//	compute driving force (iCa)
double Mahajan_fail::comp_rxa() {
    double csm = _State[8]/1000.0;// convert micro M to mM
    const double pca=0.00054;
    double za=_State[0]*2.0*_Constants[24];
    double factor1=4.0*pca*_Constants[23]*_Constants[23]/(_Constants[22]*_Constants[20]);
    double factor=_State[0]*factor1;
    double rxa;
    if(fabs(za)<0.001)
    {
        rxa=factor1*(csm*exp(za)-0.341*(_Constants[3]))/(2.0*_Constants[24]);
    }
    else
    {
        rxa=factor*(csm*exp(za)-0.341*(_Constants[3]))/(exp(za)-1.0);
    }
    return rxa;
} // comp rxa



// ------ Markovian Ca current --------------------------------
//	Markov model:All parameters have been fitted directly to
//	experimental current traces using a multidimensional current fitting
//	routine.
double Mahajan_fail::comp_icalpo(double dt) {
    const double vth=0.0;
    const double s6=8.0;
    const double taupo=1.0;
    double poinf=1.0/(1.0+exp(-(_State[0]-vth)/s6));

    double alpha=poinf/taupo;
    double beta=(1.0-poinf)/taupo;

    const double r1=0.30;
    const double r2=3.0;
    const double cat=3.0;
    double fca=1.0/(1.0+pow(double(cat/_State[7]),3));
    double s1=0.0182688*fca;
    const double s1t=0.00195;
    double k1=0.024168*fca;
    const double k2=1.03615e-4;
    const double k1t=0.00413;
    const double k2t=0.00224;
    double s2=s1*(r1/r2)*(k2/k1);
    const double s2t=s1t*(r1/r2)*(k2t/k1t);

    const double vx=-40;
    const double sx=3.0;
    double poi=1.0/(1.0+exp(-(_State[0]-vx)/sx));
    const double tau3=3.0;

    double k3=(1.0-poi)/tau3;
    double k3t=k3;

    const double vy=-40.0;
    const double sy=4.0;
    double Pr=1.0-1.0/(1.0+exp(-(_State[0]-vy)/sy));

    double recov=10.0+4954.0*exp(_State[0]/15.6);

    const double tca=78.0329;
    const double cpt=6.09365;
    double tau_ca=tca/(1.0+pow((_State[7]/cpt),4))+0.1;

    double tauca=(recov-tau_ca)*Pr+tau_ca;
    double tauba=(recov-450.0)*Pr+450.0;

    const double vyr=-40.0;
    const double syr=11.32;
    double Ps=1.0/(1.0+exp(-(_State[0]-vyr)/syr));

    double k6=fca*Ps/tauca;
    double k5=(1.0-Ps)/tauca;

    double k6t=Ps/tauba;
    double k5t=(1.0-Ps)/tauba;

    double k4=k3*(alpha/beta)*(k1/k2)*(k5/k6);
    double k4t=k3t*(alpha/beta)*(k1t/k2t)*(k5t/k6t);

    double po=1.0-_State[15]-_State[17]-_State[16]-_State[18]-_State[13]-_State[14];
    _Rates[14]= beta*_State[13]+k5*_State[17]+k5t*_State[18]-(k6+k6t+alpha)*_State[14];
    _Rates[13]=alpha*_State[14]+k2*_State[15]+k2t*_State[16]+r2*po-(beta+r1+k1t+k1)*_State[13];
    _Rates[15]=k1*_State[13]+k4*_State[17]+s1*po-(k3+k2+s2)*_State[15];
    _Rates[17]=k3*_State[15]+k6*_State[14]-(k5+k4)*_State[17];
    _Rates[16]=k1t*_State[13]+k4t*_State[18]+s1t*po-(k3t+k2t+s2t)*_State[16];
    _Rates[18]=k3t*_State[16]+k6t*_State[14]-(k5t+k4t)*_State[18];

    return po;
} // comp icalpo



//----- SERCA2a uptake current ------------------------------------
double Mahajan_fail::comp_iup() {
    const double cup = 0.5;// uptake threshold

    return _Constants[13]*_State[9]*_State[9]/(_State[9]*_State[9] + cup*cup );
}  // comp iup



// ---------leak from the SR--------------------------
double Mahajan_fail::comp_ileak(){
  const double gleak = 0.00002069; // if there, *1.5 is for test only with failing cell
    const double kj = 50.0;
    return gleak * _State[10]*_State[10]/(_State[10]*_State[10]+kj*kj) * (_State[10]*16.667-_State[9]); //vsr/vcell=0.06
} //  comp ileak



// ---------- buffer dynamics in the myoplasm -----------------------
//buffering to calmodulin and SR are instantaneous, while buffering to
//Troponin C is time dependent.These are important to have reasonable
//Ca transient.Note: we have buffering in the suB_membrane space and
//the myoplasm.
double Mahajan_fail::comp_inst_buffer(double c) {
    const double bcal=24.0;
    const double xkcal=7.0;
    const double srmax=47.0;
    const double srkd=0.6;
    const double bmem=15.0;
    const double kmem=0.3;
    const double bsar=42.0;
    const double ksar=13.0;
    double bpx=bcal*xkcal/((xkcal+c)*(xkcal+c));
    double spx=srmax*srkd/((srkd+c)*(srkd+c));
    double mempx=bmem*kmem/((kmem+c)*(kmem+c));
    double sarpx=bsar*ksar/((ksar+c)*(ksar+c));
    return 1.0/(1.0+bpx+spx+mempx+sarpx);
} // comp inst_buffer



// --------- release-load functional dependence ----------------
double Mahajan_fail::comp_Q() {
    double bv=(_Constants[19]-_Constants[26])-_Constants[18]*_Constants[19];
    double Qr=0;
    if (_State[11]>_Constants[26] && _State[11]<_Constants[19])
    {
        Qr=_State[11]-_Constants[26];
    }
    else if(_State[11]>=_Constants[19])
    {
        Qr=_Constants[18]*_State[11]+bv;
    }

    if (Qr >= _Constants[25])
        Qr = _Constants[25];

    return _State[10]*Qr/_Constants[19];
} // comp Q



double Mahajan_fail::comp_dir(double po, double Qr, double rxa, double dcj) {
    const double ay=0.05;
    double sparkV=exp(-ay*(_State[0]+30))/(1.+exp(-ay*(_State[0]+30)));
    const double gryr=2.58079;
    double spark_rate=gryr*po*fabs(rxa)*sparkV;
    return spark_rate*Qr-_State[12]*(1-_Constants[16]*_Rates[10]/_State[10])/_Constants[16];
} // comp dir



// ----------- dyadic junction dynamics ------------------------
double Mahajan_fail::comp_dcp(double po, double Qr, double rxa) {
    const double gbarsr=26841.8;// m mol/(cm C)
    const double ax=0.3576;
    const double gdyad=9000.0;// m mol/(cm C)
    double gsr=gbarsr*exp(-ax*(_State[0]+30))/(1.0+exp(-ax*(_State[0]+30)));
    double xirp=po*Qr*fabs(rxa)*gsr;

    double xicap=po*gdyad*fabs(rxa);
    const double taups=0.5;
    return xirp+xicap-(_State[7]-_State[8])/taups;
} // comp dcp


  // Runge Kutta Update
  void Mahajan_fail::UpdateStateVariables(const Real dt, const Real i_stim) {
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
  Real Mahajan_fail::Compute_Ion(Real Xi, bool userXi, Real C_m, Real dt,
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
