//-*-C++-*-
// Luo Rudy ionic model class 
#include "LR.h"

namespace voom {    
  void LuoRudy::reinitialize(Real y[]){
    for(int i=0;i<7;i++) _in[i]=y[i];
  }
  
  void LuoRudy::setGna(Real val){
    Gna = val;
  }
  
  void LuoRudy::setGsi(Real val){
    Gsi = val;
  }
  
  void LuoRudy::initialize() {    
    _in[0] = .001;      // m
    _in[1] = 0.0001;    // j
    _in[2] = 0.99;      // h
    _in[3] = .001;      // d
    _in[4] = 0.9;       // f
    _in[5] = .001;      //
    _in[6] = 0.0000001; //
    _in[7] = 0;         // Voltage

    /*
      SurfaceAreaToVolumeRatio - Geometric data from LR Paper
      r = 11um = 11 *10^-4 cm
      L = 100um = 100 * 10^-4 cm
      Volume = pi*r^2*L
      SurfaceArea = 2*pi*r^2 + 2*pi*r*L
      Xi = (2*pi*r^2 + 2*pi*r*L)/pi*r^2*L = 2*(r+L)/(rL)
    */
    Real r = 11e-4, L = 100E-4;
    _Xi = 2.*(r + L)/(r*L);
    
    Gna=12.00;
    //    Gsi=0.09; // correct value
    // Modified
    Gsi = 0.07;
    c_K0 = 5.4; // 5.4;

    //    GK=.282*sqrt(c_K0/5.4);
    // Modified 
    GK = 0.705*sqrt(c_K0/5.4);
    
    c_Nai=18.0;
    c_Na0=140.0;
    c_Ki=145.0;
    RT_F=26.73;
    E_Na=RT_F*log(c_Na0/c_Nai);
    E_K=RT_F*log((c_K0+.01833*c_Na0)/(c_Ki+.01833*c_Nai));
    E_K1=RT_F*log(c_K0/c_Ki);
    E_kp=E_K1;
    G_K1=.6047*sqrt(c_K0/5.4);

  }
  
  /*! 
    Luo Rudy model compute current as uA/uF. 
    We convert istim to get the same units and return dt/dt
  */
  Real LuoRudy::Compute_Ion(Real Xi, bool userXi, Real C_m,Real dt, Real volt,
			    Real istim){
    _in[7]=volt;
    Real tabz[14];
    //    std::cout << "iStim: " << istim/(Xi*C_m) << " ";
    //! If use LR Gate model then
    if (_useGate) {
      if(volt>-90 && volt <200 ){
	gate->get_parameters(_in[7], tabz);
      }else{
	make_parameter1(_in[7],tabz);
      }
    }
    // No Gate model requested
    else
      make_parameter1(_in[7],tabz);
    
    
    for(int i=0;i<6;i++){
      _in[i]=tabz[2*i]-(tabz[2*i]-_in[i])*exp(-dt*tabz[2*i+1]) ;
    }
    
    Real Es=7.7-13.0287*log(_in[6]);
    Real i_na=Gna*_in[0]*_in[0]*_in[0]*_in[1]*_in[2]*(_in[7]-54.4);
    Real i_si=Gsi*_in[3]*_in[4]*(_in[7]-Es);
    //    const Real L = 2.;
    //    Real alpha = (1. - 0.2*x/L);
    Real i_k= GK*tabz[12]*_in[5]*(_in[7]+77.0);  
    _in[6]=_in[6]+(-.0001*i_si+0.07*(0.0001-_in[6]))*dt;
    Real i_tot = -(tabz[13]+i_k+i_na+i_si);

    // istim is in ua/CC. istim/(Xi*C) = ua/cc /(1/cm*uF/cm^2) = uA/uF
    //std::cout << (i_tot + istim/(Xi*C_m)) << "\n";
    //if (userXi) return (i_tot + istim/(Xi*C_m));
    return (i_tot + istim/(_Xi*C_m));
  }
  
  void LuoRudy::make_parameter1(Real v,Real tab[]){
    Real beta1,beta2,beta3;
    
    Real alpha_m, beta_m,alpha_h,beta_h,alpha_j,beta_j;
    Real alpha_d, beta_d,alpha_f,beta_f,alpha_x,beta_x;
    Real alpha_k1,beta_k1,k1_inf,kp,I_k1,I_kp,I_b,I_total;
    Real m_inf,h_inf,j_inf,d_inf,f_inf,x_inf;
    Real tao_m,tao_h,tao_j,tao_d,tao_f,tao_x,x1;
    
    if (fabs(v+47.13)<0.001) {
      alpha_m=3.2;
    }else{
      alpha_m=0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
    }
    beta_m=0.08*exp(-v/11.);
    if (v>-40.) {
      alpha_h=0;
      alpha_j=0;
      beta_h=1./(0.13*(1+exp(-(v+10.66)/11.1)));
      beta_j=0.3*exp(-2.535e-7*v)/(1+exp(-0.1*(v+32)));
    }else{
      alpha_h=0.135*exp(-(80+v)/6.8);
      alpha_j=(-1.2714e5*exp(0.2444*v)-3.474e-5*exp(-0.04391*v))*(v+37.78)/(1+exp(0.311*(v+79.23)));
      beta_h=3.56*exp(0.079*v)+3.1e5*exp(0.35*v);
      beta_j=0.1212*exp(-0.01052*v)/(1+exp(-0.1378*(v+40.14)));
    } 
    
    alpha_d=.095*exp(-.01*(v-5))/(1+exp(-.072*(v-5)));
    beta_d=.07*exp(-.017*(v+44))/(1+exp(.05*(v+44)));
    alpha_f=.012*exp(-.008*(v+28))/(1+exp(.15*(v+28)));
    beta_f=.0065*exp(-.02*(v+30))/(1+exp(-.2*(v+30)));
    alpha_x=.0005*exp(.083*(v+50))/(1+exp(.057*(v+50)));
    beta_x=.0013*exp(-.06*(v+20))/(1+exp(-.04*(v+20)));
    tao_m=alpha_m+beta_m;
    m_inf=alpha_m/tao_m;
    tao_h=alpha_h+beta_h;
    h_inf=alpha_h/tao_h;
    tao_j=alpha_j+beta_j;
    j_inf=alpha_j/tao_j;
    tao_d=(alpha_d+beta_d);
    d_inf=alpha_d/(alpha_d+beta_d);
    tao_f=(alpha_f+beta_f);
    f_inf=alpha_f/(alpha_f+beta_f);
    tao_x=(alpha_x+beta_x);
    x_inf=alpha_x/(alpha_x+beta_x);
    if  (v>-100) {
      if  (fabs(v+77.)<0.001) {
	x1=0.608889;
      }else{
	x1=2.837*(exp(0.04*(v+77))-1)/((v+77)*exp(0.04*(v+35)));
      }
    }else{
      x1=1.;
    }
    alpha_k1=1.02/(1+exp(0.2385*(v-E_K1-59.215)));
    beta1=0.49124*exp(0.08032*(v-E_K1+5.476));
    beta2=exp(0.06175*(v-E_K1-594.31));
    beta3=1+exp(-0.5143*(v-E_K1+4.753));
    beta_k1=(beta1+beta2)/beta3;
    k1_inf=alpha_k1/(alpha_k1+beta_k1);
    I_k1=G_K1*k1_inf*(v-E_K1);
    kp=1./(1+exp((7.488-v)/5.98));
    I_kp=0.0183*kp*(v-E_K1);
    I_b=0.0392*(v+59.87);
    tab[0]=m_inf;
    tab[1]=tao_m;
    tab[2]=h_inf;
    tab[3]=tao_h;
    tab[4]=j_inf;
    tab[5]=tao_j;
    tab[6]=d_inf;

    tab[7]=tao_d;
    tab[8]=f_inf;
    //c     tab[10]=tao_f*10.0
    tab[9]=tao_f;
    tab[10]=x_inf;
    tab[11]=tao_x;
    tab[12]=x1;
    tab[13]=I_k1+I_kp+I_b;
    return;
  }  

  Real LuoRudy::getGamma() {
    // This works all the way
    // This gave good deformation with CompNeo Hookean
    // Works so far.
    // return 0.2222*(4. - 0.5*tanh(1.5*log(_in[6]) + 8.));
    // return 0.217*(4. - 0.6*tanh(1.5*log(_in[6]) + 9.));
    // Gamma going upto 60%
    return 0.2126*(4. - 0.7*tanh(1.4*log(_in[6]) + 9.5));
  }
};  // namespace voom
