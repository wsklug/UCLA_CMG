//-*-C++-*-
// Luo Rudy ionic model class 
#if !defined(__LRGT_h__)
#define __LRGT_h__
#include <iostream>
#include <math.h>
#include <stdlib.h>


using namespace std;
using namespace voom;

class LuoRudyGateTable  {
    
public:

    inline void makeGateTable() ;
    
    inline void get_parameters(Real v, Real tabz[]);
    inline void makeGateTableSlope();
    inline void make_parameter1(Real v,int index);
    Real **_gateTable;  //[129001][14];// 290001 = (200+90)/.001 + 1
    Real **_gateTableSlope;// [129001][14];

    //! default constructor
    LuoRudyGateTable(){
      Gna=12.00;
      Gsi=0.09; // correct value
      //Gsi =0.01; // use this for spiral
      c_K0=5.4;
      GK=.282*sqrt(c_K0/5.4);
      
      c_Nai=18.0;
      c_Na0=140.0;
      c_Ki=145.0;
      RT_F=26.73;
      E_Na=RT_F*log(c_Na0/c_Nai);
      E_K=RT_F*log((c_K0+.01833*c_Na0)/(c_Ki+.01833*c_Nai));
      E_K1=RT_F*log(c_K0/c_Ki);
      E_kp=E_K1;
      G_K1=.6047*sqrt(c_K0/5.4);
      
      _gateTable = new Real*[290001];
      
      for(int i = 0; i < 290001; ++i){
	_gateTable[i] = new Real[14];
      }
#ifdef _USE_GATE_SLOPE
      _gateTableSlope = new Real*[290000];
      for(int i = 0; i < 290000; ++i){
	_gateTableSlope[i] = new Real[14];
      }
#endif
      makeGateTable();
      
    }
         
    Real Gna,Gsi,c_K0,GK;
    Real Es;
    Real i_na,i_si,i_k,i_tot;
    Real c_Na0,c_Nai,c_Ki,RT_F;
    Real E_Na, E_K,E_K1,E_kp,G_K1; 
};



void LuoRudyGateTable::makeGateTable(){

  for(int index=0; index<290001; index++){    
    Real v = -90 + 0.001*index;
    make_parameter1(v,index);
  }

#ifdef _USE_GATE_SLOPE
  makeGateTableSlope();
#endif

  return;
}



void LuoRudyGateTable::makeGateTableSlope(){

  for(int index=0; index<290000; index++){
    for(int i =0;i<14; i++){
      _gateTableSlope[index][i] = (_gateTable[index+1][i] -_gateTable[index][i] )*1000;
    }
  }

  return;
}


void LuoRudyGateTable::get_parameters(Real v,Real tab[]){

  int index = (v + 90)*1000;

#ifdef _USE_GATE_SLOPE
  Real delv = v -(-90 + 0.001*index); 
  for(int i=0; i<14; i++){
    tab[i] = _gateTableSlope[index][i]*delv + _gateTable[index][i];
  }  
#else 
  Real delVbydT = 1000*(v -(-90 + 0.001*index)); 
  for(int i=0; i<14; i++){
    tab[i] = (_gateTable[index+1][i] -_gateTable[index][i] )*delVbydT + _gateTable[index][i];
  }
#endif
}

//void LuoRudyGateTable::make_parameter1(Real v,Real tab[]){
void LuoRudyGateTable::make_parameter1(Real v, int index){
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
    _gateTable[index][0]=m_inf;
    _gateTable[index][1]=tao_m;
    _gateTable[index][2]=h_inf;
    _gateTable[index][3]=tao_h;
    _gateTable[index][4]=j_inf;
    _gateTable[index][5]=tao_j;
    _gateTable[index][6]=d_inf;
    //c     tab[8]=tao_d*10.0
    _gateTable[index][7]=tao_d;
    _gateTable[index][8]=f_inf;
    //c     tab[10]=tao_f*10.0
    _gateTable[index][9]=tao_f;
    _gateTable[index][10]=x_inf;
    _gateTable[index][11]=tao_x;
    _gateTable[index][12]=x1;
    _gateTable[index][13]=I_k1+I_kp+I_b;
    return;
}




#endif
