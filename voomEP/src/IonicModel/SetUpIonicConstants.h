//-*-C++-*-

/*! \brief
  Set up parameters for Ionic model. We will pass these value to the 
  Ionic Model Class. Enough to Create it once and pass a pointer to
  all Ionic Class objects. Makes it memory efficient. Used by Body Class
  For consistency this resides in the Ionic Class directory
  
  Purkine Model
  Tusscher Model
  Mahajan (UCLA) Model
*/

namespace voom {
  //! Set up the Constants for Purkinje Model. For info check IonicModel class
  void SetUpPurkinjeParameters(Real **constants) {
    Real *Constants;
    Constants = new Real[63];
 
    Constants[0] = 310.0;
    Constants[1] = 96485.3415;
    Constants[2] = 8314.472;
    Constants[3] = 2.0;
    Constants[4] = 140;
    Constants[5] = 5.4;
    Constants[6] = 69.0;
    Constants[7] = 13266.5;
    Constants[8] = 0.6;
    Constants[9] = 0.2;
    Constants[10] = 0.06;
    Constants[11] = 700;     // iStim Start
    Constants[12] = 9000000; // iStim end
    Constants[13] = 500;     // iStim Period
    Constants[14] = 0.5;     // iStim duration
    Constants[15] = -4320.0; // iStim magnitude
    Constants[16] = 1.0;
    Constants[17] = 0.9;
    Constants[18] = 0.7;
    Constants[19] = 5.4;
    Constants[20] = 5.0;
    Constants[21] = 350.0;
    Constants[22] = 10.0;
    Constants[23] = 3.0;
    Constants[24] = 0.005;
    Constants[25] = 2.0;
    Constants[26] = 1140.0;
    Constants[27] = 15.0;
    Constants[28] = 2.0;
    Constants[29] = 20.0;
    Constants[30] = 1.5;
    Constants[31] = 3.0;
    Constants[32] = 3;
    Constants[33] = 0.001;
    Constants[34] = 0.001;
    Constants[35] = 0.5;
    Constants[36] = 442.2;
    Constants[37] = 0.01;
    Constants[38] = 0.01;
    Constants[39] = 0.0001;
    Constants[40] = 5.0;
    Constants[41] = 0.0001;
    Constants[42] = 1.5;
    Constants[43] = 0.188709677;
    Constants[44] = 0.045290323;
    Constants[45] = 2500.0;
    Constants[46] = 0.001;
    Constants[47] = 120.0;
    Constants[48] = 0.000246;
    Constants[49] = 1.7;
    Constants[50] = 1.6;
    Constants[51] = 120.0;
    Constants[52] = 0.000246;
    Constants[53] = 1.7;
    Constants[54] = 1.6;
    Constants[55] = 10.0;
    Constants[56] = 5000.0;
    Constants[57] = Constants[1]/( Constants[2]*Constants[0]);
    Constants[58] = ( Constants[2]*Constants[0])/Constants[1];
    Constants[59] =  Constants[9]*Constants[7];
    Constants[60] =  Constants[8]*Constants[7];
    Constants[61] =  (Constants[8]+Constants[9])*Constants[7];
    Constants[62] =  Constants[10]*Constants[7];

    *constants = Constants;
  }

  //! Constants for Tusscher Model. 
  void SetUpTusscherParameters(Real **constants) {
    Real *Constants;
    Constants = new Real[53];

    Constants[0] = 8314.472;
    Constants[1] = 310;
    Constants[2] = 96485.3415;
    Constants[3] = 0.185;
    Constants[4] = 0.016404;
    Constants[5] = 10;
    Constants[6] = 1000;
    Constants[7] = 1;
    Constants[8] = 52;
    Constants[9] = 0.03;
    Constants[10] = 5.4;
    Constants[11] = 140;
    Constants[12] = 2;
    Constants[13] = 5.405;
    Constants[14] = 0.153;
    Constants[15] = 0.392;
    Constants[16] = 14.838;
    Constants[17] = 0.00029;
    Constants[18] = 0.0000398;
    Constants[19] = 0.000592;
    Constants[20] = 0.294;
    Constants[21] = 2.724;
    Constants[22] = 1;
    Constants[23] = 40;
    Constants[24] = 1000;
    Constants[25] = 0.1;
    Constants[26] = 2.5;
    Constants[27] = 0.35;
    Constants[28] = 1.38;
    Constants[29] = 87.5;
    Constants[30] = 0.1238;
    Constants[31] = 0.0005;
    Constants[32] = 0.0146;
    Constants[33] = 0.15;
    Constants[34] = 0.045;
    Constants[35] = 0.06;
    Constants[36] = 0.005;
    Constants[37] = 1.5;
    Constants[38] = 2.5;
    Constants[39] = 1;
    Constants[40] = 0.102;
    Constants[41] = 0.0038;
    Constants[42] = 0.00025;
    Constants[43] = 0.00036;
    Constants[44] = 0.006375;
    Constants[45] = 0.2;
    Constants[46] = 0.001;
    Constants[47] = 10;
    Constants[48] = 0.3;
    Constants[49] = 0.4;
    Constants[50] = 0.00025;
    Constants[51] = 0.001094;
    Constants[52] = 0.00005468;
    
    *constants = Constants;
  }

  /*!
    Set up constants for UCLA Cell Model. Refer to Mahajan.h for description
    of the constants in the model.The pos value determines which location
    we need
         Apex Center Base
    =======================
    Epi    0     1     2
    Myo    3     4     5
    Endo   6     7     8
   */
  void SetUpMahajanParameters(Real **constants,  int pos) {
    Real *Constants = new Real[25];
    // Define Constants. RCONSTS is Other code
    Constants[0] = 136.0; //xnao    
    Constants[1] = 140.0; //xki
    Constants[2] = 5.4;   //xko
    Constants[3] = 1.8;   //cao
    Constants[4] = 182.0; //gca
    Constants[5] = 1.0;   //ica_factor
    Constants[6] = 0.04;  //0.094; //0.11;  //gtos


// Constants[7] = 0.11; 
// Constants[10] = 0.32;


    switch (pos) {
    
	  //control
    case 0: { Constants[7] = 0.11;     Constants[10] = 1.9*0.32;  break; }
    case 1: { Constants[7] = 0.11;     Constants[10] = 1.4*0.32;  break; }
    case 2: { Constants[7] = 0.11;     Constants[10] = 1.0*0.32;  break; }
    case 3: { Constants[7] = 0.11;     Constants[10] = 0.74*0.32; break; }
    case 4: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 5: { Constants[7] = 0.11;     Constants[10] = 0.35*0.32; break; }
    case 6: { Constants[7] = 0.11*.85; Constants[10] = 0.98*0.32; break; }
    case 7: { Constants[7] = 0.11*.85; Constants[10] = 0.7*0.32;  break; }
    case 8: { Constants[7] = 0.11*.85; Constants[10] = 0.5*0.32;  break; }
	
      /*
    // TRANSMURAL GRADIENTS ONLY:
    case 0: { Constants[7] = 0.11;     Constants[10] = 1.4*0.32;  break; }
    case 1: { Constants[7] = 0.11;     Constants[10] = 1.4*0.32;  break; }
    case 2: { Constants[7] = 0.11;     Constants[10] = 1.4*0.32;  break; }
    case 3: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 4: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 5: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 6: { Constants[7] = 0.11*.85; Constants[10] = 0.7*0.32;  break; }
    case 7: { Constants[7] = 0.11*.85; Constants[10] = 0.7*0.32;  break; }
    case 8: { Constants[7] = 0.11*.85; Constants[10] = 0.7*0.32;  break; }
      */

    /* 
    // APICOBASAL GRADIENTS ONLY:
    case 0: { Constants[7] = 0.11;     Constants[10] = 0.74*0.32; break; }
    case 1: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 2: { Constants[7] = 0.11;     Constants[10] = 0.35*0.32; break; }
    case 3: { Constants[7] = 0.11;     Constants[10] = 0.74*0.32; break; }
    case 4: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 5: { Constants[7] = 0.11;     Constants[10] = 0.35*0.32; break; }
    case 6: { Constants[7] = 0.11;     Constants[10] = 0.74*0.32; break; }
    case 7: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 8: { Constants[7] = 0.11;     Constants[10] = 0.35*0.32; break; }
    */
    default: 
    {
      std::cout << "Tissue Type not set correctly. Setting to Epi-Base\n";
      Constants[7] = 0.11;  Constants[10] = 0.32; 
    }
    }
	
    Constants[8] = 0.84;  //gnaca
    Constants[9] = 0.0125; //gkr
 	
    Constants[11] = 0.3;  //gk1
    Constants[12] = 1.5;  //gnak
    Constants[13] = 0.4;  //vup
    Constants[14] = 4.0;  //taus
    Constants[15] = 12.0; //gna
    Constants[16] = 30.0; //taur
    Constants[17] = 100.0; //taua
    Constants[18] = 11.3; //av
    Constants[19] = 90.0; //cstar
    Constants[20] = 308.; //T
    Constants[21] = -80.; //Vc
    Constants[22] = 8.314;//xxr
    Constants[23] =96.485;//xf
    Constants[24] = Constants[23]/(Constants[22]*Constants[20]); //frt

    *constants = Constants;
  }


void SetUpMahajan_failParameters(Real **constants,  int pos) {
    Real *Constants = new Real[25];
   //  Define Constants. RCONSTS is Other code
    Constants[0] = 136.0;//xnao
    Constants[1] = 140.0;//xki
    Constants[2] = 5.4;//xko
    Constants[3] = 1.8;//cao
    Constants[4] = 182.0;// 1000.; /182.0; /gca
    Constants[5] = 1.0;   //ica_factor
    Constants[6] = 0.04*0.65; //0.094; /0.11;  /gtos

 switch (pos) {
    case 0: { Constants[7] = 0.11;     Constants[10] = 1.9*0.32;  break; }
    case 1: { Constants[7] = 0.11;     Constants[10] = 1.4*0.32;  break; }
    case 2: { Constants[7] = 0.11;     Constants[10] = 1.0*0.32;  break; }
    case 3: { Constants[7] = 0.11;     Constants[10] = 0.74*0.32; break; }
    case 4: { Constants[7] = 0.11;     Constants[10] = 0.52*0.32; break; }
    case 5: { Constants[7] = 0.11;     Constants[10] = 0.35*0.32; break; }
    case 6: { Constants[7] = 0.11*.85; Constants[10] = 0.98*0.32; break; }
    case 7: { Constants[7] = 0.11*.85; Constants[10] = 0.7*0.32;  break; }
    case 8: { Constants[7] = 0.11*.85; Constants[10] = 0.5*0.32;  break; }
 
 /*
   // APICOBASAL GRADIENTS ONLY:
 case 0: { Constants[7] = 0.11;     Constants[10] = 1.2*0.32; break; }
 case 1: { Constants[7] = 0.11;     Constants[10] = 1.0*0.32; break; }
 case 2: { Constants[7] = 0.11;     Constants[10] = 0.8*0.32; break; }
 case 3: { Constants[7] = 0.11;     Constants[10] = 1.2*0.32; break; }
 case 4: { Constants[7] = 0.11;     Constants[10] = 1.0*0.32; break; }
 case 5: { Constants[7] = 0.11;     Constants[10] = 0.8*0.32; break; }
 case 6: { Constants[7] = 0.11;     Constants[10] = 1.2*0.32; break; }
 case 7: { Constants[7] = 0.11;     Constants[10] = 1.0*0.32; break; }
 case 8: { Constants[7] = 0.11;     Constants[10] = 0.8*0.32; break; }
 */
 default:
   {
     std::cout << "Tissue Type not set correctly. Setting to Epi-Base\n";
     
     Constants[7] = 0.11;  Constants[10] = 0.32;
   }
 }

    Constants[7] *=0.67; 
    Constants[8] = 0.84*2.0; //gnaca
    Constants[9] = 0.0125;   //gkr
    Constants[10] *= 0.6;
    Constants[11] = 0.3*0.5;//gk1
    Constants[12] = 1.5;    //gnak
    Constants[13] = 0.24;    //vup
    Constants[14] = 4.0;    //taus
    Constants[15] = 12.0;   //gna
    Constants[16] = 30.0;   //taur
    Constants[17] = 100.0;  //taua
    Constants[18] = 56.5;   //av
    Constants[19] = 60.0;   //cstar
    Constants[20] = 308.;   //T
    Constants[21] = -80.;   //Vc
    Constants[22] = 8.314;  //xxr
    Constants[23] = 96.485; //xf
    Constants[24] = Constants[23]/(Constants[22]*Constants[20]); //frt

    *constants = Constants;
  }

}
