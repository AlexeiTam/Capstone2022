void X17Simulation(){

//===================== INITIALIZING VARIABLES, ARRAYS, CONSTANTS ==================================================
//==================================================================================================================

const Int_t N = 1000;	//Number of events = 1,000
const Int_t NBins = 1000;	//for histograms
float pi = TMath::Pi();
float c = 1.0;
  
    //mass generation parameters
    float mPMean = 8023.15;		//8Be* mass (MeV)
	  float mPSigma = 3.9E-11;	//8Be* stable

	  float mC2Mean = 8005.30510;	//given base value of 8Be isotope
	  float mC2Sigma = 3.9E-11;	//8Be also stable

	  float mC1Mean = 16.98;	//X17 ==>invariant mass approx. 17 MeV
	  float mC1Sigma = 1.0; 	//will play with later?

	  float me = 0.511;	//mass of electron/positron = 0.511 MeV
  
//================================ REST FRAME: 8Be* ================================
  
  float mP[N], mC2[N];
  float EP[N], EC2[N];
  float EC[N];
  float mC1[N];
  float pC1x[N], pC1y[N], pC1z[N];
  float beta[N], gamma[N];
  float L[N][4][4];
  float PC[N][4];

//================================ REST FRAME: X17 ================================
  float PCNew[N][4];
  float ECNew[N], pCNewx[N], pCNewy[N], pCNewz[N];
  float EG11[N], EG12[N];
  float pG11x[N], pG11y[N], pG11z[N]; 
  float pG12x[N], pG12y[N], pG12z[N];
  float LNew[N][4][4];
  float PG11[N][4], PG12[N][4];

//================================ REST FRAME: 8Be* (AGAIN) ================================
  float PG11New[N][4], PG12New[N][4];
  float PGNewTotal[N][4];
  float EG11New[N], EG12New[N];
  float pG11Newx[N], pG11Newy[N], pG11Newz[N];
  float pG12Newx[N], pG12Newy[N], pG12Newz[N];
  float dot[N], cosine[N], pG11New[N], pG11New[N], pG12New[N];
  float mTotal[N], pTotal[N];
  
  
  float mPMin, mPMax;
  float EPMin, EPMax;
  float mC1Min, mC1Max;
  float EC1Min, EC1Max;
  float pC1Min, pC1Max;
  float mCNewMin, mCNewMax;
  float ECNewMin, ECNewMax;
  float EG11Min, EG11Max;
  float pG11Min, pG11Max;
  float mTotalMin, mTotalMax;
  
 
  
  //------------------------------------------------GENERATING P--------------------------------------------------------------------
  
  int i = 0;
  while(i < N){
    i++;
  }
  
    
    
}
