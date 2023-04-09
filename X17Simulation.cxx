void X17Simulation(){

//===================== INITIALIZING VARIABLES, ARRAYS, CONSTANTS ==================================================
//==================================================================================================================

const Int_t N = 1000;	//Number of events = 1,000
const Int_t NBins = 1000;	//for histograms
float pi = TMath::Pi();
float c = 1.0;
  
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
  
  
    
    
}
