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
  float dot[N], cosine[N], pG11New[N], pG12New[N];
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
	
	float mPNorm = (2*mPSigma*mPSigma)/((fabs(mPSigma))*(fabs(mPSigma))*(fabs(mPSigma)));
	float cutoff = 0.10; 	//not taking values less probable than cutoff
	mPMax = mPMean + (sqrt(((mPNorm*mPSigma)/(pi*cutoff))-(mPSigma*mPSigma)));

		TF1 *BWP = new TF1("BWP","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mPMax);
		BWP->SetParameters(pi, mPMean, mPSigma, mPNorm);

	float chance, prob;
	int mPindex = 0;

	float mPplaceholder;
		
	while(mPindex < N + 1){
	
		mPindex++;

		mPplaceholder = mPMax*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax
	
		prob = gRandom->Rndm();		//U(x)
		chance = BWP->Eval(mPplaceholder);	//P(x)	
	
		if(prob > chance) {
			mPindex = mPindex - 1;
			continue;
		}

		else if(prob < chance) {
			mP[mPindex] = mPplaceholder;
			
			if(mPindex == 0.125*N) std::cout << "12.5%" << std::endl;
			if(mPindex == 0.250*N) std::cout << "25.0%" << std::endl;
			if(mPindex == 0.375*N) std::cout << "37.5%" << std::endl;
			if(mPindex == 0.500*N) std::cout << "50.0%" << std::endl;
			if(mPindex == 0.625*N) std::cout << "62.5%" << std::endl;
			if(mPindex == 0.750*N) std::cout << "75.0%" << std::endl;
			if(mPindex == 0.875*N) std::cout << "87.5%" << std::endl;
			if(mPindex == N-1) std::cout << "PARENT GENERATION COMPLETE" << std::endl;
		}
		
			
	}

	for(int i = 0; i < N; i++){
		EP[i] = c*c*mP[i];
	}
	
	TCanvas *c1 = new TCanvas("c1","",900,900);
	c1->Divide(2,1);
	
	TH1D *hmP = new TH1D("hmP","",NBins, (*min_element(mP,mP+N)) - 100, (*max_element(mP,mP+N)) + 100);
	TH1D *hEP = new TH1D("hEP","",NBins, (*min_element(EP,EP+N)) - 100, (*max_element(EP,EP+N)) + 100);
	
	for(int i = 0; i < N; i++){
	hmP->Fill(mP[i]);
	hEP->Fill(EP[i]);
	}
	
	c1->cd(1);
	hmP->Draw();
	
	c1->cd(2);
	hEP->Draw();
	
	c1->Draw();



	/*
  int i = 0;
  while(i < N){
    i++;
  }
  */
    
    
}
