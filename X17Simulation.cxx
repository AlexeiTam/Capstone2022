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
  
  float mP[N], mC2[N], dmPC2[N];
  float EP[N], EC2[N];
  float EC1[N];
  float mC1[N];
  float pC1x[N], pC1y[N], pC1z[N];
	float pC2x[N], pC2y[N], pC2z[N];
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
  
 float mC2Max;
  
  //------------------------------------------------GENERATING P--------------------------------------------------------------------
	
	float mPNorm = (2*mPSigma*mPSigma)/((fabs(mPSigma))*(fabs(mPSigma))*(fabs(mPSigma)));
	float cutoff = 0.10; 	//not taking values less probable than cutoff
	mPMax = mPMean + (sqrt(((mPNorm*mPSigma)/(pi*cutoff))-(mPSigma*mPSigma)));

		TF1 *BWP = new TF1("BWP","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mPMax);
		BWP->SetParameters(pi, mPMean, mPSigma, mPNorm);

	float chance, prob;
	int mPindex = 0;

	float mPplaceholder;
	
		
	while(mPindex < N){


		mPplaceholder = mPMax*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax
	
		prob = gRandom->Rndm();		//U(x)
		chance = BWP->Eval(mPplaceholder);	//P(x)	
	
		if(prob > chance) {
			continue;
		}

		else if(prob < chance) {
			mP[mPindex] = mPplaceholder;
			mPindex++;
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
	
//================================= C1,C2 MASS PARAMETERS ==================================================
	
	float mC1Norm = (2*mC1Sigma*mC1Sigma)/((fabs(mC1Sigma))*(fabs(mC1Sigma))*(fabs(mC1Sigma)));
	float mC2Norm = (2*mC2Sigma*mC2Sigma)/((fabs(mC2Sigma))*(fabs(mC2Sigma))*(fabs(mC2Sigma)));
	mC1Max = mC1Mean + (sqrt(((mC1Norm*mC1Sigma)/(pi*cutoff))-(mC1Sigma*mC1Sigma)));
	mC2Max = mC2Mean + (sqrt(((mC2Norm*mC2Sigma)/(pi*cutoff))-(mC2Sigma*mC2Sigma)));

		TF1 *BWC1 = new TF1("BWC1","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mC1Max);
		BWC1->SetParameters(pi, mC1Mean, mC1Sigma, mC1Norm);
	
		TF1 *BWC2 = new TF1("BWC2","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mC2Max);
		BWC2->SetParameters(pi, mC2Mean, mC2Sigma, mC2Norm);

	float C1prob, C1chance, C2prob, C2chance;
	float mC1placeholder;
	float mC2placeholder;
	
	
	//================================= GENERATING C2 =================================
	cout << "GENERATING C2" << endl;
	int mC2index = 0;
	
	
	while(mC2index < N){


		mC2placeholder = (mP[mC2index])*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax (ORIGINAL)
	  								//now, enforce limit for 8Be to not go above 8Be*
		if(mC2placeholder > mP[mC2index]){
		cout << "EVENT #" << mC2index <<": nonphysical mC2 rejected" << endl;
			continue;
		}
	
		C2prob = gRandom->Rndm();		//U(x)
		C2chance = BWC2->Eval(mC2placeholder);	//P(x)	
	
		if(C2prob > C2chance) {
			continue;
		}

		else if(C2prob < C2chance) {
			mC2[mC2index] = mC2placeholder;
			mC2index++;
			if(mC2index == 0.125*N) std::cout << "12.5%" << std::endl;
			if(mC2index == 0.250*N) std::cout << "25.0%" << std::endl;
			if(mC2index == 0.375*N) std::cout << "37.5%" << std::endl;
			if(mC2index == 0.500*N) std::cout << "50.0%" << std::endl;
			if(mC2index == 0.625*N) std::cout << "62.5%" << std::endl;
			if(mC2index == 0.750*N) std::cout << "75.0%" << std::endl;
			if(mC2index == 0.875*N) std::cout << "87.5%" << std::endl;
			if(mC2index == N-1) std::cout << "C2 GENERATION COMPLETE" << std::endl;
		}
		
			
	}
	
	for(int i = 0; i < N; i++){
		EC2[i] = c*c*mC2[i];
		dmPC2[i] = mP[i] - mC2[i];
	}
	
	cout << "GENERATING C1" << endl;
	
	for(int i = 0; i < N; i++){
		
		EC1[i] = EP[i] - EC2[i];
	}
	
	//generating mC1
	
	int index = 0;
	while(index < N){
	
		mC1[index] = BWC1->GetRandom();
		
		if(mC1[index] > dmPC2[index]) {
		continue;	
		}
		   
		   
		   	if(index == 0.125*N) std::cout << "12.5%" << std::endl;
			if(index == 0.250*N) std::cout << "25.0%" << std::endl;
			if(index == 0.375*N) std::cout << "37.5%" << std::endl;
			if(index == 0.500*N) std::cout << "50.0%" << std::endl;
			if(index == 0.625*N) std::cout << "62.5%" << std::endl;
			if(index == 0.750*N) std::cout << "75.0%" << std::endl;
			if(index == 0.875*N) std::cout << "87.5%" << std::endl;
			if(index == N-1) std::cout << "C1 GENERATION COMPLETE" << std::endl;
		index++;
	}
	
	
	/*
	int mC1index = 0;
	while(mC1index < N){


		mC1placeholder = (dmPC2[mC1index])*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax (ORIGINAL)
	  										//now, enforce limit for 8Be to not go above 8Be*
	
		C1prob = gRandom->Rndm();		//U(x)
		C1chance = BWC1->Eval(mC1placeholder);	//P(x)	
	
		if(C1prob > C1chance) {
			continue;
		}

		else if(C1prob < C1chance) {
			mC1[mC2index] = mC1placeholder;
			mC1index++;
			if(mC1index == 0.125*N) std::cout << "12.5%" << std::endl;
			if(mC1index == 0.250*N) std::cout << "25.0%" << std::endl;
			if(mC1index == 0.375*N) std::cout << "37.5%" << std::endl;
			if(mC1index == 0.500*N) std::cout << "50.0%" << std::endl;
			if(mC1index == 0.625*N) std::cout << "62.5%" << std::endl;
			if(mC1index == 0.750*N) std::cout << "75.0%" << std::endl;
			if(mC1index == 0.875*N) std::cout << "87.5%" << std::endl;
			if(mC1index == N-1) std::cout << "C1 GENERATION COMPLETE" << std::endl;
		}
		
			
	}
	
	*/
	
	
	cout<< "SCANNING FOR NONPHYSICAL (mP, mC1, mC2)" << endl;
	for(int i = 0; i < N; i++){
	
		if(mC2[i] > mP[i]){
			cout<<"i:"<<i<<"...mP("<<mP[i]<<")...mC2("<<mC2[i]<<")...mC1("<<mC1[i]<<")"<<endl;	
		}
		if(mC1[i] > dmPC2[i]){
			cout<<"i:"<<i<<"...mP("<<mP[i]<<")...mC2("<<mC2[i]<<")...mC1("<<mC1[i]<<")"<<endl;	
		}
	}
	
	//================================= CONTINUING GENERATION in P FRAME =================================
	
	//initialize L, LNew to 0
	for(int i = 0; i < N; i++){
	
		for(int j = 0; j < 4; j++){
	
			for(int k = 0; k < 4; k++){
				L[i][j][k] = 0.0;
				LNew[i][j][k] = 0.0;
			}
		}
	
	}
	
	
  int i = 0;
  while(i < N){
	  
	  pC1y[i] = 0.0;
	  pC1z[i] = 0.0;
	  pC2y[i] = 0.0;
	  pC2z[i] = 0.0;
	  
	  pC2x[i] = 0.0;
	  //pC1x[i] = sqrt((((EC1[i])*(EC1[i]))/(c*c))-((mC1[i])*(mC1[i])*(c*c)));
	  pC1x[i] = sqrt(((EC1[i])*(EC1[i]))-((mC1[i])*(mC1[i])));
	  //cout << pC1x[i] << endl;
	  
	  //beta[i] = (((pC1x[i])/(c))/(sqrt(((mC1[i])*(mC1[i]))+(((pC1x[i])*(pC1x[i]))/(c*c)))));
		beta[i] = ((c*pC1x[i])/(EC1[i]));// use p = gamma*mv + E = gamma*m*c^2 --> beta = v/c = pc/E
		gamma[i] = (sqrt((1.0)/(1.0 - ((beta[i])*(beta[i])))));
		
		L[i][0][0] = gamma[i];
		L[i][1][1] = gamma[i];
		L[i][0][1] = (-1.0)*(gamma[i])*(beta[i]);
		L[i][1][0] = (-1.0)*(gamma[i])*(beta[i]);
		L[i][2][2] = 1.0;
		L[i][3][3] = 1.0;
		
		//generate LNew, since it just goes the other way around
		LNew[i][0][0] = gamma[i];
		LNew[i][1][1] = gamma[i];
		LNew[i][0][1] = (gamma[i])*(beta[i]);
		LNew[i][1][0] = (gamma[i])*(beta[i]);
		LNew[i][2][2] = 1.0;
		LNew[i][3][3] = 1.0;
		
	  /*
		if(beta[i] != beta[i]){
			cout << "EVENT #" << i <<": beta is NaN!" << endl;
		}
		
		if(gamma[i] != gamma[i]){
			cout << "EVENT #" << i <<": gamma is NaN!" << endl;
		}
		*/
	
		//fill 4-vectors of C1 in P Frame
		PC[i][0] = EC1[i];
		PC[i][1] = pC1x[i];
		PC[i][2] = pC1y[i];
		PC[i][3] = pC1z[i];
		
    i++;	//on to next event
  }	//end of event generation
  
	//PRINTING OUT NAN CASES
	for(int i = 0; i < N; i++){
	
		if(pC1x[i] != pC1x[i]){
		cout << "NAN EVENT #" << i << "mP(" << mP[i] << ")...mC2(" << mC2[i] << ")...mC1(" << mC1[i] << ")" << endl;	
		}
	}
	
	//======================================== BOOST: P-->C1 =========================================
	//======================================== REST FRAME: C1 ========================================
	//initialize PCNew to 0
	
	for(int i = 0; i < N; i++){
		
		for(int j = 0; j < 4; j++){
			PCNew[i][j] = 0.0;
		}
			
	}
	
	//do the boost
	for(int i = 0; i < N; i++){
	
		for(int j = 0; j < 4; j++){
		
			for(int k = 0; k < 4; k++){
			PCNew[i][j] = PCNew[i][j] + ((L[i][j][k])*(PC[i][k]));
			}
		}
	}
	
	
	//scan for NaN in PCNew
	
	for(int i = 0; i < N; i++){
	
		for(int j = 0; j < 4; j++){
		
			if(PCNew[i][j] != PCNew[i][j]){
			cout << "NaN in PCNew[" << i <<"][" << j <<"]:" << PCNew[i][0] << "..." << PCNew[i][1] << "..." <<  PCNew[i][2] << "..." <<  PCNew[i][3] << endl;	
			}
		}
	}
	
	//=================TEST VISUALIZATION/=================
	
	TCanvas *c1 = new TCanvas("c1","",900,900);
	c1->Divide(3,2);
	
	TH1D *hmP = new TH1D("hmP","hmP",NBins, (*min_element(mP,mP+N)) - 100, (*max_element(mP,mP+N)) + 100);
	TH1D *hEP = new TH1D("hEP","hEP",NBins, (*min_element(EP,EP+N)) - 100, (*max_element(EP,EP+N)) + 100);
	TH1D *hmC2 = new TH1D("hmC2","hmC2",NBins, (*min_element(mC2,mC2+N)) - 100, (*max_element(mC2,mC2+N)) + 100);
	TH1D *hEC2 = new TH1D("hEC2","hEC2",NBins, (*min_element(EC2,EP+N)) - 100, (*max_element(EC2,EC2+N)) + 100);
	TH1D *hmC1 = new TH1D("hmC1","hmC1",NBins, (*min_element(mC1,mC1+N)) - 10, (*max_element(mC1,mC1+N)) + 100);
	TH1D *hEC1 = new TH1D("hEC1","hEC1",NBins, (*min_element(EC1,EC1+N)) - 10, (*max_element(EC1,EC1+N)) + 10);
	
	for(int i = 0; i < N; i++){
	hmP->Fill(mP[i]);
	hEP->Fill(EP[i]);
	hmC1->Fill(mC1[i]);
	hEC1->Fill(EC1[i]);
	hmC2->Fill(mC2[i]);
	hEC2->Fill(EC2[i]);
	}
	
	c1->cd(1);
	hmP->Draw();
	
	c1->cd(2);
	hEP->Draw();
	
	c1->cd(3);
	hmC1->Draw();
	
	c1->cd(4);
	hEC1->Draw();
	
	c1->cd(5);
	hmC2->Draw();
	
	c1->cd(6);
	hEC2->Draw();
	
	c1->Draw();
	
	

    
    
}
