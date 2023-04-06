void Run3New(){

	//everything is in MeV, c = 1 (Natural units)

	const Int_t NEvents = 1000;	//number of decays: 1000
	float c = 1.0;	//c: speed of light = 1
	float pi = TMath::Pi();

//REST FRAME: PARENT


	float mP[NEvents], EP[NEvents], mC1[NEvents], mC2[NEvents], EC1[NEvents], EC2[NEvents];
	float pC1x[NEvents], pC1y[NEvents], pC1z[NEvents];	
	float pC2x[NEvents], pC2y[NEvents], pC2z[NEvents];

	float mPMean = 8023.15;		//8Be* mass (MeV)
	float mPSigma = 3.9E-11;	//8Be* stable

	float mC2Mean = 8005.30510;	//given base value of 8Be isotope
	float mC2Sigma = 3.9E-11;	//8Be also stable

	float mC1Mean = 16.98;	//X17 ==>invariant mass approx. 17 MeV
	float mC1Sigma = 1.0; 	//will play with later?

	float me = 0.511;	//mass of electron/positron = 0.511 MeV


	//GENERATING P-----------------------------------------------------------------------------------------------------------------------------------------------------
		//mP: in RF(P), just from Lorentzian Dist., w/ given avg. & width
			//Breit-Wigner
	
	float mPNorm = (2*mPSigma*mPSigma)/((fabs(mPSigma))*(fabs(mPSigma))*(fabs(mPSigma)));
	float cutoff = 0.10; 	//not taking values less probable than cutoff
	float mPMax = mPMean + (sqrt(((mPNorm*mPSigma)/(pi*cutoff))-(mPSigma*mPSigma)));

	TF1 *BWP = new TF1("BWP","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mPMax);
	BWP->SetParameters(pi, mPMean, mPSigma, mPNorm);

	float chance, prob;
	int mPindex = 0;

	float mPplaceholder;
	while(mPindex < NEvents + 1){
	
		mPindex++;

		mPplaceholder = mPMax*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax
	
		prob = gRandom->Rndm();		//U(x)
		chance = BWP->Eval(mPplaceholder);	//P(x)	
	
		if(prob > chance) {
			mPindex = mPindex - 1;
		}

		else if(prob < chance) {
			mP[mPindex] = mPplaceholder;
		}
	}
	//prob < chance -->accept, if not reject
	//repeat until we get all the parents we need
	
		//EP: in RF(P) from rest mass
	for(int i = 0; i < NEvents; i++){
		EP[i] = c*c*mP[i];	
	}
	
	//test:will use
	/*
	const Int_t NBins = 1000; 
	TH1D *hmP = new TH1D("hmP", "^{8}Be^{*} Mass", NBins, mPMax - 0.2, mPMax + 0.2);
	hmP->GetXaxis()->SetTitle("Mass [MeV]");
	hmP->GetYaxis()->SetTitle("Counts");
	
	for(int i = 0; i < NEvents; i++){
		hmP->Fill(mP[i]);
	}

	hmP->Draw();
	
	*/
	//GENERATING C2------------------------------------------------------------------------------------------------------------------------------------------------------
	
	float mC2Norm = (2*mC2Sigma*mC2Sigma)/((fabs(mC2Sigma))*(fabs(mC2Sigma))*(fabs(mC2Sigma)));
	//float cutoff = 0.10; 	//not taking values less probable than cutoff
	float mC2Max = mC2Mean + (sqrt(((mC2Norm*mC2Sigma)/(pi*cutoff))-(mC2Sigma*mC2Sigma)));

	TF1 *BWC2 = new TF1("BWC2","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mC2Max);
	BWC2->SetParameters(pi, mC2Mean, mC2Sigma, mC2Norm);
	float mC2placeholder;
	int mC2index = 0;
	while(mC2index < NEvents + 1){
	
		mC2index++;

		mC2placeholder = mC2Max*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax
	
		prob = gRandom->Rndm();		//U(x)
		chance = BWC2->Eval(mC2placeholder);	//P(x)	
	
		if(prob > chance) {
			mC2index = mC2index - 1;
		}

		else if(prob < chance) {
			mC2[mC2index] = mC2placeholder;
		}
	}
	
	/*
	const Int_t NBins = 1000; 
	TH1D *hmC2 = new TH1D("hmC2", "^{8}Be Mass", NBins, mC2Max - 0.2, mC2Max + 0.2);
	hmC2->GetXaxis()->SetTitle("Mass [MeV]");
	hmC2->GetYaxis()->SetTitle("Counts");
	
	for(int i = 0; i < NEvents; i++){
		hmC2->Fill(mC2[i]);
	}

	hmC2->Draw();
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HOW TO GET MAX.VALUE FROM AN ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	cout << "mC2Max:" << *max_element(mC2,mC2+NEvents) << endl;
	*/
	
		//EC2: C2 = 8Be approx. at rest ==> EC2 approx. rest mass
	
	for(int i = 0; i < NEvents; i++){
		EC2[i] = c*c*mC2[i];	
	}
	
	//GENERATING C1----------------------------------------------------------------------------------------------------------------------------------------------------------

		//mC1
	float mC1Norm = (2*mC1Sigma*mC1Sigma)/((fabs(mC1Sigma))*(fabs(mC1Sigma))*(fabs(mC1Sigma)));
	//float cutoff = 0.10; 	//not taking values less probable than cutoff
	float mC1Max = mC1Mean + (sqrt(((mC1Norm*mC1Sigma)/(pi*cutoff))-(mC1Sigma*mC1Sigma)));

	TF1 *BWC1 = new TF1("BWC1","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mC1Max);
	BWC1->SetParameters(pi, mC1Mean, mC1Sigma, mC1Norm);
	float mC1placeholder;
	int mC1index = 0;
	while(mC1index < NEvents + 1){
	
		mC1index++;

		mC1placeholder = mC1Max*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax
	
		prob = gRandom->Rndm();		//U(x)
		chance = BWC1->Eval(mC1placeholder);	//P(x)	
	
		if(prob > chance) {
			mC1index = mC1index - 1;
		}

		else if(prob < chance) {
			mC1[mC1index] = mC1placeholder;
		}
	}
	
	
	/*
	const Int_t NBins = 1000; 
	TH1D *hmC1 = new TH1D("hmC1", "X17 Mass", NBins, 0.0, mC1Max + 5.0);
	hmC1->GetXaxis()->SetTitle("Mass [MeV]");
	hmC1->GetYaxis()->SetTitle("Counts");
	
	for(int i = 0; i < NEvents; i++){
		hmC1->Fill(mC1[i]);
	}

	hmC1->Draw();
	*/
	
		//EC1: from energy conservation; gap btwn P & C2
	for(int i = 0; i < NEvents; i++){
		EC1[i] = fabs(EP[i] - EC2[i]);	//!!check later if we get issues w/ EP > EC2
	}
		//pC1: use energy conservation, relativistic eqn.
	for(int i = 0; i < NEvents; i++){
		pC1x[i] = c*sqrt();
		//pC1 = c*sqrt( (mP-mC2)^2 - (mC1)^2 )
	}
	
	


}
