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
	
	//test:
	const Int_t NBins = 100; 
	TH1D *hmP = ("hmP", "^{8}Be^{*} Mass", NBins, -0.2, mPMax + 0.2);
	hmP->GetXaxis()->SetTitle("Mass [MeV]");
	hmP->GetYaxis()->SetTitle("Counts");
	
	hmP->Draw();
	



}
