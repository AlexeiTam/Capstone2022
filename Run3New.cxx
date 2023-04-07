void Run3New(){

	//everything is in MeV, c = 1 (Natural units)

	const Int_t NEvents = 1000;	//number of decays: 1000
	const Int_t NBins = 1000;
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
	
		//pC2: at rest
	
		for(int i = 0; i < NEvents; i++){
			pC2x[i] = 0.0;
			pC2y[i] = 0.0;
			pC2z[i] = 0.0;
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
		EC1[i] = EP[i] - EC2[i];	//!!check later if we get issues w/ EP > EC2
	}
		//pC1: use energy conservation, relativistic eqn.
	for(int i = 0; i < NEvents; i++){
		//pC1x[i] = c*sqrt(((mP[i] - mC1[i])*(mP[i] - mC1[i])) - ((mC1[i])*(mC1[i])));
		pC1x[i] = sqrt((((EC1[i])*(EC1[i]))/(c*c)) - ((c*mC1[i])*(c*mC1[i])));	//E^2 = p^2 + m^2 --> p = sqrt( E^2 - m^2 );
		pC1y[i] = 0.0;
		pC1z[i] = 0.0;
		//pC1 = c*sqrt( (mP-mC2)^2 - (mC1)^2 )
	}
	
	//RF(P): beta, gamma
	float beta[NEvents];
	float gamma[NEvents];
	
	for(int i = 0; i < NEvents; i++){
		
		beta[i] = (((pC1x[i])/(c))/(sqrt(((mC1[i])*(mC1[i]))+(((pC1x[i])*(pC1x[i]))/(c*c)))));
		gamma[i] = (sqrt((1.0)/(1.0 - ((beta[i])*(beta[i])))));
	
	}
	
	//RF(P) --> RF(C1)
	float L[NEvents][4][4];

	//initialize as zero
	for(int i = 0 ; i < NEvents; i++) {

		for(int j = 0; j < 4; j++) {

			for(int k = 0; k < 4; k++) {

				L[i][j][k] = 0.0;
			}
		}
	}

	//fill w/ lorentz transformation values
	for(int i = 0; i < NEvents; i++) {

		L[i][0][0] = gamma[i];
		L[i][1][1] = gamma[i];
		L[i][0][1] = (gamma[i])*(beta[i]);
		L[i][1][0] = (gamma[i])*(beta[i]);

		L[i][2][2] = 1.0;
		L[i][3][3] = 1.0;

	}

	//4-vectors
	float PC[NEvents][4];	//of C1 in RF(P)
	
	for(int i = 0; i < NEvents; i++){
		
		PC[i][0] = EC1[i];
		PC[i][1] = pC1x[i];
		PC[i][2] = pC1y[i];
		PC[i][3] = pC1z[i];
		
	}
	
	//REST FRAME: C1 ============================================================================================================================================
	
	//4-vector in RF(C1)
	
	float PCNew[NEvents][4];
	
	//initialize at zero
	for(int i = 0; i < NEvents; i++){
		
		for(int j = 0; j < 4; j++){
		
			PCNew[i][j] = 0.0;
			
		}
	}
	
	//perform Lorentz Boost
	
	for(int i = 0; i < NEvents; i++) {
		
		PCNew[i][0] = L[i][0][0]*PC[i][0] + L[i][0][1]*PC[i][1] + L[i][0][2]*PC[i][2] + L[i][0][3]*PC[i][3];
		PCNew[i][1] = L[i][1][0]*PC[i][0] + L[i][1][1]*PC[i][1] + L[i][1][2]*PC[i][2] + L[i][1][3]*PC[i][3];
		PCNew[i][2] = L[i][2][0]*PC[i][0] + L[i][2][1]*PC[i][1] + L[i][2][2]*PC[i][2] + L[i][2][3]*PC[i][3];
		PCNew[i][3] = L[i][3][0]*PC[i][0] + L[i][3][1]*PC[i][1] + L[i][3][2]*PC[i][2] + L[i][3][3]*PC[i][3];
	}
	
	
	//C1
	
	float ECNew[NEvents];	//new values for C1 in RF(C1)
	float mCNew[NEvents];

	float pCNewx[NEvents];
	float pCNewy[NEvents];
	float pCNewz[NEvents];

	float pNew[NEvents];

	//read values from 4-vectors
	for(int i = 0; i < NEvents; i++) {

		ECNew[i] = PCNew[i][0];
		pCNewx[i] = PCNew[i][1];
		pCNewy[i] = PCNew[i][2];
		pCNewz[i] = PCNew[i][3];

		pNew[i] = sqrt(((pCNewx[i])*(pCNewx[i]))+((pCNewy[i])*(pCNewy[i]))+((pCNewz[i])*(pCNewz[i])));	
		
		cout << pCNewx[i] << "..." << pCNewy[i] << "..." << pCNewz[i] << endl; 
	}
	
	//mCNew: mC in RF(C1)
	//E^2 = p^2c^2 + m^2c^4 --> m = sqrt((E^2/c^4) + (p^2/c^2);
	for(int i = 0 ; i < NEvents; i++) {
		
		mCNew[i] = sqrt((((ECNew[i])*(ECNew[i]))/(c*c*c*c))-(((pNew[i])*(pNew[i]))/(c*c)));
		

	}
	
	
	//G11, G12
	
	//e^+, e^- basically the same
	float mG[NEvents];
	float EG[NEvents];

	for(int i = 0; i < NEvents; i++) {

		mG[i] = me;	//electron & positron are VERY stable
		EG[i] = 0.5*ECNew[i];

	}
	
	//individual values of G11 and G12
	float pG[NEvents];
	float pG11x[NEvents];
	float pG11y[NEvents];
	float pG11z[NEvents];

	float pG21x[NEvents];
	float pG21y[NEvents];
	float pG21z[NEvents];

	float EG11[NEvents];
	float EG21[NEvents];
	float mG11[NEvents];
	float mG21[NEvents];
	
	float SignChance;
	float Sign;
	//pG, mG, EG arrays
	
	//test
	cout << "test line 349" << endl;
	for(int i = 0; i < NEvents; i++) {


		pG[i] = c*sqrt((0.25*(mCNew[i])*(mCNew[i]))-((me*me)));		//EG = 0.5 EC1' (+) E^2 = p^2c^2 + m^2c^4 [G] (+) E' = m'c^2 [C1]
		//cout << (0.25*(mCNew[i])*(mCNew[i]))-((me*me)) << endl;
		
		
		
		
		pG11x[i] = 0.0;
		pG11z[i] = 0.0;
		pG21x[i] = 0.0;
		pG21z[i] = 0.0;
		
		//Sign of G11,G12 momentum
		SignChance = gRandom->Rndm();
		
		if(SignChance < 0.5){
			Sign = 1.0;
		}
		else if(SignChance > 0.5){
			Sign = -1.0;
		}
		else if(SignChance == 0.5){
			Sign = -1.0;
		}
			
		
		pG11y[i] = Sign*pG[i];	//all momentum in y-dir.
		pG21y[i] = -1.0*(pG11y[i]);	//conservation of momentum
	
		EG11[i] = EG[i];
		EG21[i] = EG[i];

		mG11[i] = mG[i];
		mG21[i] = mG[i];
	
	}
	
	//LCtoP: BOOST C1 --> P
	
	float LCtoP[NEvents][4][4];

		//initialize to 0;
		for(int i = 0; i < NEvents; i++) {
			
			for(int j = 0; j < 4; j++) {

				for(int k = 0; k < 4; k++) {
					
					LCtoP[i][j][k] = 0.0;
				}
			}
		}


		//now, beta(in RF(C1)) = -1.0*beta(in RF(P))
		for(int i = 0; i < NEvents; i++) {

			LCtoP[i][0][0] = gamma[i];
			LCtoP[i][1][1] = gamma[i];
			LCtoP[i][0][1] = -1.0*(gamma[i])*(beta[i]);
			LCtoP[i][1][0] = -1.0*(gamma[i])*(beta[i]);

			LCtoP[i][2][2] = 1.0;
			LCtoP[i][3][3] = 1.0;

		}
	
		//generate new pG11 in RF(P)
		
		float PG11[NEvents][4];		//4-vector of g11 in RF(C1)
		float PG21[NEvents][4];		//4-vector of g21 in RF(C1)
	
		float PG11New[NEvents][4];	//4-vector of g11 in RF(P)
		float PG21New[NEvents][4];	//4-vector of g11 in RF(P)
		
		//pGNewSum: sum of e+, e- 4-vectors in RF(P)
		float PGNewSum[NEvents][4];
	
	
		//fill PG11 & PG21, initialize PG11New as 0
		for(int i = 0; i < NEvents; i++) {

			PG11[i][0] = EG11[i];
			PG11[i][1] = pG11x[i];
			PG11[i][2] = pG11y[i];
			PG11[i][3] = pG11z[i];

			PG21[i][0] = EG21[i];
			PG21[i][1] = pG21x[i];
			PG21[i][2] = pG21y[i];
			PG21[i][3] = pG21z[i];
		}

		for(int i = 0; i < NEvents; i++) {

			for(int j = 0; j <4; j++) {

				PG11New[i][j] = 0.0;
				PG21New[i][j] = 0.0;
			}
		}
	
		for(int i = 0; i < NEvents; i++) {

			PG11New[i][0] = LCtoP[i][0][0]*PG11[i][0] + LCtoP[i][0][1]*PG11[i][1] + LCtoP[i][0][2]*PG11[i][2] + LCtoP[i][0][3]*PG11[i][3];
			PG11New[i][1] = LCtoP[i][1][0]*PG11[i][0] + LCtoP[i][1][1]*PG11[i][1] + LCtoP[i][1][2]*PG11[i][2] + LCtoP[i][1][3]*PG11[i][3];
			PG11New[i][2] = LCtoP[i][2][0]*PG11[i][0] + LCtoP[i][2][1]*PG11[i][1] + LCtoP[i][2][2]*PG11[i][2] + LCtoP[i][2][3]*PG11[i][3];
			PG11New[i][3] = LCtoP[i][3][0]*PG11[i][0] + LCtoP[i][3][1]*PG11[i][1] + LCtoP[i][3][2]*PG11[i][2] + LCtoP[i][3][3]*PG11[i][3];

			
			PG21New[i][0] = LCtoP[i][0][0]*PG21[i][0] + LCtoP[i][0][1]*PG21[i][1] + LCtoP[i][0][2]*PG21[i][2] + LCtoP[i][0][3]*PG21[i][3];
			PG21New[i][0] = LCtoP[i][0][0]*PG21[i][0] + LCtoP[i][0][1]*PG21[i][1] + LCtoP[i][0][2]*PG21[i][2] + LCtoP[i][0][3]*PG21[i][3];
			PG21New[i][1] = LCtoP[i][1][0]*PG21[i][0] + LCtoP[i][1][1]*PG21[i][1] + LCtoP[i][1][2]*PG21[i][2] + LCtoP[i][1][3]*PG21[i][3];
			PG21New[i][1] = LCtoP[i][1][0]*PG21[i][0] + LCtoP[i][1][1]*PG21[i][1] + LCtoP[i][1][2]*PG21[i][2] + LCtoP[i][1][3]*PG21[i][3];

		}
	
	//summing up 4-vectors
		
		for(int i = 0; i < NEvents; i++) {
		
			PGNewSum[i][0] = PG11New[i][0] + PG21New[i][0];
			PGNewSum[i][1] = PG11New[i][1] + PG21New[i][1];
			PGNewSum[i][2] = PG11New[i][2] + PG21New[i][2];
			PGNewSum[i][3] = PG11New[i][3] + PG21New[i][3];
		}
		
		float mTotal[NEvents];	//the money plot
		
		float pTotal[NEvents];
		
		for(int i = 0; i < NEvents; i++){
		
			pTotal[i] = sqrt(((PGNewSum[i][1])*(PGNewSum[i][1]))+((PGNewSum[i][2])*(PGNewSum[i][2]))+((PGNewSum[i][3])*(PGNewSum[i][3])));
		}
		
		for(int i = 0; i < NEvents; i++){
		
			mTotal[i] = sqrt((((PGNewSum[i][0])/(c*c))*((PGNewSum[i][0])/(c*c)))+(((pTotal[i])/(c))*((pTotal[i])/(c))));
			
		}

		
		//calculating Theta
		//u * v = uv*cos(theta) --> cos(theta) = (u*v)/(uv) --> theta = acos(...)
		//all in RF(P)

		float dot[NEvents]; //dot product pG11*pG21
		float pG11[NEvents]; //|pG11| = sqrt( pG11x*pG11x + ...)
		float pG21[NEvents]; //|pG21| = sqrt( pG21x*pG21x + ...)

		float cosine[NEvents];
		float theta[NEvents];

		for(int i = 0; i < NEvents; i++) {

			dot[i] = (((PG11New[i][1])*(PG21New[i][1]))+((PG11New[i][2])*(PG21New[i][2]))+((PG11New[i][3])*(PG21New[i][3])));
			pG11[i] = sqrt(((PG11New[i][1])*(PG11New[i][1]))+((PG11New[i][2])*(PG11New[i][2]))+((PG11New[i][3])*(PG11New[i][3])));
			pG21[i] = sqrt(((PG21New[i][1])*(PG21New[i][1]))+((PG21New[i][2])*(PG21New[i][2]))+((PG21New[i][3])*(PG11New[i][3])));

			cosine[i] = (dot[i])/(pG11[i]*pG21[i]);
			theta[i] = acos(cosine[i]);	//in radians
			theta[i] = (theta[i])*((180.0)/(TMath::Pi()));	//conversion: radians -- degrees
		}

	//VISUALIZATION========================================================================================================================
	
	//canvas: split based on RF
		TCanvas *cP = new TCanvas("cP","REST FRAME:PARENT", 1500, 1500);
		TCanvas *cC1 = new TCanvas("cC1", "REST FRAME:C1", 1500, 1500);
		
		TCanvas *c1 = new TCanvas("c1","Important Plots; RF(P)", 1500, 1500);

		cP->Divide(3,2);
		cC1->Divide(3,2);
		c1->Divide(2,1);
	
	TH1D *hmP = new TH1D("hmP","^{8}Be^{*} Mass", NBins, (*min_element(mP,mP+NEvents)) - 10.0, (*max_element(mP,mP+NEvents)) + 10.0 );
	TH1D *hEP = new TH1D("hEP","^{8}Be^{*} Energy", NBins, (*min_element(EP,EP+NEvents)) - 10.0, (*max_element(EP,EP+NEvents)) + 10.0 );
	TH1D *hmC1 = new TH1D("hmC1","X17 Mass", NBins, (*min_element(mC1,mC1+NEvents)) - 10.0, (*max_element(mC1,mC1+NEvents)) + 10.0 );
	//TH1D *hmC2 = new TH1D("hmC1","C2 Mass Distribution", NBins, (*min_element(vmC2.begin(), vmC1.end())) - 5.0; (*max_element(vmC2.begin(), vmC1.end())) + 5.0);
	TH1D *hEC1 = new TH1D("hEC","X17 Energy", NBins, (*min_element(EC1,EC1+NEvents)) - 10.0, (*max_element(EC1,EC1+NEvents)) + 10.0 );
	TH1D *hpC1 = new TH1D("hpC1","X17 Momentum", NBins, (*min_element(pC1x,pC1x+NEvents)) - 10.0, (*max_element(pC1x,pC1x+NEvents)) + 10.0 );
	//TH1D *hpC2 = new TH1D("hpC2","Parent Mass Distribution", NBins, mPmin - 5.0; mPmax + 5.0);
	//

	TH1D *hmCNew = new TH1D("hmCNew","X17 Mass (RF:X17)", NBins,(*min_element(mCNew,mCNew+NEvents)) - 10.0, (*max_element(mCNew,mCNew+NEvents)) + 10.0 );
	TH1D *hECNew = new TH1D("hECNew","X17 Energy (RF:C1)", NBins, -10.0, (*max_element(ECNew,ECNew+NEvents)) + 10.0 );
	TH1D *hmG11 = new TH1D("hmG11","e^{+},e^{-} Mass", NBins, -10.0, (*max_element(mG11,mG11+NEvents)) + 10.0 );
	//TH1D *hmC2 = new TH1D("hmC1","C2 Mass", NBins, (*min_element(vmC2.begin(), vmC1.end())) - 5.0; (*max_element(vmC2.begin(), vmC1.end())) + 5.0);
	TH1D *hEG11 = new TH1D("hEG11","e^{+},e^{-} Energy", NBins, -10.0, (*max_element(EG11,EG11+NEvents)) + 10.0 );
	TH1D *hpG11 = new TH1D("hpG11","e^{+},e^{-} Momentum", NBins, (*min_element(pG11y,pG11y+NEvents))-10.0, (*max_element(pG11y,pG11y+NEvents)) + 10.0 );
	//TH1D *hpC2 = new TH1D("hpC2","Parent Mass Distribution", NBins, mPmin - 5.0; mPmax + 5.0);

	TH1D *hTheta = new TH1D("hTheta", "Angular Deflection", NBins, -5.0, 185.0);
		
	//determine min,max values of me+e-
	TH1D *hmTotal = new TH1D("hmTotal", "Invariant Mass Sum", NBins, (*min_element(mTotal,mTotal+NEvents)) - 10.0, (*max_element(mTotal,mTotal+NEvents)) + 10.0 );
	
	//filling histograms
	
	for(int i = 0 ; i < NEvents; i++) {

		hmP->Fill(mP[i]);
		hEP->Fill(EP[i]);
		hmC1->Fill(mC1[i]);
		hEC1->Fill(EC1[i]);
		hpC1->Fill(pC1x[i]);


		hmCNew->Fill(mCNew[i]);
		hECNew->Fill(ECNew[i]);
		hmG11->Fill(mG11[i]);
		hEG11->Fill(EG11[i]);
		hpG11->Fill(pG11y[i]);

		hTheta->Fill(theta[i]);
		hmTotal->Fill(mTotal[i]);

	}

	hmP->GetXaxis()->SetTitle("Mass [\frac{MeV}{c^2}]");
	hmP->GetYaxis()->SetTitle("Counts");
	hEP->GetXaxis()->SetTitle("Energy [MeV]");
	hEP->GetYaxis()->SetTitle("Counts");
	hmC1->GetXaxis()->SetTitle("Mass [\frac{MeV}{c^2}]");
	hmC1->GetYaxis()->SetTitle("Counts");
	hEC1->GetXaxis()->SetTitle("Energy [MeV]");
	hEC1->GetYaxis()->SetTitle("Counts");
	hpC1->GetXaxis()->SetTitle("Momentum [\frac{MeV}{c}]");
	hpC1->GetYaxis()->SetTitle("Counts");

	hmCNew->GetXaxis()->SetTitle("Mass [\frac{MeV}{c^2}]");
	hmCNew->GetYaxis()->SetTitle("Counts");
	hECNew->GetXaxis()->SetTitle("Energy [MeV]");
	hECNew->GetYaxis()->SetTitle("Counts");
	hmG11->GetXaxis()->SetTitle("Mass [\frac{MeV}{c^2}]");
	hmG11->GetYaxis()->SetTitle("Counts");
	hEG11->GetXaxis()->SetTitle("Energy [MeV]");
	hEG11->GetYaxis()->SetTitle("Counts");
	hpG11->GetXaxis()->SetTitle("Momentum [\frac{MeV}{c}]");
	hpG11->GetYaxis()->SetTitle("Counts");

	hTheta->GetXaxis()->SetTitle("#Theta [^{#circ}]");
	hTheta->GetYaxis()->SetTitle("Counts");
		
	hmTotal->GetXaxis()->SetTitle("m_{e^{+}e^{-}} [\frac{MeV}{c^2}]");
	hmTotal->GetYaxis()->SetTitle("Counts");
	
	hmP->SetFillColor(0);
	hmP->SetLineWidth(2);
	hEP->SetFillColor(1);
	hmP->SetLineWidth(2);
	hmC1->SetFillColor(2);
	hmC1->SetLineWidth(2);
	hEC1->SetFillColor(3);
	hEC1->SetLineWidth(2);
	hpC1->SetFillColor(4);
	hpC1->SetLineWidth(2);

	hmCNew->SetFillColor(5);
	hmCNew->SetLineWidth(2);
	hECNew->SetFillColor(6);
	hECNew->SetLineWidth(2);
	hmG11->SetFillColor(7);
	hmG11->SetLineWidth(2);
	hEG11->SetFillColor(8);
	hEG11->SetLineWidth(2);
	hpG11->SetFillColor(9);
	hpG11->SetLineWidth(2);

	hTheta->SetFillColor(10);
	hTheta->SetLineWidth(2);
		
	hmTotal->SetFillColor(11);
	hmTotal->SetLineWidth(2);

	cP->cd(1);
	hmP->Draw();

	cP->cd(2);
	hEP->Draw();
	
	cP->cd(3);
	hmC1->Draw();
	
	cP->cd(4);
	hEC1->Draw();
	
	cP->cd(5);
	hpC1->Draw();

	c1->cd(1);
	hTheta->Draw();
		
	c1->cd(2);
	hmTotal->Draw();
	
	cC1->cd(1);
	hmCNew->Draw();
	
	cC1->cd(2);
	hECNew->Draw();

	cC1->cd(3);
	hmG11->Draw();

	cC1->cd(4);
	hEG11->Draw();

	cC1->cd(5);
	hpG11->Draw();

	

	cP->Draw();
	cC1->Draw();
	c1->Draw();


	

}
