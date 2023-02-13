void FirstSim() {
	//!! ED: c=1;
	const Int_t NEvents = 10000;	//Number of simulations
	float c = 1.0;
//PARENT

	//Generating mP -- Mass of Parent
	//then generate beta-->Calculate gamma = 1/sqrt(1-beta^2)
	float mPMean = 5.0;
	float mPSigma = 0.50;
	float mP[NEvents];
	float EP[NEvents];	//constrained:E^2 = p^2c^2 + mP^2c^4


		
		//
		//beta-->constrained: [0,1]	//for simplicity
	for(int i = 0; i < NEvents; i++) {

		mP[i] = gRandom->Gaus(mPMean, mPSigma);
		EP[i] = sqrt((mP[i])*(mP[i])*c*c*c*c);

	}

//P-->C1 + C2
	
	float betaC[NEvents];	//betaC1 = -1.0*betaC2	//RF(P): p(c1) = - p(c2)	
	float gammaC[NEvents];	//gammaC1 = gammaC2, for same reasons as above
	float pC1x[NEvents];
	float pC1y[NEvents];
	float pC1z[NEvents];
	
	float pC2x[NEvents];
	float pC2y[NEvents];
	float pC2z[NEvents];

	float mC[NEvents];	//mC1 = mC2	//for simplicity
	
		//beta-->constrained: [0,1]	//for simplicity
	for(int i = 0; i < NEvents; i++) {

		betaC[i] = gRandom->Rndm();
		gammaC[i] = sqrt((1.0)/(1.0 - (betaC[i])*(betaC[i])));

	}

		//generating pC1, pC2 with pC1 = -pC2
		//start: all null
	for(int i = 0; i < NEvents; i++) {

		pC1x[i] = 0.0;
		pC1y[i] = 0.0;
		pC1z[i] = 0.0;

		pC2x[i] = 0.0;
		pC2y[i] = 0.0;
		pC2z[i] = 0.0;
	}

		//CHOOSE z-direction s.t. pC1 and PC2 are purely in z-axis
	for(int i = 0; i < NEvents; i++) {

		pC1z[i] = 
	}


}
