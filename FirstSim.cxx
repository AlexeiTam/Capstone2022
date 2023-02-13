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
	float EC1[NEvents];
	float EC2[NEvents];
		//beta-->constrained: [0,1]	//for simplicity
	for(int i = 0; i < NEvents; i++) {

		betaC[i] = gRandom->Rndm();
		gammaC[i] = sqrt((1.0)/(1.0 - (betaC[i])*(betaC[i])));
		mC[i] = gRandom->Gaus(0.5*mPMean, 0.5*mPSigma);	
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
		//(?)p = gamma*mv = gamma*mc*beta
		//!!ED: generate p, let E be determined by constraint
	for(int i = 0; i < NEvents; i++) {

		pC1z[i] = gRandom->Rndm();
		pC2z[i] = -1.0*pC1z[i];
	}

		//constraint: E(i)^2 = p(i)^2c^2 + m(i)^2c^4

	for(int i = 0; i < NEvents; i++) {

		EC1[i] = sqrt((pC1z[i])*(pC1z[i])*(c*c) + (mC[i])*(mC[i])*(c*c*c*c) );	//|p| = pz ; as px = py = 0
		EC2[i] = sqrt((pC2z[i])*(pC2z[i])*(c*c) + (mC[i])*(mC[i])*(c*c*c*c) );	//|p| = pz ; as px = py = 0
	}


		//Generate L[NEvents][rows][columns]
		//
	

}
