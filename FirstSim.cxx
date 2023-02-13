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

	//!!Move mC generation here
	//
		//beta-->constrained: [0,1]	//for simplicity
	for(int i = 0; i < NEvents; i++) {

		betaC[i] = gRandom->Rndm();	//fix: constrained by E^2 = p^2 + m^2
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

		//CHOOSE x-direction s.t. pC1 and PC2 are purely in x-axis
		//(?)p = gamma*mv = gamma*mc*beta
		//!!ED: generate p, let E be determined by constraint
	for(int i = 0; i < NEvents; i++) {

		pC1x[i] = gRandom->Rndm();//!!E^2 = p^2 + m^2
		pC2x[i] = -1.0*pC1x[i];
	}

		//constraint: E(i)^2 = p(i)^2c^2 + m(i)^2c^4

	for(int i = 0; i < NEvents; i++) {

		EC1[i] = sqrt((pC1x[i])*(pC1x[i])*(c*c) + (mC[i])*(mC[i])*(c*c*c*c) );	//|p| = pz ; as px = py = 0
		EC2[i] = sqrt((pC2x[i])*(pC2x[i])*(c*c) + (mC[i])*(mC[i])*(c*c*c*c) );	//|p| = pz ; as px = py = 0
	}


		//Generate L[NEvents][rows][columns]
		//and momentum 4-vectors P[NEvents][rows]
	float L[NEvents][4][4];
	float PC1[NEvents][4];
	float PC2[NEvents][4];

	//initialize as 0 matrix
	for(int i = 0; i < NEvents; i++) { //Event Scan

		for(int j = 0; j < 4; j++) {

			for(int k = 0; k < 4; k++) {

				L[i][j][k] = 0.0;
			}
		}

	}	


	//fill with appropriate gamma, beta values for each event
	//
	for(int i = 0; i < NEvents; i++) {

		L[i][0][0] = gammaC[i];
		L[i][1][1] = gammaC[i];
		L[i][0][1] = -1.0*(gammaC[i])*(betaC[i]);
		L[i][1][0] = -1.0*(gammaC[i])*(betaC[i]);

		L[i][2][2] = 1.0;
		L[i][3][3] = 1.0;
	}

	//fill P-vectors:
	//
	for(int i = 0; i < NEvents; i++) {

		PC1[i][0] = EC1[i];
		PC1[i][1] = pC1x[i];
		PC1[i][2] = pC1y[i];
		PC1[i][3] = pC1z[i];

		PC2[i][0] = EC2[i];
		PC2[i][1] = pC2x[i];
		PC2[i][2] = pC2y[i];
		PC2[i][3] = pC2z[i];
	}

//BOOST: P->C1

	float PCprime[NEvents][4];	//PCprime = PC1 in RF(C1)
	
	//initialize as 0
	for(int i = 0; i < NEvents; i++) {

		for(int j = 0; j < 4; j++) {

			PCprime[i][j] = 0.0;
		}
	}

	//perform LT:
	for(int k = 0; k < NEvents; k++) {
		
		for(int i = 0 ; i < 4; i ++) {

			for(int j = 0; j < 4; j++) {

				PCprime[k][i] = PCprime[k][i] + (L[k][i][j])*PC1[k][j];
			}
		}
		
	}



//RF(C1):
//

	//CHECK::
	//
	const Int_t NBins = 1000;
	float xmin = mPMean - (5.0*mPSigma); 
	float xmax = mPMean + (5.0*mPSigma); 
	TH1D *hmP = new TH1D("hmP","mP Distribution", NBins, xmin, xmax);

		

}
