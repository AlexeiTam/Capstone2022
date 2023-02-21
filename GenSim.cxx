void GenSim(){
	//GOAL: simulate general two-step decay: P->C1 + C2; C1->g11 + g12
	//ED: Natural Units; i.e, c = 1; can change to c = 3e-8 later
	const Int_t NEvents = 10000; //number of simulations: 10,000
	float c = 1.0;

//REST FRAME: PARENT----------------------------------------------------------------------------------------
		//generate mP-->gives EP as EP = mP*c^2
		//generate mP by Breit-Wigner
		//generate mC,EP by BW (avg 0.25*P, width 0.25*(width of P))
		//generate pC:
			//choose x-dir. of RF(P) s.t. pC1y = pC1z = 0, |pC| = pCx
			//by conservation, pC1x = -1.0*pC2x
	
	std::cout << "REST FRAME: PARENT" << std::endl;
	std::cout << "generating mP..." << std::endl;

	//GENERATNG mP, EP, mC, EC
	//parameters, arrays
	float mP[NEvents];
	float EP[NEvents];
	float mC[NEvents];
	float EC[NEvents];
	float mC1[NEvents];
	float mC2[NEvents];
	float EC1[NEvents];
	float EC2[NEvents];

	float mPMean = 5.0;	//can be made initializable
	float mPSigma = 0.50;	//same as above
	float mCMean = 0.25*mPMean;
	float mCSigma = 0.25*mCSigma;


	//generation: mP, EP
	for(int i = 0; i < NEvents; i++) {

		//ED: gaussian instead of BW; work out later
		mP[i] = gRandom->Gaus(mPMean, mPSigma);
		EP[i] = (mP[i])*c*c;

		mC[i] = gRandom->Gaus(0.25*mPMean, 0.25*mPSigma);
		mC1[i] = mC[i];
		mC2[i] = mC[i];

		EC[i] = gRandom->Gaus(0.25*c*c*mPMean, 0.25*mPSigma*c*c);	//fix::get avg of EP?
		EC1[i] = EC[i];
		EC2[i] = EC[i];

			//if(i == 0.125*NEvents) std::cout << "12.5%" << std::endl;
			//if(i == 0.250*NEvents) std::cout << "25.0%" << std::endl;
			//if(i == 0.375*NEvents) std::cout << "37.5%" << std::endl;
			//if(i == 0.500*NEvents) std::cout << "50.0%" << std::endl;
			//if(i == 0.625*NEvents) std::cout << "62.5%" << std::endl;
			//if(i == 0.750*NEvents) std::cout << "75.0%" << std::endl;
			//if(i == 0.875*NEvents) std::cout << "87.5%" << std::endl;
			//if(i == NEvents-1) std::cout << "PARENT MASS GENERATION COMPLETE" << std::endl;
	}

	//generation: mC, EC, pC1, pC2
	//??p-->beta?
	//

	float beta[NEvents];
	float gamma[NEvents];	//related to pC in RF(P)

	float pC1x[NEvents];
	float pC1y[NEvents];
	float pC1z[NEvents];

	float pC2x[NEvents];
	float pC2y[NEvents];
	float pC2z[NEvents];

	float chance; 	//flip a coin; add chance that pC1 can be negative

	//generate pC: y and z = 0; x constrained by E^2 = p^2c^2 + m^2c^4
	for(int i = 0; i < NEvents; i++) {

		chance = gRandom->Rndm();
		pC1x[i] = sqrt(((EC[i]*EC[i])/(c*c))+(mC[i]*mC[i]*c*c));
		if(chance < 0.5) pC1x[i] = 1.0*pC1x[i];

		pC2x[i] = -1.0*pC1x[i];

		pC1y[i] = 0.0;
		pC1z[i] = 0.0;
		pC2y[i] = 0.0;
		pC2z[i] = 0.0;
	
	}


	//generate beta, gamma
	for(int i = 0; i < NEvents; i++) {

		beta[i] = ((pC1x[i])/(c))/(sqrt((mC1[i]*mC1[i])+((pC1x[i]*pC1x[i])/(c*c))));
		gamma[i] = sqrt((1.0)/(1.0 - (beta[i]*beta[i])));

	}

	//generate Boost: L[Event][x][y][z]
	//initialize all as zero
	//then, copy/paste from prev. code
	
	float L[NEvents][4][4];

	for(int i = 0 ; i < NEvents; i++) {

		for(int j = 0; j < 4; j++) {

			for(int k = 0; k < 4; k++) {

				L[i][j][k] = 0.0;
			}
		}
	}


	//??Use reverse lorentz transformation?
	//generate based on gamma, beta
	for(int i = 0; i < NEvents; i++) {

		L[i][0][0] = gamma[i];
		L[i][1][1] = gamma[i];
		L[i][0][1] = (gamma[i])*(beta[i]);
		L[i][1][0] = (gamma[i])*(beta[i]);

		L[i][2][2] = 1.0;
		L[i][3][3] = 1.0;

	}

	//REST FRAME: C1------------------------------------------------------------------------------------------------
	
	//BOOKMARK: perform boost on PC	

	float PC[NEvents][4];	//4-vector of C1 in RF(P)
	float PCNew[NEvents][4];	//4-vector of C1 in RF(C1)

	for(int i = 0; i < NEvents; i++) {

		PC[i][0] = EC1[i];
		PC[i][1] = pC1x[i];
		PC[i][2] = pC1y[i];
		PC[i][3] = pC1z[i];
		
	}

	for(int i = 0; i < NEvents; i++) {


		for(int j = 0; j < 4; j++) {

			PCNew[i][j] = 0.0;
		}
	}


	for(int i = 0; i < NEvents; i++) {
		
		for(int j = 0; j < 4; j++) {

			for(int k = 0; k < 4; k++) {

				PCNew[i][j] = PCNew[i][j] + (L[i][j][k])*(PC[i][k]); 
			}
		}

	}

	//fill ECnew, generate mCnew: from ECnew^2 = EC^2 = mCnew*c^2
	
	float ECNew[NEvents];
	float mCNew[NEvents];

	float pCNewx[NEvents];
	float pCNewy[NEvents];
	float pCNewz[NEvents];

	float pNew[NEvents];

	for(int i = 0; i < NEvents; i++) {

		PCNew[i][0] = ECNew[i];
		PCNew[i][1] = pCNewx[i];
		PCNew[i][2] = pCNewy[i];
		PCNew[i][3] = pCNewz[i];

		pNew[i] = sqrt(((pCNewx[i])*(pCNewx[i]))+((pCNewy[i])*(pCNewy[i]))+((pCNewz[i])*(pCNewz[i])));	
	}

	//find mCNew: mC in RF(C1)
	//E^2 = p^2c^2 + m^2c^4 --> m = sqrt((E^2/c^4) + (p^2/c^2);
	for(int i = 0 ; i < NEvents; i++) {
		
		mCNew[i] = sqrt((((ECNew[i])*(ECNew[i]))/(c*c*c*c))+(((pNew[i])*(pNew[i]))/(c*c)));

	}

//old testing
//cout << "pC1x:" << pC1x[0] << "..." << "pCNewx:" << pCNewx[0] << endl;
//cout << "pC1y:" << pC1y[0] << "..." << "pCNewy:" << pCNewy[0] << endl;
//cout << "pC1z:" << pC1z[0] << "..." << "pCNewz:" << pCNewz[0] << endl;

	//C1-->g11 + g21
	
	//mG = mG11 = mG12
	
		//first, determine avg, width of mCNew
		float mCNewMean = 0.0;
		float mCNewSigma = 0.0;

		float ECNewMean;
		float ECNewSigma;

		for(int i = 0; i < NEvents; i++) {

			mCNewMean = mCNewMean + (mCNew[i])/NEvents;
			ECNewMean = ECNewMean + (ECNew[i])/NEvents;
		}

		for(int i = 0; i < NEvents; i++) {

			mCNewSigma = mCNewSigma + ((mCNew[i] - mCNewMean)*(mCNew[i] - mCNewMean))/(NEvents - 1);
			ECNewSigma = ECNewSigma + ((ECNew[i] - ECNewMean)*(ECNew[i] - ECNewMean))/(NEvents - 1);
		}

		mCNewSigma = sqrt(mCNewSigma);
		ECNewSigma = sqrt(ECNewSigma);
	
	//generating mG::	
	float mG[NEvents];

	for(int i = 0; i < NEvents; i++) {

		mG[i] = gRandom->Gaus(mCNewMean, mCNewSigma);
		EG[i] = gRandom->Gaus(ECNewMean, ECNewSigma);

	}

	//generating pG: 
	//constrained: EG^2 = pG^2c^2 + mG^2c^4
	//pG11 = - pG21 --> pG11x = - pG21x; same for y and z
	//??is this right?
		//pick y s.t. pG || y in RF(C1)
		//then pG11x = pG21 x = pG11z = pG21z = 0.0
	
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
	
	//filling pG, individual EG, mG arrays
	for(int i = 0; i < NEvents; i++) {


		pG[i] = sqrt((((EG[i])*(EG[i]))/(c*c))+((mG[i])*(mG[i])*c*c));
	
		EG11[i] = EG[i];
		EG21[i] = EG[i];

		mG11[i] = mG[i];
		mG21[i] = mG[i];
	
	}

	//filling individual pG's
	
	for(int i = 0; i < NEvents; i++) {

		pG11x[i] = 0.0;
		pG11z[i] = 0.0;
		pG21x[i] = 0.0;
		pG21z[i] = 0.0;

	}

	for(int i = 0; i < NEvents; i++) {

		pG11y[i] = pG[i];
		pG21y[i] = -1.0*(pG11y[i]);

	}

	//BOOST: RF(C1) --> RF(P) --------------------------------------------------------------------------------------
	

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
			LCtoP[i][0][1] = (gamma[i])*(beta[i]);
			LCtoP[i][1][0] = (gamma[i])*(beta[i]);

			LCtoP[i][2][2] = 1.0;
			LCtoP[i][3][3] = 1.0;

		}


		//generate new pG11 in RF(P)
		
		float PG11[NEvents][4];		//4-vector of g11 in RF(C1)
		float PG11New[NEvents][4];	//4-vector of g11 in RF(P)

		//fill pG11, initialize pG11New as 0
		for(int i = 0; i < NEvents; i++) {

			PG11[i][0] = EG11[i];
			PG11[i][1] = pG11x[i];
			PG11[i][2] = pG11y[i];
			PG11[i][3] = pG11z[i];

		}

		for(int i = 0; i < NEvents; i++) {

			for(int j = 0; j <4; j++) {

				PG11New[i][j] = 0.0;
			}
		}

		//PG11New = LCtoP*PG11
		for(int i = 0; i < NEvents; i++) {

			for(int j = 0; j < 4; j++) {

				for(int k = 0; k < 4; k++) {

					PG11New[i][j] = PG11New[i][j] + (LCtoP[i][j][k])*(PG11[i][k]);
				}
			}
		}






}
