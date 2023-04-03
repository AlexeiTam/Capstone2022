void Run2GenSim(){
	
	//SAME AS GenSim.cxx, but with improvements (Run2Improvements.txt)
	//GOAL: simulate general two-step decay: P->C1 + C2; C1->g11 + g12

	//!!Everything is working on MeV
	//mass = MeV/c^2
	//momentum = MeV/c
	
	//ED: Natural Units; i.e, c = 1; can change to c = 3e-8 later
	const Int_t NEvents = 1000; //number of simulations: 10,000
	float c = 1.0;
	float pi = TMath::Pi();

//REST FRAME: PARENT----------------------------------------------------------------------------------------
		//generate mP-->gives EP as EP = mP*c^2
		//generate mP by Breit-Wigner
		//generate EP by BW rest mass-energy: E = mc^2
		//generate pC:
			//choose x-dir. of RF(P) s.t. pC1y = pC1z = 0, |pC| = pCx
			//by conservation, pC1x = -1.0*pC2x
	
	std::cout << "REST FRAME: PARENT" << std::endl;
	std::cout << "generating Parent, C1 masses & energies..." << std::endl;

	//GENERATNG mP, EP, mC, EC
	//parameters, arrays
	std::vector<float> mP;	//float mP[NEvents];
	std::vector<float> EP;	//float EP[NEvents];
	std::vector<float> mC;	//float mC[NEvents];
	std::vector<float> EC;	//float EC[NEvents];
	std::vector<float> mC1;	//float mC1[NEvents];
	std::vector<float> mC2;	//float mC2[NEvents];
	std::vector<float> EC1;	//float EC1[NEvents];
	std::vector<float> EC2;	//float EC2[NEvents];

	//P = Be*
	float mPMean = 8022.2851;	//?? issues w/ m(Be*) determined by p(X17), random value; rn using m(Be) + m(X17)
	float mPSigma = 3.9E-11;	//Be8 is very stable
	
	//C2 = Be
	float mC2Mean = 8005.30510;	//base value
	float mC2Sigma = 3.9E-11;	//Be8 is very stable

	//C1 = X17
	float mC1Mean = 16.98;	//X17 ==> invariant mass about 17 MeV
	float mC1Sigma = 1.0;	//will play with later?

	//G11, G12 = e^+, e^-
	
	float me = 0.511;	//mass of electron = 0.511 MeV
	//float meSigma = 0.0;	//very stable, no spread

	//float mCMean = 0.25*mPMean;
	//float mCSigma = 0.25*mCSigma;

	//float fmP, fEP, fmC, fEC, fmC1, fmC2, fEC1, fEC2;
	//generation: mP, EP

//Breit-Wigner Dist.

	//more parameters
	float norm = (2*mPSigma*mPSigma)/((fabs(mPSigma))*(fabs(mPSigma))*(fabs(mPSigma)));	
	//normalization const.
	
	float cutoff = 0.10;
	float mPMax = mPMean + (sqrt(((norm*mPSigma)/(pi*cutoff))-(mPSigma*mPSigma))); 
	//not taking m with P < cutoff
	
	TF1 *BW = new TF1("BW","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mPMax);
	BW->SetParameters(pi, mPMean, mPSigma, norm);

	//generating P:
	
	float chance, prob;
	int index = 0;

	float mPplaceholder;	//placeholder
	while(index < NEvents + 1) {

		index++;
		
		mPplaceholder = mPMax*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax

		prob = gRandom->Rndm();		//choose some #, associated with x
		chance = BW->Eval(mPplaceholder);		//threshold for x to be accepted

		if(prob > chance) {
			index = index -1;
		}
		//prob > chance --> reject

		else if(prob < chance) {
			mP.emplace_back(mPplaceholder);
		}
		//prob < chance --> accept
		
		//do it again, until we get all the parents we need
	
	}



//GENERATING C1, C2
	
	for(int i = 0; i < NEvents; i++){
	
		//C2 mass & energy
		
		float mC2chance, mC2prob;
		int mC2index = 0;
		float mC2norm = (2*mC2Sigma*mC2Sigma)/((fabs(mC2Sigma))*(fabs(mC2Sigma))*(fabs(mC2Sigma)));
		float mC2Max = mC2Mean + (sqrt(((mC2norm*mC2Sigma)/(pi*cutoff))-(mC2Sigma*mC2Sigma)));

	float mC2placeholder;	//placeholder
	while(mC2index < NEvents + 1) {

		mC2index++;
		
		mC2placeholder = mC2Max*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax

		prob = gRandom->Rndm();		//choose some #, associated with x
		chance = BW->Eval(mC2placeholder);		//threshold for x to be accepted

		if(prob > chance) {
			index = index -1;
		}
		//prob > chance --> reject

		else if(prob < chance) {
			mC2.emplace_back(mC2placeholder);
		}
		//prob < chance --> accept
		
		EC2.emplace_back((mC2.at(i))*c*c);	//P: at rest --> C2, approx. at rest
		
		//C1: mass & energy
		
		EC1.emplace_back((EP.at(i))-(EC2.at(i)));	//EP = EC1 + EC2 (Energy conservation) --> EC1 = EP - EC2
	
		float mC1chance, mC1prob;
		int mC1index = 0;
		float mC1norm = (2*mC1Sigma*mC1Sigma)/((fabs(mC1Sigma))*(fabs(mC1Sigma))*(fabs(mC1Sigma)));
		float mC1Max = mC1Mean + (sqrt(((mC1norm*mC1Sigma)/(pi*cutoff))-(mC1Sigma*mC1Sigma)));

	float mC1placeholder;	//placeholder
	while(mC1index < NEvents + 1) {

		mC1index++;
		
		mC1placeholder = mC1Max*(gRandom->Rndm());	//select a candidate for mP, from 0 < mP < mPMax

		prob = gRandom->Rndm();		//choose some #, associated with x
		chance = BW->Eval(mC1placeholder);		//threshold for x to be accepted

		if(prob > chance) {
			index = index -1;
		}
		//prob > chance --> reject

		else if(prob < chance) {
			mC1.emplace_back(mC1placeholder);
		}
		//prob < chance --> accept
	}
		
	//C1,C2: momenta:
		
	std::vector<float> beta;
	std::vector<float> gamma;

	std::vector<float> pC1x;
	std::vector<float> pC1y;
	std::vector<float> pC1z;
	
	std::vector<float> pC2x;
	std::vector<float> pC2y;
	std::vector<float> pC2z;	
		
	for(int i = 0; i < NEvents; i++){	
		//!!momenta only in x-axis
	pC1y.emplace_back(0.0);
	pC1z.emplace_back(0.0);
	pC2y.emplace_back(0.0);
	pC2z.emplace_back(0.0);
		
	pC2x.emplace_back(0.0);	//!! 8Be* @ rest --> 8Be approx. at rest
	pC1x.emplace_back(c*sqrt((((mP.at(i))-(mC2.at(i)))*((mP.at(i))-(mC2.at(i))))-((mC1.at(i))*(mC1.at(i)))));	//Energy conservation + Relativistic Eqn.
		//pC1 = c*sqrt((mP - mC2)^2 - (mC1)^2)
		
	}

	//generate beta, gamma
	
	for(int i = 0; i < NEvents; i++) {

		beta.emplace_back(((pC1x.at(i))/(c))/(sqrt(((mC1.at(i))*(mC1.at(i)))+(((pC1x.at(i))*(pC1x.at(i)))/(c*c)))));
		gamma.emplace_back(sqrt((1.0)/(1.0 - ((beta.at(i))*(beta.at(i))))));

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

		L[i][0][0] = gamma.at(i);
		L[i][1][1] = gamma.at(i);
		L[i][0][1] = (gamma.at(i))*(beta.at(i));
		L[i][1][0] = (gamma.at(i))*(beta.at(i));

		L[i][2][2] = 1.0;
		L[i][3][3] = 1.0;

	}

	//REST FRAME: C1------------------------------------------------------------------------------------------------
	
	//BOOKMARK: perform boost on PC	

	float PC[NEvents][4];
	float PCNew[NEvents][4];

	//float PC[NEvents][4];	//4-vector of C1 in RF(P)
	//float PCNew[NEvents][4];	//4-vector of C1 in RF(C1)
	for(int i = 0; i < NEvents; i++) {

		PC[i][0] = EC1.at(i);
		PC[i][1] = pC1x.at(i);
		PC[i][2] = pC1y.at(i);
		PC[i][3] = pC1z.at(i);
		
		//PC.push_back({EC1.at(i), pC1x.at(i), pC1y.at(i), pC1z.at(i)});
	}

	for(int i = 0; i < NEvents; i++) {

			for(int j = 0; j < 4; j++){
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

	//??May take out this loop

	//float placeholder2 = 0.0;

	//for(int i = 0; i < NEvents; i++) {
		
		//for(int j = 0; j < 4; j++) {

		//	placeholder2 = 0.0;
		//	for(int k = 0; k < 4; k++) {

		//		placeholder2 = placeholder2 + ((L[i][j][k])*(PC[i].at(k)))

		//	}
		//}

		
		
	}

	
	//fill ECnew, generate mCnew: from ECnew^2 = EC^2 = mCnew*c^2
	
	
	float ECNew[NEvents];
	float mCNew[NEvents];

	float pCNewx[NEvents];
	float pCNewy[NEvents];
	float pCNewz[NEvents];

	float pNew[NEvents];

	for(int i = 0; i < NEvents; i++) {

		ECNew.emplace_back(PCNew[i][0]);
		pCNewx[i] = PCNew[i][1];
		pCNewy[i] = PCNew[i][2];
		pCNewz[i] = PCNew[i][3];

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
		//float mCNewMean = 0.0;
		//float mCNewSigma = 0.0;

		//float ECNewMean;
		//float ECNewSigma;

		//for(int i = 0; i < NEvents; i++) {

		//	mCNewMean = mCNewMean + (mCNew[i])/NEvents;
		//	ECNewMean = ECNewMean + (ECNew[i])/NEvents;
		//}

		//for(int i = 0; i < NEvents; i++) {

			//mCNewSigma = mCNewSigma + ((mCNew[i] - mCNewMean)*(mCNew[i] - mCNewMean))/(NEvents - 1);
			//ECNewSigma = ECNewSigma + ((ECNew[i] - ECNewMean)*(ECNew[i] - ECNewMean))/(NEvents - 1);
		//}

		//mCNewSigma = sqrt(mCNewSigma);
		//ECNewSigma = sqrt(ECNewSigma);
	
	//generating mG::	
	float mG[NEvents];
	float EG[NEvents];

	for(int i = 0; i < NEvents; i++) {

		mG[i] = meMean;	//electron & positron are VERY stable
		EG[i] = 0.5*ECNew.at(i);

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


		pG[i] = c*sqrt((0.25*(mCNew[i])*(mCNew[i]))-((me*me)));
	
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
		

		//fill pG11, initialize pG11New as 0
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

		//PG11New = LCtoP*PG11
		//for(int i = 0; i < NEvents; i++) {

			//for(int j = 0; j < 4; j++) {

			//	for(int k = 0; k < 4; k++) {

			//		PG11New[i][j] = PG11New[i][j] + (LCtoP[i][j][k])*(PG11[i][k]);
			//		PG21New[i][j] = PG21New[i][j] + (LCtoP[i][j][k])*(PG21[i][k]);
			//	}
			//}
		//}

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
		std::vector<float> vmTotal;
		float pTotal[NEvents];
		
		for(int i = 0; i < NEvents; i++){
		
			pTotal[i] = sqrt(((pGNewSum[i][1])*(pGNewSum[i][1]))+((pGNewSum[i][1])*(pGNewSum[i][1]))+((pGNewSum[i][1])*(pGNewSum[i][1])));
		}
		
		for(int i = 0; i < NEvents; i++){
		
			mTotal[i] = sqrt((((PGNewSum[i][0])/(c*c))*((PGNewSum[i][0])/(c*c)))-(((pTotal[i])/(c))*((pTotal[i])/(c))));
			
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

		//test

		cout << "RF:P" << endl;
		//cout << "i...pG11[i]" << endl;
		cout << "i...pg11x[i]...pg11y[i]...pg11z[i]...cosine[i]...theta[i]" << endl;
		for(int i = 0; i < NEvents; i++) {

			//cout << i << "..." << pG11[i] << endl;
			cout << i << "..." << PG11New[i][1] << "..." << PG11New[i][2] << "..." << PG11New[i][3] << "..." << cosine[i] << "..." << theta[i] << endl;

		}
	//VISUALIZATION :D  --------------------------------------------------------------------------

		//note to self; maybe just use vectors in the first place
	
	const Int_t NBins = 300;	//
	std::vector<float> vmP;
	std::vector<float> vEP;
	std::vector<float> vmC1;
	std::vector<float> vmC2;
	std::vector<float> vEC1;
	std::vector<float> vpC1;
	std::vector<float> vpC2;

	std::vector<float> vmCNew;
	std::vector<float> vECNew;
	std::vector<float> vmG11;
	std::vector<float> vEG11;
	std::vector<float> vmG21;
	std::vector<float> vpG11;
	std::vector<float> vpG21;

	std::vector<float> vEG11New;
	std::vector<float> vmG11New;
	std::vector<float> vpG11xNew;
	std::vector<float> vpG11yNew;
	std::vector<float> vpG11zNew;
		
	std::vector<float> vmTotal;	//the money plot :)

		//canvas: split based on RF
		TCanvas *cP = new TCanvas("cP","REST FRAME:PARENT", 1500, 1500);
		TCanvas *cC1 = new TCanvas("cC1", "REST FRAME:C1", 1500, 1500);
		
		TCanvas *c1 = new TCanvas("c1","The Money Plot; RF(P)", 1500, 1500);

		cP->Divide(3,2);
		cC1->Divide(3,2);

	for(int i = 0; i < NEvents; i++) {

		vmP.emplace_back(mP[i]);
		vEP.emplace_back(EP[i]);
		vmC1.emplace_back(mC1[i]);
		vmC2.emplace_back(mC2[i]);
		vEC1.emplace_back(EC1[i]);
		vpC1.emplace_back(pC1x[i]);
		vpC2.emplace_back(pC2x[i]);

		vmCNew.emplace_back(mCNew[i]);
		vECNew.emplace_back(EC[i]);
		vmG11.emplace_back(mG11[i]);
		vEG11.emplace_back(EG11[i]);
		vmG21.emplace_back(mG21[i]);
		vpG11.emplace_back(pG11[i]);
		vpG21.emplace_back(pG21[i]);
		
		vmTotal.emplace_back(mTotal[i]);
		

	}

	float mPmin = *min_element(vmP.begin(), vmP.end());
	float mPmax = *max_element(vmP.begin(), vmP.end());

	TH1D *hmP = new TH1D("hmP","Parent Mass", NBins, mPmin - 5.0, mPmax + 5.0);
	TH1D *hEP = new TH1D("hEP","Parent Energy", NBins, (*min_element(vEP.begin(), vEP.end())) - 5.0, (*max_element(vEP.begin(), vEP.end())) + 5.0);
	TH1D *hmC1 = new TH1D("hmC1","C1 Mass", NBins, (*min_element(vmC1.begin(), vmC1.end())) - 5.0, (*max_element(vmC1.begin(), vmC1.end())) + 5.0);
	//TH1D *hmC2 = new TH1D("hmC1","C2 Mass Distribution", NBins, (*min_element(vmC2.begin(), vmC1.end())) - 5.0; (*max_element(vmC2.begin(), vmC1.end())) + 5.0);
	TH1D *hEC1 = new TH1D("hEC","C1 Energy", NBins, (*min_element(vEC1.begin(), vEC1.end())) - 5.0, (*max_element(vEC1.begin(), vEC1.end())) + 5.0);
	TH1D *hpC1 = new TH1D("hpC1","Parent Mass", NBins, (*min_element(vpC1.begin(), vpC1.end())) - 5.0, (*max_element(vpC1.begin(), vpC1.end())) + 5.0);
	//TH1D *hpC2 = new TH1D("hpC2","Parent Mass Distribution", NBins, mPmin - 5.0; mPmax + 5.0);
	//

	TH1D *hmCNew = new TH1D("hmCNew","C1 Mass (RF:C1)", NBins, (*min_element(vmCNew.begin(), vmCNew.end())) - 5.0, (*max_element(vmCNew.begin(), vmCNew.end())) + 5.0);
	TH1D *hECNew = new TH1D("hEP","C1 Energy (RF:C1)", NBins, (*min_element(vECNew.begin(), vECNew.end())) - 5.0, (*max_element(vECNew.begin(), vECNew.end())) + 5.0);
	TH1D *hmG11 = new TH1D("hmG11","G11 Mass", NBins, (*min_element(vmG11.begin(), vmG11.end())) - 5.0, (*max_element(vmG11.begin(), vmG11.end())) + 5.0);
	//TH1D *hmC2 = new TH1D("hmC1","C2 Mass", NBins, (*min_element(vmC2.begin(), vmC1.end())) - 5.0; (*max_element(vmC2.begin(), vmC1.end())) + 5.0);
	TH1D *hEG11 = new TH1D("hEG11","G11 Energy Distribution", NBins, (*min_element(vEG11.begin(), vEG11.end())) - 5.0, (*max_element(vEG11.begin(), vEG11.end())) + 5.0);
	TH1D *hpG11 = new TH1D("hpG11","G11 Momentum", NBins, (*min_element(vpG11.begin(), vpG11.end())) - 5.0, (*max_element(vpG11.begin(), vpG11.end())) + 5.0);
	//TH1D *hpC2 = new TH1D("hpC2","Parent Mass Distribution", NBins, mPmin - 5.0; mPmax + 5.0);

	TH1D *hTheta = new TH1D("hTheta", "Angular Deflection", NBins, -5.0, 185.0);
		
	//determine min,max values of me+e-
	TH1D *hmTotal = new TH1D("hmTotal", "Invariant Mass Sum", NBins, *min_element(vmTotal.begin(), vmTotal.end()), *max_element(vmTotal.begin(), vmTotal.end()));

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

	hmP->GetXaxis()->SetTitle("M");
	hmP->GetYaxis()->SetTitle("Counts");
	hEP->GetXaxis()->SetTitle("E");
	hEP->GetYaxis()->SetTitle("Counts");
	hmC1->GetXaxis()->SetTitle("m_C");
	hmC1->GetYaxis()->SetTitle("Counts");
	hEC1->GetXaxis()->SetTitle("E_C");
	hEC1->GetYaxis()->SetTitle("Counts");
	hpC1->GetXaxis()->SetTitle("p_C");
	hpC1->GetYaxis()->SetTitle("Counts");

	hmCNew->GetXaxis()->SetTitle("m_C (REST)");
	hmCNew->GetYaxis()->SetTitle("Counts");
	hECNew->GetXaxis()->SetTitle("E_C (REST)");
	hECNew->GetYaxis()->SetTitle("Counts");
	hmG11->GetXaxis()->SetTitle("m_{G11}");
	hmG11->GetYaxis()->SetTitle("Counts");
	hEG11->GetXaxis()->SetTitle("E_{G11}");
	hEG11->GetYaxis()->SetTitle("Counts");
	hpG11->GetXaxis()->SetTitle("p_{G11,y}");
	hpG11->GetYaxis()->SetTitle("Counts");

	hTheta->GetXaxis()->SetTitle("#Theta");
	hTheta->GetYaxis()->SetTitle("Counts");
		
	hmTotal->GetXaxis()->SetTitle("m_{e^{+}e^{-}}");
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
