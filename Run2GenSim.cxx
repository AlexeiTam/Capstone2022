void Run2GenSim(){
	//SAME AS GenSim.cxx, but with imporvements (Run2Improvements.txt)
	//GOAL: simulate general two-step decay: P->C1 + C2; C1->g11 + g12

	//ED: Natural Units; i.e, c = 1; can change to c = 3e-8 later
	const Int_t NEvents = 1000; //number of simulations: 10,000
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

	float mPMean = 17.0;	//can be made initializable
	float mPSigma = 2.0;	//same as above
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

		ECNew[i] = PCNew[i][0];
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
	float EG[NEvents];

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

		//canvas: split based on RF
		TCanvas *cP = new TCanvas("cP","REST FRAME:PARENT", 1500, 1500);
		TCanvas *cC1 = new TCanvas("cC1", "REST FRAME:C1", 1500, 1500);

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
	
	hmP->SetFillColor(0);
	hmP->SetLineWidth(5);
	hEP->SetFillColor(1);
	hmP->SetLineWidth(5);
	hmC1->SetFillColor(2);
	hmC1->SetLineWidth(5);
	hEC1->SetFillColor(3);
	hEC1->SetLineWidth(5);
	hpC1->SetFillColor(4);
	hpC1->SetLineWidth(5);

	hmCNew->SetFillColor(5);
	hmCNew->SetLineWidth(5);
	hECNew->SetFillColor(6);
	hECNew->SetLineWidth(5);
	hmG11->SetFillColor(7);
	hmG11->SetLineWidth(5);
	hEG11->SetFillColor(8);
	hEG11->SetLineWidth(5);
	hpG11->SetFillColor(9);
	hpG11->SetLineWidth(5);

	hTheta->SetFillColor(10);
	hTheta->SetLineWidth(5);

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

	cP->cd(6);
	hTheta->Draw();
	
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


}