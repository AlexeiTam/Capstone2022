void GenSim(){
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

	float mPMean = 5.0;	//can be made initializable
	float mPSigma = 0.50;	//same as above
	float mCMean = 0.25*mPMean;
	float mCSigma = 0.25*mCSigma;


	//generation: mP, EP
	for(int i = 0; i < NEvents; i++) {

		//ED: gaussian instead of BW; work out later
		mP[i] = gRandom->Gaus(mPMean, mPSigma);
		EP[i] = (mP[i])*c*c;


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

}
