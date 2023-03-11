void BreitWignerTest() {

	int NEvents = 10000;
	int NBins = 100;

	float mean = 17.0;
	float width = 2.0;

	float norm = (2*width*width)/((fabs(width))*(fabs(width))*(fabs(width)));  //normalization const.

	TCanvas *c1 = new TCanvas("c1","",900,900);
	c1->Divide(2,1);

	std::vector<float> U;
	std::vector<float> vBW;


	float pi = TMath::Pi();
	float cutoff = 0.10;
	float xmax = mean + (sqrt(((norm*width)/(pi*cutoff))-(width*width)));	//end when P(x) = cutoff
	
	TF1 *BW = new TF1("BW","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, xmax);

	BW->SetParameters(pi, mean, width, norm);

	c1->cd(1);
	BW->Draw();

	float x;

	float BWMax = BW->GetMaximum();

	cout << "BWMax:" << BWMax << endl;

	float chance, prob;
//	float dx = xmax/NBins;	//slice of space for probability density

	int index = 0;

	//test
//	cout << "i...x...prob...chance" << endl;
	
//	int it;
	while(index < NEvents+1) {
		
		index++;

//		it = index;
		x = 25.0*(gRandom->Rndm());
		
		prob = gRandom->Rndm();

		chance = BW->Eval(x);

		//coin = (BWMax)*gRandom->Rndm();	//??

//		cout << vBW.size() << "..." << x << "..." << prob << "..." << chance << endl;
		
		if(prob > chance) { 
			index = index - 1;	
		}
		else if(prob < chance) {
			vBW.emplace_back(x);
		}
	}


	float Size = vBW.size() ;
	
	TH1 *hist = new TH1D("hist","Breit-Wigner Dist.", NBins, 0, xmax);

//test	
//	std::cout << "i...vBW[i]" << std::endl;

//	for(int i = 0; i < Size; i++) {

//		std::cout << i << "..." << vBW.at(i) << std::endl;

//	}

	float placeholder;
	for(int i = 0; i < Size; i++) {

		placeholder = vBW.at(i);
		hist->Fill( placeholder );

	}

	c1->cd(2);
	hist->Draw();
	

	




}

