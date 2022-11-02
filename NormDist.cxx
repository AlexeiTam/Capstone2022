void NormDist(int NEvents = 1000, float scale = 2.0, float shift = 1.0) {
	
	//GOAL: make a normal distribution
	
	//initializing values
	int NBins = 100;
	//int NEvents = 1000;
	float xmin = scale*(-1.1);
	float xmax = scale*(1.1);
	float x;

	//histogram
	TH1D *h1 = new TH1D("h1","Normal Distribution", NBins, xmin, xmax);
	
	//filling histogram, generating distribution
	for(int i = 0; i < NEvents; i++) {
	
	x = scale*(gRandom->Rndm()) - shift;
	h1->Fill(x);		
	}

	//drawing hist
	h1->SetFillColor(2);
	h1->Draw();

	}
