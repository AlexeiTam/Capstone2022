void LorentzianDist(int NEvents = 100, float a = 0.0, float b = 1.0, float scale = 1.0, float shift = 0.0) {

	//GOAL: make a gaussian distribution
	
	//initializing values
	int NBins = 100;
	float xmin = scale*(-6.1);
	float xmax = scale*(6.1);
	float x;

	//histogram
	TH1D *h1 = new TH1D("h1","Lorentzian Distribution", NBins, xmin, xmax);
	
	//filling histogram, generating distribution
	for(int i = 0; i < NEvents; i++) {
	
	x = scale*(gRandom->CauchyDist(mean,sigma)) + shift;
	h1->Fill(x);		
	}

	//drawing hist
	h1->SetFillColor(4);
	h1->Draw();

	}
