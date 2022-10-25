void GausDist(int NEvents = 100, float mean = 0.0, float sigma = 1.0, float scale = 1.0, float shift = 0.0) {

	//initializing values
	int NBins = 100;
	float xmin = scale*(-6.1);
	float xmax = scale*(6.1);
	float x;

	//histogram
	TH1D *h1 = new TH1D("h1","Gaussian Distribution", NBins, xmin, xmax);
	
	//filling histogram, generating distribution
	for(int i = 0; i < NEvents; i++) {
	
	x = scale*(gRandom->Gaus(mean,sigma)) + shift;
	h1->Fill(x);		
	}

	//drawing hist
	h1->Draw();

	}
