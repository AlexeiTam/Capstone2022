void LinBiasDist(){

	//GOAL:bias a uniform distribution to make a linear-looking function

	//initializing variables, constants
	int NEvents = 1000;
	int NBins = 100;

	float xmin = 0.0;	//give range of x
	float xmax = 10.0;

	float m = 1.0;	//linear slope
	float b = 0.0;	//y-intercept
	
	//linear function
	TF1* f1 = new TF1("fit","[0]*x + [1]", xmin, xmax);
	f1->SetParameters(m,b);

	//histogram
	TH1D* h1 = new TH1D("h1","", NBins, xmin, xmax);

	float xRight[NBins];

	//finding right edges of bins
	for(int i = 0; i < NBins; i++){
	
		xRight[i] = (i+1)*((xmax - xmin)/(NBins));
		std::cout << xRight[i] << std::endl;
	}
}
