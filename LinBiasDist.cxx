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

	float xb[NBins+1];	//boundaries of bins
	float dx = (xmax - xmin)/(NBins);	//width of bins

	float X[NBins];	//midpoints of bins
	float fX[NBins]; //value of func. at midpt. of bin

	//filling xb
	for(int i = 0; i < NBins+1; i++){
	
		xb[i] = i*dx;
		std::cout << xb[i] << std::endl;
	}

	//filling X, fX
	
	for(int i = 0; i < NBins; i++) {

		X[i] = xb[i] + (0.5)*(dx);

		fX[i] = f1->Eval(X[i]);
	}

	float x;	//variable to fill histogram
	float RejectChance;	//


	for(int i = 0; i < NEvents; i++) {

		x = gRandom->Rndm();
		
		for(int j = 0; j < NBins; j++) {
			
			if( <  && < ){
			
			}	
			else continue;

		}
	}

}
