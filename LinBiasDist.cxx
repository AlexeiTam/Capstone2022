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

	float RejectChance[NBins];	//if chance < RejectChance, reject
					//else, accept
	
	//filling xb
	for(int i = 0; i < NBins+1; i++){
	
		xb[i] = i*dx;
		std::cout << xb[i] << std::endl;
	}

	//filling X, fX, RejectChance
	
	for(int i = 0; i < NBins; i++) {

		X[i] = xb[i] + (0.5)*(dx);

		fX[i] = f1->Eval(X[i]);

	}

	for(int i = 0; i < NBins; i++) {

		RejectChance[i] = (fX[NBins - 1] - fX[i])/(fX[NBins - 1]);

	}



	float x;	//variable to fill histogram
	//float fx;	//function value at x
	float chance;	//random number assigned to each (x,fx)


	for(int i = 0; i < NEvents; i++) {

		x = gRandom->Rndm();
		//fx = f1->Eval(x);
		chance = gRandom->Rndm();

		for(int j = 0; j < NBins; j++) {	//bin scan
			
			if( x < xb[j+1]  && xb[j] < x ){

				if(chance < RejectChance[j]){
				h1->Fill(x);
				}
				else continue;	
			
			}	
			else continue;

		}
	}

	//visualization
	h1->SetLineColors(0);
	h1->SetFillColor(1);
	h1->Draw();


}
