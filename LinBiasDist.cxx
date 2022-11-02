void LinBiasDist(int NEvents = 1000){

	//GOAL:bias a uniform distribution to make a linear-looking function

	//initializing variables, constants
	int NBins = 100;

	float xmin = 0.0;	//give range of x
	float xmax = 10.0;

	float m = (NEvents)/(xmax - xmin);	//linear slope
	//test
	std::cout << "m=" << m << std::endl;
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
		//test
		//std::cout << xb[i] << std::endl;
	}

	//filling X, fX, RejectChance
	
	for(int i = 0; i < NBins; i++) {

		X[i] = xb[i] + (0.5)*(dx);

		fX[i] = f1->Eval(X[i]);

	}

	for(int i = 0; i < NBins; i++) {

		RejectChance[i] = (fX[NBins - 1] - fX[i])/(fX[NBins - 1]);

		//test
		//std::cout << "X[" << i << "]:" << X[i] <<"..." << "RejectChance[" << i << "]:" << RejectChance[i] << std::endl;
	}





	float x;	//variable to fill histogram
	//float fx;	//function value at x
	float chance;	//random number assigned to each (x,fx)


	for(int i = 0; i < NEvents; i++) {

		x = xmax*(gRandom->Rndm());
		//fx = f1->Eval(x);

		for(int j = 0; j < NBins; j++) {	//bin scan
			
			if( x < xb[j+1]  && xb[j] < x ){
				
				chance = 1.0*(gRandom->Rndm());

				if(chance > RejectChance[j]){
				
					h1->Fill(x);
				
				} else continue;
				
			
			} else continue;	
			

		}
	}

	//visualization
	h1->SetLineColor(1);
	h1->SetFillColor(5);
	
	h1->GetXaxis()->SetTitle("x");
	h1->GetYaxis()->SetTitle("Counts");
	h1->Draw();
	

}
