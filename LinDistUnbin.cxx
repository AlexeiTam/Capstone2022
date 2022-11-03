void LinDistUnbinPRELIM(){

	//Make similar to RooFit; gen. random outputs and based on div. from model output, give chance of rejection
	//Error bars based on Gaussian? or Poisson?
	//Preliminary, testing a lot

	
	int NEvents = 1000;
	float xmin = 0.0;
	float xmax = 10.0;


	TF1* f1 = new TF1("f1","[0]*x + [1]", xmin, xmax);


}
