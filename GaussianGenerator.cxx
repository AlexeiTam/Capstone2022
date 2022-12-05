void GaussianGenerator() {

	//GOAL: Generate gaussian distribution using inverse transformation sampling
	
	 TCanvas *c1=new TCanvas("c1","",900,900);
	 c1->Divide(3,1);

	int NBins = 1000;	
	int NEvents = 500000;	//NEvents = 500,000
	
	TF1* f1 = new TF1("f1","erf(x)",-100,100);
	
	

	std::vector<float> vU;
	std::vector<float> vG;

	float x;

	for(int i = 0; i < NEvents; i++) {

	x = 2.0*(gRandom->Rndm()) - 1.0;

	vU.emplace_back(x);

	}


	TH1 *h1 = new TH1D("h1","Uniform Distribution", NBins, -1.1, 1.1);

	for(int i = 0; i < NEvents; i++) {
		
		h1->Fill(vU.at(i));

	}

	c1->cd(1);
	h1->SetFillColor(3);
	h1->GetYaxis()->SetTitle("Counts");
	h1->GetXaxis()->SetTitle("y");

	h1->Draw();

	float y;

	for(int i = 0; i < NEvents; i++) {

		y = f1->GetX(vU.at(i));
		vG.emplace_back(y);

	}

	TH1 *hG = new TH1D("hG","Gaussian Distribution", NBins, *min_element(vG.begin(), vG.end()), *max_element(vG.begin(), vG.end()));
	
	for(int i = 0; i < NEvents; i++) {
		
		hG->Fill(vG.at(i));
	}

	c1->cd(3);
	hG->SetFillColor(7);
	hG->GetXaxis()->SetTitle("x");
	hG->GetYaxis()->SetTitle("Counts");
	hG->Draw();

	
	TF1* func = new TF1("func","erf(x)",*min_element(vG.begin(), vG.end()),*max_element(vG.begin(), vG.end()));
	c1->cd(2);
	func->GetXaxis()->SetTitle("x");
	func->GetYaxis()->SetTitle("y");
	func->Draw();


	

}
