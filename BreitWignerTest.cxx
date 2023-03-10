void BreitWignerTest() {

	int NEvents = 10;
	int NBins = 100;

	float mean = 17.0;
	float width = 2.0;

	float norm = (2*width*width)/((fabs(width))*(fabs(width))*(fabs(width)));  //normalization const.

	TCanvas *c1 = new TCanvas("c1","",900,900);

	std::vector<float> U;
	std::vector<float> vBW;


	float pi = TMath::Pi();
	float xmax = mean + (sqrt(((norm*width)/(pi*0.10))-(width*width)));
	TF1 *BW = new TF1("BW","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, xmax);

	BW->SetParameters(pi, mean, width, norm);

	BW->Draw();
	float x;

	float BWMax = BW->GetMaximum();

	cout << "BWMax:" << BWMax << endl;

	float chance, prob, coin;
	float dx = 0.01;	//slice of space for probability density

	int index = 0;

	//test
	cout << "i...x...coin...chance" << endl;
	//while(index<NEvents+1) {
	for(int i = 0; i < NEvents; i++) {	
		
		//index++;

		x = 25.0*(gRandom->Rndm());
		
		prob = BW->Eval(x);

		chance = (prob)*(prob)*dx;

		coin = (BWMax*BWMax)*gRandom->Rndm();	//??

		cout << i << "..." << x << "..." << coin << "..." << chance << endl;
		
		//if(coin < chance) { vBW.emplace_back(x); }
		//else if(coin > chance){
			
		//index = index - 1;
		//continue;
		//}

	}




}

