void HistBoundaryTest(){

	float v1min = -0.0;
	float v1max = 0.0;
	std::vector<float> v1;
	
	TH1D *hmP = new TH1D("hmP","^{8}Be^{*} Mass", NBins, v1min, v1max);
	
	for(int i = 0; i < NEvents; i++){
		
		v1.emplace_back(gRandom->Rndm());
		if(v1.at(i) > v1max){
			v1max = v1.at(i);
		}

		cout << v1.at(i) << "..." << v1max << endl;
	}

	int index = 0;

	while(index < v1.Size()){

		hmP->Fill(v1.at(index));
		cout << v1.at(index) << endl;
	}

	hmP->Update();
	hmP->Draw();

}
