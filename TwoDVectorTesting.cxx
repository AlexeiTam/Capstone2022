void TwoDVectorTesting(){

	const Int_t NEvents = 10;

	std::vector<float> E1;
	std::vector<float> px1;
	std::vector<float> py1;
	std::vector<float> pz1;

	std::vector<float> E2;
	std::vector<float> px2;
	std::vector<float> py2;
	std::vector<float> pz2;

	std::vector<float> vP1;
	std::vector<vector<float>> P1;
	std::vector<vector<float>> P2;


	std::vector<float> v1;
	std::vector<vector<float>> V;

	for(int i = 0; i < NEvents; i++){

		E1.emplace_back(gRandom->Rndm());
		px1.emplace_back(1.0 + gRandom->Rndm());
		py1.emplace_back(2.0 + gRandom->Rndm());
		pz1.emplace_back(3.0 + gRandom->Rndm());

	}

	for(int i = 0; i < NEvents; i++){

		//P1[i][0] = E1.at(i);
	       	//P1[i][1] = px1.at(i);
		//P1[i][2] = py1.at(i);
		//P1[i][3] = pz1.at(i);

		vP1 = {E1.at(i), px1.at(i), py1.at(i), pz1.at(i)};
		P1.push_back({E1.at(i), px1.at(i), py1.at(i), pz1.at(i)});
	}

	std::cout << "i...E1[i]...px1[i]...py1[i]...pz1[i]" << std::endl;

	for(int i = 0; i < NEvents; i++){

		std::cout << i << "..." << E1.at(i) << "..." << px1.at(i) << "..." << py1.at(i) << "..." << pz1.at(i) << std::endl;

	}

	std::cout << "-----------------" << std::endl;

	std::cout << "i...P1[i][1]...P1[i][2]...P1[i][3]..." << std::endl;

	for(int i = 0; i < NEvents; i++){

		std::cout << i << "..." << P1[i][0] << "..." << P1[i][1] << "..." << P1[i][2] << "..." << P1[i][3] << std::endl;

	}

	//DOT PRODUCT?
	
	std::cout << "DOT PRODUCT:" << std::endl;
	
	std::vector<vector<float>> L 
	{
		{-2.0, 0.0, 0.0, 0.0},
		{0.0, -2.0, 0.0, 0.0},
		{0.0, 0.0, 1.0, 0.0 },
		{0.0, 0.0, 0.0, 1.0}
	};
	std::vector<float> p {1.0, 1.0, 1.0, 1.0};

	std::vector<float> pnew;

	float placeholder = 0.0;
	for(int i = 0; i < 4; i++) {
		placeholder = 0.0;
		for(int j = 0; j < 4; j++) {

			placeholder = placeholder + (L[i].at(j))*(p.at(j));
		}

		pnew.emplace_back(placeholder);

	}

	std::cout << "pnew = { -2.0, -2.0, 1.0, 1.0}:" << std::endl;
	std::cout << "actual = {" << pnew.at(0) << "," << pnew.at(1) << "," << pnew.at(2) << "," << pnew.at(3) << std::endl;
	


	
}
