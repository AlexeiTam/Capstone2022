void LorentzBoost(float beta = 0.5, float px = 2.0, float py = 1.0, float pz = 1.0, float E = 1.0) {

	//beta - v/c
	std::cout << "INPUT: Beta, Px, Py, Pz, E" << std::endl;
	
	std::vector<float> v1;
	std::vector<float> v2;

	float L[4][4];	//4x4 lorentz Boost in X direction
			//L[row][column]
	
	float gamma = sqrt((1.0)/(1.0 - (beta)*(beta)));

	std::cout << "gamma:" << gamma << std::endl;

	std::cout << "beta:" << beta << std::endl;
	for (int i = 0; i < 4; i++) {

		for(int j = 0; j < 4; j++) {

			L[i][j] = 0.0;
		}
	}

	//fix:
	L[0][0] = gamma;
	L[1][1] = gamma;
	L[0][1] = -1.0*gamma*beta;
	L[1][0] = -1.0*gamma*beta;

	L[2][2] = 1.0;
	L[3][3] = 1.0;
	
	//print matrix:
	for(int i = 0; i < 4; i++){

		std::cout << L[i][0] << "|" << L[i][1] << "|" << L[i][2] << "|" << L[i][3] << std::endl;
	}
	
	//Filling P-vector:
	
	float P[4];

	P[0] = px;
	P[1] = py;
	P[2] = pz;
	P[3] = E;

	float PNew[4] = {0,0,0,0};

	for(int i = 0; i < 4; i++) {
	
		for(int j = 0; j < 4; j++) {

			PNew[i] = PNew[i] + (L[i][j])*P[j];
		}
		
	}

	std::cout << "Px:"  << P[0] << "... Py:"  << P[1] << "... Pz:" << P[2] << "...E:" << P[3] << std::endl;

	
	std::cout << "Px':"  << PNew[0] << "... Py':"  << PNew[1] << "... Pz':" << PNew[2] << "...E':" << PNew[3] << std::endl;


}
