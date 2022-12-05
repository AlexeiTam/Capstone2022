void LorentzBoostBASIC(float beta = 0.5, float px = 2.0, float py = 1.0, float pz = 1.0, float E = 1.0) {

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

	P[0] = E;
	P[1] = px;
	P[2] = py;
	P[3] = pz;

	float PNew1[4] = {0,0,0,0};
	float PNew2[4] = {0,0,0,0};


	for(int i = 0; i < 4; i++) {
	
		for(int j = 0; j < 4; j++) {

			PNew1[i] = PNew1[i] + (L[i][j])*P[j];
		}
		
	}


	for(int i = 0; i < 4; i++) {

		PNew2[0] = L[0][0]*P[0] + L[0][1]*P[1] + L[0][2]*P[2] + L[0][3]*P[3];
		PNew2[1] = L[1][0]*P[0] + L[1][1]*P[1] + L[1][2]*P[2] + L[1][3]*P[3];
		PNew2[2] = L[2][0]*P[0] + L[2][1]*P[1] + L[2][2]*P[2] + L[2][3]*P[3];
		PNew2[3] = L[3][0]*P[0] + L[3][1]*P[1] + L[3][2]*P[2] + L[3][3]*P[3];
	}

	std::cout << "E:"  << P[0] << "... Px:"  << P[1] << "... Py:" << P[2] << "...Pz:" << P[3] << std::endl;

	
	std::cout << "E':"  << PNew1[0] << "... Px':"  << PNew1[1] << "... Py':" << PNew1[2] << "...Pz':" << PNew1[3] << std::endl;


	std::cout << "E'':"  << PNew2[0] << "... Px'':"  << PNew2[1] << "... Py'':" << PNew2[2] << "...Pz'':" << PNew2[3] << std::endl;
}
