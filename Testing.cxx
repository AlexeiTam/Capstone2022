void Testing(float xmin = 0.5, float xmax = 0.8) {

	float x[10];

	for(int i = 0; i < 10; i++) {

		x[i] = gRandom->Rndm();
		std::cout << "x[" << i << "]:" << x[i] << std::endl;
		if( x >> xmin && x << xmax){
			std::cout << "Mark!" << std::endl;	
		}
	}
}
