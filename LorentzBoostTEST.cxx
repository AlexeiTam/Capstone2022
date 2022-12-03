void LorentzBoostTEST() {

	
	ROOT::Math::PxPyPzEVector v1(1.1, 2.2, 3.3, 4.4);


	float px, py, pz, E;

	px = v1.Px();
	py = v1.Py();
	pz = v1.Pz();
	E = v1.E();

	//print 4-vector
	
	std::cout << "px:" << px << "..." << "py:" << py << "..." << "pz:" << pz << "..." << "E:" << E << std::endl;
	

	float beta = 0.5; 	//beta = v/c
	ROOT::Math::BoostX Mat(beta);
	
	ROOT::Math::PxPyPzEVector v2 = Mat(v1);

	
	float px2, py2, pz2, E2;


	std::cout << "px':" << px2 << "..." << "py':" << py2 << "..." << "pz':" << pz2 << "..." << "E':" << E2 << std::endl;

}
