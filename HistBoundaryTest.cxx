void HistBoundaryTest(){

	std::vector<float> v1;
	for(int i = 0; i < NEvents; i++){
		
		v1.emplace_back(i*2);

	}	

	int index = 0;

	while(index < v1.Size()){

		cout << v1.at(index) << endl;
	}


}
