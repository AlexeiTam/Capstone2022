void VectorTesting() {

std::vector<float> v1;

float x;

for(int i = 0; i < 10; i++) {

    x = gRandom->Rndm();
    v1.emplace_back(x);
}

for(int i = 0; i < 10; i++) {

    std::cout << v1.at(i) << std::endl;
}



}