void BreitWignerTest() {

	int NEvents = 1000;
	int NBins = 100;

	float mean = 17.0;
	float width = 2.0;

	float pi = TMath::Pi();
	TF1 *BW = new TF1("BW","((1.0/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, 20); 
	BW->SetParameters(pi, mean, width);

	BW->Draw();

}

