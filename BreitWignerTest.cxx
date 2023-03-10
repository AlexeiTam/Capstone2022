#include "TMath.h"

void BreitWignerTest() {

	float pi = TMath::Pi();
	//float infinity = std::numeric_limits<float>::infinity();
	float m = 17.0;	//related to resonance
	float b = 2.0;	//related to width

	TF1 *IntBW = new TF1("IntBW", "0.5 + ((1.0/pi)*(atan((x-[0])/([1]))))" , -100000, 100000);

	IntBW->SetParameters(m,b);



	int NEvents = 1000;	//1,000
	int NBins = 100;

	std::vector<float> U;	//uniform dist

	std::vector<float> f;	//final dist

	std::vector<float> vBW;

	float fBW
	float fU, ff;
	for(int i = 0; i < NEvents; i++) {

		fU = gRandom->Rndm();
	U.emplace_back(fU);

		ff = IntBW->GetX(U.at(i));
	f.emplace_back(ff);

		fBW = gRandom->Cauchy(m);
	}

	float Umin = *min_element(U.begin(), U.end());
	float Umax = *max_element(U.begin(), U.end());
	
	float fmin = *min_element(f.begin(), f.end());
	float fmax = *max_element(f.begin(), f.end());

	TH1 *h1 = new TH1D("h1", "Breit-Wigner Distribution", NBins, fmin - 1.0, fmax + 1.0);

	TH1 *h2 = new TH1D("h1", "Uniform Distribution", NBins, Umin - 1.0, Umax + 1.0);

	TCanvas *c1 = new TCanvas("c1","",900, 900);
	c1->Divide(3,1);

	float placeholder;
	for(int i = 0; i < NEvents; i++) {

		placeholder = U.at(i);
		h1->Fill(f.at(i));
		h2->Fill(placeholder);
	}

	c1->cd(1);
	h2->Draw();

	c1->cd(2);
	IntBW->Draw();

	c1->cd(3);
	h1->Draw();

	

}
