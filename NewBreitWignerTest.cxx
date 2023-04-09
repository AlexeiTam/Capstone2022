void NewBreitWignerTest(){
  
   float mC1Mean = 16.98;	//X17 ==>invariant mass approx. 17 MeV
	  float mC1Sigma = 1.0; 	//will play with later?
float pi = TMath::Pi();
  float cutoff = 0.10;

float mC1Norm = (2*mC1Sigma*mC1Sigma)/((fabs(mC1Sigma))*(fabs(mC1Sigma))*(fabs(mC1Sigma)));
float mC1Max = mC1Mean + (sqrt(((mC1Norm*mC1Sigma)/(pi*cutoff))-(mC1Sigma*mC1Sigma)));

		TF1 *BWC1 = new TF1("BWC1","((([3])/[0])*(([2])/(((x-[1])*(x-[1]))+([2]*[2]))))", 0.0, mC1Max);
		BWC1->SetParameters(pi, mC1Mean, mC1Sigma, mC1Norm);
  
  float x[10000];
  
  for(int i = 0; i < 10000; i++){
   x[i] = BWC1->GetRandom(); 
  }

  TH1D *h1 = new TH1D("h1","Breit-Wigner",1000,(*min_element(x,x+10000)) - 5, (*max_element(x,x+10000)) + 5);
  for(int i = 0; i < 10000; i++){
   h1->Fill(x[i]);
  }
  
  h1->Draw();
  
}
