void SolidAngleDist(int NEvents = 100000){
	
	//GOAL: make a distribution uniform in solid angle

	//canvas
	
    TCanvas *c1=new TCanvas("c1","Angular Distribution",900,900);
    c1->Divide(3,1);


  //initializing values and arrays
  int NBins = 100;
  
  float theta;
  float phi;
  float x;

  std::vector<float> vTheta;
  std::vector<float> vPhi;

  float pi = TMath::Pi();
  
  //histogram,
  //x-axis is phi, y-axis is theta
    TH2D *h1 = new TH2D("h1", "Uniform Distribution", NBins, -0.1, 2.0*pi + 0.1, NBins, - 0.1, pi + 0.1);
    h1->GetXaxis()->SetTitle("#phi");
    h1->GetYaxis()->SetTitle("#theta");
  
  for(int i = 0; i < NEvents; i++){
   phi =  2.0*pi*(gRandom->Rndm());
   x = 2.0*(gRandom->Rndm()) - 1.0;
   theta = TMath::ACos(x); 
  

   vTheta.emplace_back(theta);
   vPhi.emplace_back(phi);


   h1->Fill(phi,theta);
  }
  
  //Draw
  c1->cd(1);
  h1->Draw();

  c1->cd(2);
  h1->Draw("LEGO1");

  TGraph2D *gr = new TGraph2D();
  float gTheta;
  float gPhi;

  float cosPhi;
  float sinPhi;
  float cosTheta;
  float sinTheta;

  float Size = (vTheta.size()) - 1;
  for(int i = 0; i < Size; i++) {

	  gTheta = vTheta.at(i);
	  gPhi = vPhi.at(i);
	
	  cosPhi = TMath::Cos(gPhi);
	  sinPhi = TMath::Sin(gPhi);

	  cosTheta = TMath::Cos(gTheta);
	  sinTheta = TMath::Sin(gTheta);

	  gr->AddPoint(cosPhi*sinTheta, sinPhi*sinTheta, cosTheta);
  }

  c1->cd(3);

  gr->GetXaxis()->SetTitle("x");
  gr->GetYaxis()->SetTitle("y");
  gr->GetZaxis()->SetTitle("z");
  gr->SetMarkerStyle(3);
  gr->Draw();

}
