void AngleDist(){
	//canvas
	
	
    TCanvas *c1=new TCanvas("c1","Angular Distribution",900,900);
    c1->Divide(2,1);

  //initializing values and arrays
  int NEvents = 10000;
  int NBins = 1000;
  
  float theta;
  float phi;

  float pi = TMath::Pi();
  
  //histogram,
  //x-axis is phi, y-axis is theta
    TH2D *h1 = new TH2D("h1", "Uniform Distribution", NBins, -0.1, 2.0*pi + 0.1, NBins, -1.0*pi - 0.1, pi + 0.1);
    h1->GetXaxis()->SetTitle("#phi");
    h1->GetYaxis()->SetTitle("#theta");
  
  for(int i = 0; i < NEvents; i++){
   theta =  2.0*pi*(gRandom->Rndm()) - pi;
   phi = 2.0*pi*(gRandom->Rndm());
    h1->Fill(theta,phi);
  }
  
  //Draw
  c1->cd(1);
  h1->Draw();

  c1->cd(2);
  h1->Draw("LEGO1");
}
