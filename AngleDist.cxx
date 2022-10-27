void AngleDist(){
  
  //initializing values and arrays
  int NEvents = 10000;
  int NBins = 1000;
  
  float theta;
  float phi;

  float pi = TMath::Pi();
  
  //histogram,
  //x-axis is phi, y-axis is theta
    TH2D *h1 = new TH2D("h1", "Uniform Distribution", NBins, -0.1, 2.0*pi + 0.1, -1.0*pi - 0.1, pi + 0.1);
    h1->GetXaxis()->SetTitle("#phi");
    h1->GetYaxis()->SetTitle("#theta");
  
  for(int i = 0; i < NEvents; i++){
   theta =  2.0*pi*(gRandom->Rndm()) - pi;
   phi = 2.0*pi*(gRandom->Rndm());
    h1->Fill(theta,phi);
  }
  
  h1->Draw();
}
