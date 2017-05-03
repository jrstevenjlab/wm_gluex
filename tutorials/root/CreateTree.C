void CreateTree() 
{
  // An all-too-simple physics simulation that writes events to a ROOT Tree.
  // 03-Jun-2001 by Bill Seligman for his root course

  // Note: in a more realistic simulation, an event generator such as
  // TPythia might be used to simulate the even generation, and a
  // Monte Carlo such as Geant4 might be used to simulate the
  // detector.

  const Int_t numberEvents = 100000;

  Int_t event = 0;
  Float_t ebeam;
  Float_t energy;
  Float_t theta;
  Float_t phi;
  Float_t px; 
  Float_t py; 
  Float_t pz; 
  Float_t zv;
  Float_t chi2;
  
  TFile f("tree.root","recreate");
  TTree* t1 = new TTree("tree1","Reconstruction ntuple");
  t1->Branch("event",&event,"event/I");
  t1->Branch("ebeam",&ebeam,"ebeam/F");
  t1->Branch("px",&px,"px/F");
  t1->Branch("py",&py,"py/F");
  t1->Branch("pz",&pz,"pz/F");
  t1->Branch("zv",&zv,"zv/F");
  t1->Branch("chi2",&chi2,"chi2/F");
  
  for (event=0; event<numberEvents; event++) 
    {
      zv   = gRandom->Exp(10.);
      ebeam  = gRandom->Gaus(150.,.15);
      energy = ebeam - zv*10./50. - gRandom->Gaus(1.,.2);
      theta = TMath::Abs(gRandom->Gaus(0.,.05));
      phi   = gRandom->Uniform(2.*TMath::Pi());

      px   = energy * TMath::Sin(theta) * TMath::Cos(phi);
      py   = energy * TMath::Sin(theta) * TMath::Sin(phi);
      pz   = energy * TMath::Cos(theta);

      if (theta > 0.15)
	// At large angles, we can't fit tracks properly.
	chi2 = gRandom->Gaus(10.,3.);
      else
	chi2 = gRandom->Gaus(1.,.1);

      // Add the variables to the tree.
      t1->Fill();
    }   
  t1->Print();
  t1->Write();
}
