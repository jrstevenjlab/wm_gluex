void CreateSubdirectories() 
{
  // This macro creates sub-directories (or sub-folders) within an
  // output file.  Each sub-directory is populated with a bunch of
  // histograms.  The histograms themselves are rather dull (just
  // gaussian distributions and whatnot); they are used in a student
  // exercise to read, fit, and plot parameters from large numbers of
  // histograms.

  const Int_t numberEvents = 100000;

  // Create a fresh output file; overwrite any previous version.
  TFile* outputFile = new TFile("folders.root","recreate");

  // Create a sub-directory within that file.
  TDirectory* FallaDirectory 
    = outputFile->mkdir("example1","histograms used to make the plot with the silly esoteric pun");

  // Make that sub-directory the default location for any new objects (like histograms).
  FallaDirectory->cd();
  
  Int_t numberFallaHistograms = 6;

  for (Int_t i=0; i != numberFallaHistograms; ++i) 
    {
      // We're going to construct a histogram based on the following
      // properties:
      Double_t time = 0.5*pow(10,i-1);
      Double_t mean = 3.*log(time + 100.);
      Double_t sigma = 0.1*mean;

      Double_t lowerHistogramEdge = mean - 3.0*sigma;
      Double_t upperHistogramEdge = mean + 3.0*sigma;

      // Convert the time to a string.
      std::ostringstream os;
      os << time;
      TString timeString = os.str();

      // For the histogram name, use the index number.
      std::ostringstream os1;
      os1 << i;
      TString timeScaled = os1.str();

      // Include the time in the histogram title and the description.
      TString histogramTitle = "hist" + timeScaled;

      TString histogramDescription = "Value at time = " + timeString + " seconds";

      TH1* histogram = new TH1D(histogramTitle,
				histogramDescription,
				100, // number of bins
				lowerHistogramEdge,
				upperHistogramEdge);
      histogram->SetXTitle("n_{ions}");
      histogram->SetYTitle("number of events");

      // Fill the histogram with a gaussian.
      TF1* gaussian = new TF1("gaussian","gaus",lowerHistogramEdge,upperHistogramEdge);
      gaussian->SetParameters(1.0,mean,sigma);
      histogram->FillRandom("gaussian",numberEvents);

      histogram->Write();

      // When you create something with "new", you want to get rid of
      // it with "delete".
      delete gaussian;
      delete histogram;
    }   

  // Do this above again, but this time for a more "serious" example.

  TDirectory* DistributionDirectory 
    = outputFile->mkdir("example2","Intended to teach how to use an n-tuple to access histograms by name");
  DistributionDirectory->cd();
  
  Int_t numberDistributionHistograms = 40;

  // Something a bit different: instead of asking them to derive the
  // x-value from the histogram ID, supply the histogram->x
  // correspondence in an n-tuple.

  TTree* t1 = new TTree("histogramList","A 'map' of histograms and associated values of x.");
  Double_t x;
  Int_t histNumber = 0;
  t1->Branch("histNumber",&histNumber, "histNumber/I");
  t1->Branch("x",&x,"x/D");

  // For random numbers.
  TRandom3 random;

  for (Int_t i=0; i != numberDistributionHistograms; ++i) 
    {
      // Generate a pair of random numbers according to a normal
      // distribution.
      Double_t random1,random2;
      random.Rannor(random1,random2);

      // We're going to construct a histogram based on the following
      // properties:
      x = Double_t(i) / 2. + (random2 + 0.5)/numberDistributionHistograms;
      Double_t ey = 0.005;
      Double_t y = TMath::Landau(x, 0., 1.);

      Double_t lowerHistogramEdge = y - 4.0*ey;
      Double_t upperHistogramEdge = y + 4.0*ey;

      // Put a little randomness into the mean of the gaussian.
      y += random1 * ey;

      // Add a random amount to the histogram number, so the students
      // can't predict the value of x from it.
      histNumber += 3 + Int_t( 5. * random.Rndm() );

      // For the histogram name, use the histogram number.
      std::ostringstream os;
      os << "hist" << histNumber;
      TString name = os.str();

      TH1* histogram = new TH1D(name,
				"Number of events for this x",
				100, // number of bins
				lowerHistogramEdge,
				upperHistogramEdge);
      histogram->SetXTitle("x [arbitrary units]");
      histogram->SetYTitle("number of events");

      // Fill the histogram with a gaussian.
      TF1* gaussian = new TF1("gaussian","gaus",lowerHistogramEdge,upperHistogramEdge);
      gaussian->SetParameters(1.0,y,ey);
      histogram->FillRandom("gaussian",numberEvents);

      histogram->Write();

      // Write the entry into the n-tuple.
      t1->Fill();

      // When you create something with "new", you want to get rid of
      // it with "delete".
      delete gaussian;
      delete histogram;
    }
  // Write the n-tuple.
  t1->Write();
  delete t1;

  // One more time, this time with a large number of histograms and no
  // "guiding n-tuple."  Instead, we'll put the value of x in the
  // histogram description.

  // To be even more realistic, we'll create two sets of histograms:
  // before cuts, and after.  The student will only be asked to
  // include the values after cuts in the plot.

  TDirectory* PolyDirectory 
    = outputFile->mkdir("example3","Intended to teach the students how to iterate through a directory");
  PolyDirectory->cd();
  
  Int_t numberPolyHistograms = 250;
  histNumber = 1000;

  std::vector<Int_t> histNumbers;

  for (Int_t i=0; i != numberPolyHistograms; ++i) 
    {
      // Generate a pair of random numbers according to a normal
      // poly.
      Double_t random1,random2;
      random.Rannor(random1,random2);

      // We're going to construct a histogram based on the following
      // properties:
      x = (random.Rndm() - 0.5) * 4.0;
      Double_t ey = 0.35;
      Double_t y = 1.0*x + 2.0*x*x - x*x*x;

      Double_t lowerHistogramEdge = y - 4.0*ey;
      Double_t upperHistogramEdge = y + 4.0*ey;

      // Put a little randomness into the mean of the gaussian.
      y += random1 * ey;

      // Add a random amount to the histogram number, so the students
      // can't predict the value of x from it.  Save the number to
      // re-use it in the next loop.
      histNumber += 17 + Int_t( 10. * random.Rndm() );
      histNumbers.push_back( histNumber );

      // For the histogram name, use the histogram number.
      std::ostringstream os;
      os << "plotBeforeCuts" << histNumber;
      TString name = os.str();

      // Put the value of x in the histogram's descriptive text.
      TString title;
      title += x; // Look at the description of operator+=() in TString.

      TH1* histogram = new TH1D(name,
				title,
				100, // number of bins
				lowerHistogramEdge,
				upperHistogramEdge);
      histogram->SetXTitle("x #left[#frac{fath}{Ha}#right]");
      histogram->SetYTitle("number of events");

      // Fill the histogram with a gaussian.
      TF1* gaussian = new TF1("gaussian","gaus",lowerHistogramEdge,upperHistogramEdge);
      gaussian->SetParameters(1.0,y,ey);
      histogram->FillRandom("gaussian",numberEvents*2);

      histogram->Write();

      // When you create something with "new", you want to get rid of
      // it with "delete".
      delete gaussian;
      delete histogram;
    }

  for (Int_t i=0; i != numberPolyHistograms; ++i) 
    {
      // Generate a pair of random numbers according to a normal
      // poly.
      Double_t random1,random2;
      random.Rannor(random1,random2);

      // We're going to construct a histogram based on the following
      // properties:
      x = (random.Rndm() - 0.5) * 4.0;
      Double_t ey = 0.5;
      Double_t y = 1.0*x + 2.0*x*x - x*x*x;

      Double_t lowerHistogramEdge = y - 4.0*ey;
      Double_t upperHistogramEdge = y + 4.0*ey;

      // Put a little randomness into the mean of the gaussian.
      y += random1 * ey;

      // Use the same sequence of histogram numbers as in the previous
      // loop.
      histNumber = histNumbers[i];

      // For the histogram name, use the histogram number.
      std::ostringstream os;
      os << "plotAfterCuts" << histNumber;
      TString name = os.str();

      // Put the value of x in the histogram's descriptive text.
      TString title;
      title += x; // Look at the description of operator+=() in TString.

      TH1* histogram = new TH1D(name,
				title,
				100, // number of bins
				lowerHistogramEdge,
				upperHistogramEdge);
      histogram->SetXTitle("x #left[#frac{fath}{Ha}#right]");
      histogram->SetYTitle("number of events");

      // Fill the histogram with a gaussian.
      TF1* gaussian = new TF1("gaussian","gaus",lowerHistogramEdge,upperHistogramEdge);
      gaussian->SetParameters(1.0,y,ey);
      histogram->FillRandom("gaussian",numberEvents);

      histogram->Write();

      // When you create something with "new", you want to get rid of
      // it with "delete".
      delete gaussian;
      delete histogram;
    }

  outputFile->Close();
}
