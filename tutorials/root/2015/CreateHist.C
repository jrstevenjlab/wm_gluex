{
// Macro to create histogram for 06-Jun-2001 ROOT class
// 01-Jun-2001 Bill Seligman

// Open a file.
TFile* myFile = new TFile("histogram.root","RECREATE","Histograms for ROOT class");

// Create a simple histogram of a gaussian distribution.
TF1* myFunc1 = new TF1("func1","gaus");
myFunc1->SetParameters(5.0, 6.0, 2.5);
TH1F* myHist1 = new TH1F("hist1","Function to be fit",100,0.,15.);
myHist1->FillRandom("func1",10000);

// I didn't figure out the following code by myself.  I created a
// histogram that looked the way I wanted it, used "Save as canvas.C',
// and looked at the code in the saved file.
myHist1->GetXaxis()->SetTitle("Gaussian");
myHist1->SetMarkerColor(9);
myHist1->SetMarkerStyle(21);
myHist1->SetMarkerSize(0.7);                                                           

// Draw the error bars.
myHist1->SetOption("e1");

myHist1->Write();


// Create another histogram of a slightly more complex function.
TF1* myFunc2 = new TF1("func2","gaus(0)+gaus(3)");
myFunc2->SetParameters(10.0, 3.0, 1.5, 6.0, 10.0, 3.0 );
TH1F* myHist2 = new TH1F("hist2","Another function to be fit",100,0.,15.);
myHist2->FillRandom("func2",10000);

myHist2->GetXaxis()->SetTitle("Double Gaussian");
myHist2->SetMarkerColor(4);
myHist2->SetMarkerStyle(21);
myHist2->SetMarkerSize(0.7);                   
myHist2->SetOption("e1");
myHist2->Write();

myFile->Close();
}
