from ROOT import TFile, gDirectory, TH2D, TCanvas

# Open the file. Note that the name of your file outside this class
# will probably NOT be experiment.root.

myfile = TFile( 'experiment.root' )

# Retrieve the n-tuple of interest. In this case, the n-tuple's name is
# "tree1". You may have to use the TBrowser to find the name of the
# n-tuple that someone gives you.
mychain = gDirectory.Get( 'tree1' )
entries = mychain.GetEntriesFast()

### The Set-up code goes here.
# Create a 2D histogram
myHist = TH2D("hist2D","chi2 vs ebeam",100,0,20,100,149,151)
myHist.GetXaxis().SetTitle("chi2")
myHist.GetYaxis().SetTitle("ebeam [GeV]")
###

for jentry in xrange( entries ):
   # Get the next tree in the chain and verify.
   ientry = mychain.LoadTree( jentry )
   if ientry < 0:
      break

   # Copy next entry into memory and verify.
   nb = mychain.GetEntry( jentry )
   if nb <= 0:
      continue

   # Use the values directly from the tree. This is an example using a
   # variable "vertex". This variables does not exist in the example
   # n-tuple experiment.root, to force you to think about what you're
   # doing.  
   # myValue = mychain.vertex 
   # myHist.Fill(myValue)

   ### The Loop code goes here.
   chi2 = mychain.chi2
   ebeam = mychain.ebeam
   myHist.Fill(chi2,ebeam)
   ###
   
### The Wrap-up code goes here
myCanvas = TCanvas()
myHist.Draw()
myCanvas.Draw()
###
