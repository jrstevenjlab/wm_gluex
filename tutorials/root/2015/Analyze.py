from ROOT import TFile, gDirectory
# You probably also want to import TH1D and TCanvas,
# unless you're not making any histograms.
from ROOT import TH1D, TCanvas

# Open the file. Note that the name of your file outside this class
# will probably NOT be experiment.root.

myfile = TFile( 'experiment.root' )

# Retrieve the n-tuple of interest. In this case, the n-tuple's name is
# "tree1". You may have to use the TBrowser to find the name of the
# n-tuple that someone gives you.
mychain = gDirectory.Get( 'tree1' )
entries = mychain.GetEntriesFast()

### The Set-up code goes here.
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
   ###
   
### The Wrap-up code goes here
###
