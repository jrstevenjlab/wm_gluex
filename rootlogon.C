// This is the file rootlogon.C                                                 
{
  printf("Load MyStyle $WM_GLUEX/rootlogon.C \n");

  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

  // no stats or title
  myStyle->SetOptStat(00000000);
  myStyle->SetOptTitle(0);

  // from ROOT plain style                                                      
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetStatColor(0);
  myStyle->SetPalette(1);

  // set axis labels
  myStyle->SetLabelSize(0.05,"xyz"); // size of axis value font
  myStyle->SetTitleSize(0.06,"xyz"); // size of axis title font
  myStyle->SetTitleOffset(1.05,"y");
  myStyle->SetTitleOffset(1.0,"x");

  // default canvas positioning                                                 
  //myStyle->SetCanvasDefX(900);
  //myStyle->SetCanvasDefY(20);
  //myStyle->SetCanvasDefH(550);
  //myStyle->SetCanvasDefW(540);
  myStyle->SetPadBottomMargin(0.13);
  myStyle->SetPadTopMargin(0.02);
  myStyle->SetPadLeftMargin(0.14);
  myStyle->SetPadRightMargin(0.1);

  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);

  myStyle->SetFrameBorderMode(0);

  gROOT->SetStyle("MyStyle"); //uncomment to set this style                        

  if(getenv("ROOT_ANALYSIS_HOME"))       
  	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/programs/MakePROOFPackage/SETUP.C");
}
