/* Plotting several different parallel coordinates requires re-using
 * several commands, so this file contains my unique class for these
 * funcitons.
 *
 * In future, may split this into header/implementation file if
 * paraCoords see more use beyond RadPlotter
 * 
 * ROOT BUGS
 * 1. AddRange to an axis (in AddSelection function) will always force
 *    the first selection to be blue, and any subsequent selection to 
 *    have the color of the previous selection
 *    i.e. AddSelection(*,*,*,kViolet) AddSelection(*,*,*,kGreen) will
 *    make the first one be kBlue, and the second one be kViolet
 * 2. Do NOT use SetCurrentMin/Max for changing range of axis, use
 *    SetCurrentLimits(min,max). This avoids issue where lines do not
 *    match the histogram locations on each axis.
 * 3. SetDotsSpacing(i) for i>0 causes lines on the first axis near
 *    its limits to not connnect with the axis
 */

class MyParaCoord {
public:
  TParallelCoord* para_;
  MyParaCoord() {
    // "ParaCoord" is recognized by ROOT to grab a parallel coordinate
    // object on the gPad
    para_ = (TParallelCoord*)gPad->
      GetListOfPrimitives()->FindObject("ParaCoord");
  }
  MyParaCoord(TParallelCoord* p) : para_(p) {}
  
  /* Rewrites the variable names with some custom filters depending
   on the name. Avoids overlapping, and makes the phase differnces
   understandable.
  */
  void CleanupAxesNames(std::map<TString, 
			std::map<TString,int>> phaseMap = {}) {
    std::string name;
    std::string refAmp;
    int start, end;
    double min, max;

    TList *varList = para_->GetVarList();
    for(TObject *var : *varList) {
      TParallelCoordVar *paraVar = (TParallelCoordVar*)var;    
      name = paraVar->GetName();      
      // remove labels (phase difference case handled below)
      if(name.find("->") != std::string::npos &&
	 name.find("[") == std::string::npos) {
	end = name.find("->");
	name = name.substr(0,end);	
	paraVar->SetName(name.c_str());
      }
      // rename phase differences to be clearer
      if(name.find("[") != std::string::npos &&
	 phaseMap.size() != 0) {           
	refAmp = name.substr(0,5); // get reference amp name
	
	// grab the index and make it an int
	start = name.find("[") + 1;
	end = name.find("]") - 1;
	int i = std::stoi(name.substr(start,end));
	
	// rename axis to phi_amp1 \n amp2
	for(auto const& itrOut : phaseMap) {
	  if(itrOut.first != refAmp) continue;
	  for(auto const& itrIn : itrOut.second) {
	    if(itrIn.second == i) 
	      name = "#splitline{#phi_" + refAmp + "}{"
		+ itrIn.first + "}";
	  }
	}	
	paraVar->SetName(name.c_str());
      } 
    }    
    return;
  }

  /* Adds a selection range to variable "varName", from rangeStart to
   * rangeEnd
   */
  void AddSelection(TString varName, double rangeStart, 
		    double rangeEnd, Color_t myColor) {
    TParallelCoordVar* axis = (TParallelCoordVar*)para_->
      GetVarList()->FindObject(varName);
    axis->AddRange(new TParallelCoordRange(axis, rangeStart, rangeEnd));
    para_->AddSelection("");
    para_->GetCurrentSelection()->SetLineColor(myColor);
    return;
  }

  void RoundRange(int decimal) {
    TList *varList = para_->GetVarList();
    double multiplier = std::pow(10.0, decimal);
    double min, max;
    double minRound, maxRound;
    for(TObject *var : *varList) {
      TParallelCoordVar *paraVar = (TParallelCoordVar*)var;    
      min = paraVar->GetCurrentMin();
      max = paraVar->GetCurrentMax();
      
      minRound = std::floor(min*multiplier) / multiplier;
      maxRound = std::ceil(max*multiplier) / multiplier;

      paraVar->SetCurrentLimits(minRound, maxRound);
    }
  }
};
