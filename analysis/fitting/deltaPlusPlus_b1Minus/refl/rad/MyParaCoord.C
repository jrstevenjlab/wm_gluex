/* Plotting several different parallel coordinates requires re-using
 * several commands, so this file contains my unique class for these
 * functions.
 *
 * In future, try inheriting this from TParallelCoord class so the 
 * "para_" member variable does not need to be called seperately
 * 
 * ROOT BUGS
 * 1. AddRange to an axis (in AddSelection function) will always force
 *    the first selection to be blue, and any subsequent selection to 
 *    have the color of the previous selection
 *      i.e. AddSelection(*,*,*,kViolet) AddSelection(*,*,*,kGreen) 
 *      will make the first be kBlue, and the second one be kViolet
 *      and kGreen is ignored, unless another selection is added
 *      This is known issue. See my forum post for example:
 *      https://root-forum.cern.ch/t/
 *      tparallelcoord-first-selection-forced-to-blue/53754
 * 2. Do NOT use SetCurrentMin/Max for changing range of axis, use
 *    SetCurrentLimits(min,max). Otherwise the lines will not adjust
 *    to match the histogram locations on each axis. This is used in
 *    the RoundRange function
 * 3. SetDotsSpacing(i) for i>0 causes the lines to the first axis
 *    near its limits to not connnect with the axis
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
  
  /* Rewrites the variable names with my custom filters depending
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
      // remove end labels
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

  /* Adds a colored selection range to variable "varName", from 
   * rangeStart to rangeEnd
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

  /* The limits of a variable  on a default ParaCoord are the exact 
   * decimal values for that variable, which is ugly. This function 
   * rounds those limits to a decimal place
   */
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
