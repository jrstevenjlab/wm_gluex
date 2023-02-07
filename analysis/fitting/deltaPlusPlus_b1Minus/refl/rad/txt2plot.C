#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <bits/stdc++.h>
#include <string>

#include "TH1I.h"
#include "TH2I.h"

void txt2plot(){

	const int nBins = 6;
	double edges[nBins+1];
	for(int i = 0; i < nBins+1; i++){
		edges[i] = 1.0 + 0.1*i;
	}

	const int maxWaves = 4;
	std::string waveNames[maxWaves] = {"0m", "1p", "1m", "2m"};

	

	TH1F* hWavePos[maxWaves];
	TH1F* hWaveNeg[maxWaves];
	TH1F* hTotalIntensity = new TH1F("hTotalIntensity", ";M_{#omega#pi^{-}} (GeV);Intensity", nBins, 1.0, 1.6);
	hTotalIntensity->SetLineColor(kBlack);
	hTotalIntensity->SetMarkerColor(kBlack);
	hTotalIntensity->SetMarkerStyle(20);
	hTotalIntensity->SetMinimum(0);

	TH1F* h1p1mPhaseDiff = new TH1F("h1p1mPhaseDiff", ";M_{#omega#pi^{-}} (GeV);Phase Diff", nBins, 1.0, 1.6);
	string phaseDiffLead = "PHASE DIFF omegapi::ImagPosSign::1pps omegapi::ImagPosSign::1mmp ";

	for(int j = 0; j < maxWaves; j++){
		hWavePos[j] = new TH1F(Form("hWavePos_%s", waveNames[j].c_str()), ";M_{#omega#pi^{-}} (GeV);Fit Fraction", nBins, 1.0, 1.6);
		hWaveNeg[j] = new TH1F(Form("hWaveNeg_%s", waveNames[j].c_str()), ";M_{#omega#pi^{-}} (GeV);Fit Fraction", nBins, 1.0, 1.6);
		cout << "creating histogram " << j << endl;
	}

	for(int i = 0; i < nBins; i++){
		string spath = Form("/volatile/halld/home/aschertz/ampToolsFits/deltaPlusPlus_b1Minus/allPeriods/PARA_0/gpu_refl_1p1m/Mass_%0.3f_%0.3f/t_0.15_0.30/omegapi_fitPars.txt", edges[i], edges[i+1]);
		
		cout << spath << endl;	

		std::string delimiter = " +- ";

		string input3;
		ifstream ftext3(spath);

		while(getline(ftext3, input3)){
			if(input3.find(phaseDiffLead) != std::string::npos){
				cout << input3 << endl;
					
				size_t pos = 0;
				std::string token[2];
				double param[2];
				input3.erase(0, 65);	
				while((pos = input3.find(delimiter)) != std::string::npos){
					token[0] = Form("%s", input3.substr(0, pos).c_str());
					param[0] = atof(token[0].c_str());
					if(param[0] < -3.14)
						param[0] += 2*TMath::Pi();
					else if(param[0] > 3.14)
						param[0] += -2*TMath::Pi();
					input3.erase(0, pos + delimiter.length());
					token[1] = Form("%s", input3.c_str());
					param[1] = atof(token[1].c_str());
				}
				h1p1mPhaseDiff->SetBinContent(i+1, param[0]);
				h1p1mPhaseDiff->SetBinError(i+1, param[1]);
			}
		}

		for(int j = 0; j < maxWaves; j++){
			//cout << "j = " << j << endl; 
			//if(!hWavePos[j]) continue;
			double intensity;
			string inten = "TOTAL EVENTS = ";
			string reflP = Form("FIT FRACTION (coherent sum) PosRefl %s = ", waveNames[j].c_str());
			string reflN = Form("FIT FRACTION (coherent sum) NegRefl %s = ", waveNames[j].c_str());
			
			string input;
			string input1;
			string input2;
			ifstream ftext1(spath);
			ifstream ftext2(spath);
			ifstream ftext(spath);
			if(!ftext) continue;
		
			while(getline(ftext1, input1)){
				if(input1.find(inten) != std::string::npos){
			//		cout << input1 << endl;
					size_t pos = 0;
					std::string temp;
					input1.erase(0, 15);
					while((pos = input1.find(delimiter)) != std::string::npos){
						temp = Form("%s", input1.substr(0, pos).c_str());
						cout << temp << endl;
						intensity = atof(temp.c_str());
						input1.erase(0, pos + delimiter.length());
					}
			//		cout << inten << intensity << endl;
					hTotalIntensity->SetBinContent(i+1, intensity);
				}
			}
			while(getline(ftext, input)){
				if(input.find(reflP) != std::string::npos){	
			//		cout << input << endl;
					 
					size_t pos = 0;
					std::string token[2];
					double param[2];
					input.erase(0, 41);
					while((pos = input.find(delimiter)) != std::string::npos){
						token[0] = Form("%s", input.substr(0, pos).c_str());
						param[0] = atof(token[0].c_str());
						input.erase(0, pos + delimiter.length());
						token[1] = Form("%s", input.c_str());
						param[1] = atof(token[1].c_str());
					}
			
			//		cout << "PosRefl " << waveNames[j] << " = " << param[0] << endl;
			//		cout << "+- " << param[1] << endl;
					hWavePos[j]->SetBinContent(i+1, param[0]*intensity);
					hWavePos[j]->SetBinError(i+1, param[1]*intensity);
				}
			}
			while(getline(ftext2, input2)){
				if(input2.find(reflN) != std::string::npos){	
			//		cout << input2 << endl;
					 
					size_t pos = 0;
					std::string token[2];
					double param[2];
					input2.erase(0, 41);
					while((pos = input2.find(delimiter)) != std::string::npos){
						token[0] = Form("%s", input2.substr(0, pos).c_str());
						param[0] = atof(token[0].c_str());
						input2.erase(0, pos + delimiter.length());
						token[1] = Form("%s", input2.c_str());
						param[1] = atof(token[1].c_str());
					}
			
		//			cout << "NegRefl " << waveNames[j] << " = " << param[0] << endl;
		//			cout << "+- " << param[1] << endl;
					hWaveNeg[j]->SetBinContent(i+1, param[0]*intensity);
					hWaveNeg[j]->SetBinError(i+1, param[1]*intensity);
				}
			}
		}
	}


	TCanvas* c1 = new TCanvas("c1", "c1", 1);
	c1->Divide(2,2);

	for(int j = 0; j < maxWaves; j++){
		c1->cd(j+1);
		hTotalIntensity->Draw();
		hWavePos[j]->SetMinimum(0);
		hWavePos[j]->Draw("same");
		hWavePos[j]->SetLineColor(kRed);
		hWavePos[j]->SetMarkerColor(kRed);
		hWavePos[j]->SetMarkerStyle(20);
		hWaveNeg[j]->SetMarkerStyle(kOpenCircle);
		hWaveNeg[j]->SetMarkerColor(kBlue);
		hWaveNeg[j]->SetLineColor(kBlue);
		hWaveNeg[j]->Draw("same");
	}

	TCanvas* c2 = new TCanvas("c2", "c2", 1);
	c2->cd(1);
	h1p1mPhaseDiff->Draw();

	return;

}
