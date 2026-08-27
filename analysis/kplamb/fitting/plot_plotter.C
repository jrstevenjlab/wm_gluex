void plot_plotter(TString name = "khyperon_plot.root", TString dir = "./", TString reac = "") {

        int rebin = 2;

        gStyle->SetOptStat(0);

        const int maxPlots = 7;
        TString plotNames[maxPlots] = {"MHyperon","cosThetaHyp","phiHyp","cosThetaX","cosThetaY","cosThetaZ","Phi"};
        const int maxAmps = 1;
        TString ampNames[maxAmps] = {""};
        TString ampTitles[maxAmps] = {"Fit Result"};
        TString ampDrawOpt[maxAmps] = {"h"};
        int ampColors[maxAmps] = {1};

        TFile *f = TFile::Open(dir+"/"+name);
        f->ls();

        TCanvas *cc = new TCanvas("cc","cc",800,600);
        cc->Divide(3,3);
        TCanvas *ee = new TCanvas("ee","ee",800,600);
        ee->Divide(3,3);

        double textSize = 0.10;
        TLegend *leg1 = new TLegend(0.1, 0.1, 0.5, 0.9);
        leg1->SetEntrySeparation(0.01);
        leg1->SetNColumns(1);
        leg1->SetColumnSeparation(1.0);
        leg1->SetMargin(0.2);
        leg1->SetFillColor(0);
        leg1->SetTextSize(textSize);
        leg1->SetBorderSize(0);
        leg1->SetLineColor(kWhite);

        for(int i=0; i<maxPlots; i++) {
                //cout<<(reac+plotNames[i]).Data()<<endl;

                TH1F *hdat = (TH1F*)f->Get(reac+plotNames[i]+"dat");
                hdat->Rebin(rebin);
                TH1F *hbkg = (TH1F*)f->Get(reac+plotNames[i]+"bkgnd");
                hbkg->Rebin(rebin);

                hdat->SetLineColor(kBlack);
                hdat->SetMinimum(0);
                cc->cd(i+1);
                hdat->SetMarkerStyle(20);
                hdat->SetMarkerSize(0.5);

                if(i==0) {
                        leg1->AddEntry(hdat, "GlueX Data", "ep");
                        hdat->GetXaxis()->SetRangeUser(1.05,1.15);
                        hdat->Draw();
                }
                else {
                        hdat->Draw();
                }

                for(int j=0; j<maxAmps; j++) {
                    cout<<(reac+plotNames[i]+"acc"+ampNames[j]).Data()<<endl;
                    TH1F *hacc = (TH1F*)f->Get(reac+plotNames[i]+"acc"+ampNames[j]);
                    TH1F *hgen = (TH1F*)f->Get(reac+plotNames[i]+"gen"+ampNames[j]);
                    hacc->Rebin(rebin);
                    hgen->Rebin(rebin);

                    TH1F *heff = (TH1F*)hacc->Clone(reac+plotNames[i]+"eff");
                    heff->SetMarkerStyle(20);
                    heff->SetMarkerSize(0.5);
                    heff->SetMinimum(0);
                    heff->Divide(hgen);

                    if(j==0) {
                        hacc->SetFillColor(31);
                        hacc->SetLineColor(31);
                        if(i==0) {
                                leg1->AddEntry(hacc, ampTitles[j], "f");
                                hacc->Draw("same" + ampDrawOpt[j]);
                                //hbkg->Draw("same" + ampDrawOpt[j]);
                        }
                        else {
                                hacc->Draw("same" + ampDrawOpt[j]);
                                //hbkg->Draw("same" + ampDrawOpt[j]);
                        }

                        ee->cd(i+1);
                        heff->Draw("ep");

                        cc->cd(i+1);
                    }
                    else {
                        hacc->SetLineColor(ampColors[j]);
                        hacc->SetMarkerColor(ampColors[j]);
                        hacc->SetMarkerSize(0.5);
                        hacc->SetMarkerStyle(20);
                        if(i==0) {
                                leg1->AddEntry(hacc, ampTitles[j], ampDrawOpt[j]);
                                hacc->Draw("same" + ampDrawOpt[j]);
                        }
                        else hacc->Draw("same" + ampDrawOpt[j]);
                    }
                }

                // re-draw data and background on top of fit result
                hdat->Draw("same");

                hbkg->SetFillColor(kRed+2);
                hbkg->SetLineColor(kRed+2);
                if(i==0) leg1->AddEntry(hbkg, "Non-#Lambda Sideband", "f");
                hbkg->Draw("same h");
        }

        cc->cd(9);
        leg1->Draw();

        cc->Print(dir+"/fit.pdf");
        ee->Print(dir+"/fit_eff.pdf");

        return;
}

