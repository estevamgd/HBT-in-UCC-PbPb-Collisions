#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TH1D.h"
#include <TStyle.h>
#include "TLegend.h"
#include <TText.h>
#include <TBenchmark.h>
#include <TLine.h>
#include <TFitResult.h>

#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"
#include "../include/data_func.h"

void fitting_sr() {
    ROOT::EnableImplicitMT();
    gBenchmark->Start("hsimple");
    auto poolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << poolSize << std::endl;

    // Setting canvases
    TCanvas *c1 = new TCanvas("c1", "", 1920, 1080);

    TCanvas *canvases[] = {c1};
    int numCanvases = 1;

    // Open the ROOT file containing histograms
    TFile *fr = TFile::Open("data/50_70_qinv_normqinv_sr.root", "READ");
    if (!fr || fr->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve histograms
    TH1D *h1 = (TH1D *)fr->Get("h1");
    TH1D *h2 = (TH1D *)fr->Get("h2");
    TH1D *h1_normalized = (TH1D *)fr->Get("h1_normalized");
    TH1D *h2_normalized = (TH1D *)fr->Get("h2_normalized");
    TH1D *sr = (TH1D *)fr->Get("sr");

    if (!h1 || !h2 || !h1_normalized || !h2_normalized || !sr) {
        std::cerr << "Error: Histograms not found!" << std::endl;
        fr->Close();
        return;
    }

    TH1D *histograms[] = {h1, h2, h1_normalized, h2_normalized, sr};
    int numHistograms = 5;
    
    // Horizontal line from x1 to x2 at y
    double y = 1.;  // Adjust this to your needs
    TLine *line = new TLine(sr->GetXaxis()->GetXmin(), y, sr->GetXaxis()->GetXmax(), y);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(1);

    // Fixing some colors
    sr->SetLineWidth(1);
    sr->SetLineColor(kBlack);
    sr->SetMarkerStyle(4);
    sr->SetMarkerSize(0.8);

    // Setting fits
    TF1 *fit_exp = new TF1("fit_exp", func1_exp, 0.0, 1.0, 4);
    fit_exp->SetParameters(1.0, 1.0, 4.0, 0.0); 
    fit_exp->SetParName(0,"Const");
    fit_exp->SetParName(1,"#lambda");
    fit_exp->SetParName(2,"R (fm)");
    fit_exp->SetParName(3,"#epsilon");
    fit_exp->SetLineColor(kBlue); 
    fit_exp->SetLineWidth(2); 

    TF1 *fit_gauss = new TF1("f_gauss", func2_gauss, 0.0, 1.0, 4);
	fit_gauss->SetParameters(1.0, 1.0, 4.0, 0.0);
	fit_gauss->SetParName(0,"Const");
	fit_gauss->SetParName(1,"#lambda"); 
	fit_gauss->SetParName(2,"R (fm)");
	fit_gauss->SetParName(3,"#epsilon");
    fit_gauss->SetLineColor(kRed); 
	fit_gauss->SetLineWidth(2);

    TF1 *fit_levy = new TF1("f_levy", func3_levy, 0.0, 1.0, 5);
	fit_levy->SetParameters(1.0, 1.0, 4.0, 0.0, 4.0);
	fit_levy->SetParName(0,"Const");
	fit_levy->SetParName(1,"#lambda");
	fit_levy->SetParName(2,"R (fm)");
	fit_levy->SetParName(3,"#epsilon");
	fit_levy->SetParName(4,"#aplha");
    fit_levy->SetLineColor(kGreen); 
	fit_levy->SetLineWidth(2);

    // Fitting
    TFitResultPtr res_exp, res_gauss, res_levy;
    res_exp = sr->Fit(fit_exp, "S R");
    res_gauss = sr->Fit(fit_gauss, "S R");
    res_levy = sr->Fit(fit_levy, "S R");

    // Getting information
    Double_t chi2_exp = res_exp->Chi2();
    Int_t ndf_exp = res_exp->Ndf();

    Double_t chi2_gauss = res_exp->Chi2();
    Int_t ndf_gauss = res_exp->Ndf();

    Double_t chi2_levy = res_levy->Chi2();
    Int_t ndf_levy = res_levy->Ndf();

    // Calculating p-value
    Double_t p_value_exp = TMath::Prob(chi2_exp, ndf_exp);
    Double_t p_value_gauss = TMath::Prob(chi2_gauss, ndf_gauss);
    Double_t p_value_levy = TMath::Prob(chi2_levy, ndf_levy);

    // Adding labels
    sr->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV");
    sr->GetXaxis()->SetTitle("q_{inv} [GeV]");
    sr->GetYaxis()->SetTitle("Single Ratio SS/OS");

    // Adding legend
    TLegend *h1_legend = new TLegend(0.7, 0.4, 0.9, 0.9); 
    h1_legend->AddEntry((TObject*)0, "50-70%", "");
    h1_legend->AddEntry(sr, "Data");
    h1_legend->AddEntry(fit_exp, "Exponential Fit", "l");
    h1_legend->AddEntry((TObject*)0, Form("  R = %.2f #pm %.2f", fit_exp->GetParameter(2), fit_exp->GetParError(2)), "");
    h1_legend->AddEntry((TObject*)0, Form("  #lambda = %.2f #pm %.2f", fit_exp->GetParameter(1), fit_exp->GetParError(1)), "");
    h1_legend->AddEntry((TObject*)0, Form("  p-value = 0.0"), "");
    
    h1_legend->AddEntry(fit_gauss, "Gaussian Fit", "l");
    h1_legend->AddEntry((TObject*)0, Form("  R = %.2f #pm %.2f", fit_gauss->GetParameter(2), fit_gauss->GetParError(2)), "");
    h1_legend->AddEntry((TObject*)0, Form("  #lambda = %.2f #pm %.2f", fit_gauss->GetParameter(1), fit_gauss->GetParError(1)), "");
    h1_legend->AddEntry((TObject*)0, Form("  p-value = 0.0"), "");
    
    h1_legend->AddEntry(fit_levy, "Levy Fit", "l");
    h1_legend->AddEntry((TObject*)0, Form("  R = %.2f #pm %.2f", fit_levy->GetParameter(2), fit_levy->GetParError(2)), "");
    h1_legend->AddEntry((TObject*)0, Form("  #lambda = %.2f #pm %.2f", fit_levy->GetParameter(1), fit_levy->GetParError(1)), "");
    h1_legend->AddEntry((TObject*)0, Form("  p-value = 0.0"), "");

    // Setting y range to 0.95<y<1.6
    sr->GetYaxis()->SetRangeUser(0.95, 1.6);
    
    // Removing statistics box
    sr->SetStats(0);

    // Drawing the single ratio
    c1->cd(); sr->Draw(); 
    fit_exp->Draw("SameL"); 
    fit_gauss->Draw("SameL"); 
    fit_levy->Draw("SameL"); 
    line->Draw("same");
    h1_legend->Draw();

    ///*    // comment/uncomment to save/show in TBrowser the image
    // Saving image
    const char *path = "./imgs/final/";
    const char *prefix = "final-fitting-sr";
    const char *file_type = "png";
    save_canvas_images(canvases, numCanvases, path, prefix, file_type);
    
    // Closing program
    delete fit_exp;
    delete fit_gauss;
    delete fit_levy;
    close_program(canvases, numCanvases, histograms, numHistograms, fr);  
    //*/    // comment/uncomment to save/show in TBrowser the image

    gBenchmark->Show("hsimple");
}