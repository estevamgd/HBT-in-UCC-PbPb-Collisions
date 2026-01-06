#ifndef FITSINGLERATIO_H
#define FITSINGLERATIO_H
#include <iostream>
#include <string>
#include <vector>
#include <filesystem> 
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h" 
#include "TSystem.h" 
#include "../include/my_func.h" 
#include "../include/data_func.h"
#include "Math/MinimizerOptions.h"


void fitSingleRatio(const char* inputFileName, const char* histName,
    const char* outputPrefix,
    const char* headerLabel = "N/A",
    double fitMin = 0.0, double fitMax = 1.0,
    double plotXMin = 0.0, double plotXMax = 0.6){
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
std::cout << "Default minimizer: "
          << ROOT::Math::MinimizerOptions::DefaultMinimizerType()
          << std::endl;
          
ROOT::EnableImplicitMT();
gStyle->SetOptStat(0);     
gStyle->SetPalette(kColorPrintableOnGrey);

TFile *inputFile = TFile::Open(inputFileName, "READ");
if (!inputFile || inputFile->IsZombie()) {
std::cerr << "Error: Could not open input file '" << inputFileName << "'!" << std::endl;
return;
}

TH1D *histToFit = (TH1D *)inputFile->Get(histName);
if (!histToFit) {
std::cerr << "Error: Histogram '" << histName << "' not found in file '" << inputFileName << "'!" << std::endl;
inputFile->Close();
delete inputFile;
return;
}

histToFit = (TH1D*)histToFit->Clone("histToFit_clone");

TString xTitle = "q_{inv} (GeV)";
if (TString(outputPrefix).Contains("qlcms"))
    xTitle = "q_{LCMS} (GeV)";

TCanvas *c1 = new TCanvas("c_fit_canvas", "Correlation Fit", 1200, 900);
c1->SetLeftMargin(0.12);
c1->SetBottomMargin(0.12);

histToFit->GetXaxis()->SetTitle(xTitle);
histToFit->GetYaxis()->SetTitle("C(q) = SS / OS");
histToFit->GetXaxis()->SetTitleSize(0.045);
histToFit->GetYaxis()->SetTitleSize(0.045);

histToFit->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
histToFit->GetYaxis()->SetRangeUser(0.9, 2.1);

histToFit->SetMarkerStyle(20);
histToFit->SetMarkerSize(1.0);
histToFit->SetMarkerColor(kBlack);
histToFit->SetLineColor(kBlack);

// Fits
TF1 *fitExp = new TF1("fitExp", FitExp, fitMin, fitMax, 4);
fitExp->SetParNames("Const", "#lambda", "R (fm)", "#epsilon");
fitExp->SetParameters(1.0, 0.5, 5.0, 0.0);
fitExp->SetLineColor(gStyle->GetColorPalette(50));
fitExp->SetLineWidth(2);

TF1 *fitGauss = new TF1("fitGauss", FitGauss, fitMin, fitMax, 4);
fitGauss->SetParNames("Const", "#lambda", "R (fm)", "#epsilon");
fitGauss->SetParameters(1.0, 0.5, 5.0, 0.0);
fitGauss->SetLineColor(gStyle->GetColorPalette(150));
fitGauss->SetLineWidth(2);

TF1 *fitLevy = new TF1("fitLevy", FitLevy, fitMin, fitMax, 5);
fitLevy->SetParNames("Const", "#lambda", "R (fm)", "#epsilon", "#alpha");
fitLevy->SetParameters(1.0, 0.5, 5.0, 0.0, 1.2);
fitLevy->SetParLimits(1, 0., 1.0);
fitLevy->SetParLimits(4, 1., 2.0);
fitLevy->SetLineColor(gStyle->GetColorPalette(220));
fitLevy->SetLineWidth(2);

TF1 *fitLevy2 = new TF1("fitLevy2", FitLevy2, fitMin, fitMax, 3);
fitLevy2->SetParNames("#lambda", "R (fm)",  "#alpha");
fitLevy2->SetParameters(0.5, 5.0, 1.2);
fitLevy2->SetParLimits(0, 0.0, 1.0);
fitLevy2->SetParLimits(2, 1.0, 2.0);
fitLevy2->SetLineColor(gStyle->GetColorPalette(20));
fitLevy2->SetLineWidth(2);

//TFitResultPtr resExp   = histToFit->Fit(fitExp,   "S R E M");
//TFitResultPtr resGauss = histToFit->Fit(fitGauss, "S R E M +");
TFitResultPtr resLevy  = histToFit->Fit(fitLevy,  "S R E M +");
TFitResultPtr resLevy2 = histToFit->Fit(fitLevy2, "S R E M +");

histToFit->GetXaxis()->SetRangeUser(plotXMin, plotXMax); 
histToFit->Draw("E1 P");

TLine *line = new TLine(plotXMin, 1.0, plotXMax, 1.0); 

line->SetLineColor(kGray + 2);
line->SetLineStyle(kDashed);
line->SetLineWidth(2);
line->Draw("SAME");

TLegend *legend = new TLegend(0.55, 0.40, 0.88, 0.88);
legend->SetTextFont(42);
legend->SetTextSize(0.028);
legend->SetBorderSize(0);
legend->SetFillStyle(0);
legend->AddEntry(histToFit, "Data", "lep");

//legend->AddEntry(fitExp, "Exponential Fit", "l");
//legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitExp->GetParameter(0), fitExp->GetParError(0)), "");
//legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitExp->GetParameter(1), fitExp->GetParError(1)), "");
//legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resExp->Chi2(), resExp->Ndf()), "");
//legend->AddEntry((TObject*)0, Form(" p-value = %.5f", resExp->Prob()), "");
//
//legend->AddEntry(fitGauss, "Gaussian Fit", "l");
//legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitGauss->GetParameter(0), fitGauss->GetParError(0)), "");
//legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitGauss->GetParameter(1), fitGauss->GetParError(1)), "");
//legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resGauss->Chi2(), resGauss->Ndf()), "");
//legend->AddEntry((TObject*)0, Form(" p-value = %.5f", resGauss->Prob()), "");

legend->AddEntry(fitLevy, "Levy Fit", "l");
legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitLevy->GetParameter(1), fitLevy->GetParError(1)), "");
legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitLevy->GetParameter(2), fitLevy->GetParError(2)), "");
legend->AddEntry((TObject*)0, Form(" #alpha = %.2f #pm %.2f", fitLevy->GetParameter(4), fitLevy->GetParError(4)), "");
legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resLevy->Chi2(), resLevy->Ndf()), "");
legend->AddEntry((TObject*)0, Form(" p-value = %.5f", resLevy->Prob()), "");

legend->AddEntry(fitLevy2, "Levy Fit", "l");
legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitLevy2->GetParameter(0), fitLevy2->GetParError(0)), "");
legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitLevy2->GetParameter(1), fitLevy2->GetParError(1)), "");
legend->AddEntry((TObject*)0, Form(" #alpha = %.2f #pm %.2f", fitLevy2->GetParameter(2), fitLevy2->GetParError(2)), "");
legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resLevy2->Chi2(), resLevy2->Ndf()), "");
legend->AddEntry((TObject*)0, Form(" p-value = %.5f", resLevy2->Prob()), "");

legend->Draw();

drawCMSHeaders("#bf{CMS} #it{Preliminary}", headerLabel);

const char* imagePath = "./imgs/test/fit_correlation/";
const char* histoPath = "./data/fit_correlation/";

TCanvas *canvases[] = {c1};
TH1D *histograms[] = {histToFit};
TLegend *legends[] = {legend};

save_histograms(histograms, 1, histoPath, outputPrefix);
save_canvas_images(canvases, 1, imagePath, outputPrefix, "png");
save_canvas_images(canvases, 1, imagePath, outputPrefix, "pdf");

close_program(canvases, 1, histograms, 1, legends, 1, inputFile);

delete fitExp;
delete fitGauss;
delete fitLevy;
delete fitLevy2;
delete line;
}

void processRatio(const std::string& searchPath, 
    const char* observable,   // qinv | qlcms
    const char* controlVar,
    double selMin,
    double selMax,
    const char* typeTag,      // Cor | Uncor
    const char* histName,
    double fitMin,
    double fitMax,
    double plotXMin,
    double plotXMax){
std::cout << "\n--- Processing " << observable << " (" << typeTag << ") ---\n";

TString filePattern = TString::Format(
    "%s_SingleRatio_%s_%s_%.0f-%.0f*.root",
    observable, typeTag, controlVar, selMin, selMax
);

TString inputFile = findFile(searchPath + filePattern.Data());

if (inputFile.IsNull()) {
    std::cerr << "File not found: " << filePattern << std::endl;
    return;
}

TString outputPrefix = TString::Format(
    "fit_%s_%s_%s_%.0f-%.0f",
    observable, typeTag, controlVar, selMin, selMax
);

TString header = TString::Format(
    "%s: %.0f-%.0f (%s)",
    controlVar, selMin, selMax, typeTag
);

fitSingleRatio(inputFile.Data(), histName,
               outputPrefix.Data(), header.Data(),
               fitMin, fitMax, plotXMin, plotXMax);
}

#endif