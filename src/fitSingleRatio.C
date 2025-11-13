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

#include "../include/data_func.h"
#include "../include/my_func.h"

/**
 * @brief Opens a ROOT file, fits a specified histogram, plots the results,
 * and saves all outputs using custom helper functions and predefined paths.
 *
 * @param inputFileName The path to the input ROOT file.
 * @param histName The name of the TH1D histogram to fit.
 * @param outputPrefix The base name for all output files (e.g., "fit_sr_50-70").
 * @param centralityLabel A string to display in the legend (e.g., "50-70%").
 * @param fitMin The lower bound of the fit range for q_inv in GeV.
 * @param fitMax The upper bound of the fit range for q_inv in GeV.
 */
void fitSingleRatio(const char* inputFileName, const char* histName,
                               const char* outputPrefix,
                               const char* centralityLabel = "N/A",
                               double fitMin = 0.0, double fitMax = 1.0) {

    ROOT::EnableImplicitMT();
    auto poolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << poolSize << std::endl;

    gStyle->SetOptStat(0);      // Disable statistics box
    gStyle->SetPalette(kLake);

    // 1. Open the ROOT file
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

    // Clone to avoid modifying the original
    histToFit = (TH1D*)histToFit->Clone("histToFit_clone");

    // 2. Create a canvas and style the histogram
    TCanvas *c1 = new TCanvas("c_fit_canvas", "Correlation Fit", 1200, 900);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    histToFit->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV");
    histToFit->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
    histToFit->GetYaxis()->SetTitle("C(q_{inv}) = (Same-Sign) / (Opposite-Sign)");
    histToFit->GetXaxis()->SetTitleSize(0.045);
    histToFit->GetYaxis()->SetTitleSize(0.045);
    histToFit->GetXaxis()->SetRangeUser(0.0, 0.6); // Visual range
    histToFit->GetYaxis()->SetRangeUser(0.9, 2.1); // Visual range

    histToFit->SetMarkerStyle(20);
    histToFit->SetMarkerSize(1.0);
    histToFit->SetMarkerColor(kBlack);
    histToFit->SetLineColor(kBlack);

    // 3. Define the fitting functions
    TF1 *fitExp = new TF1("fitExp", func1_exp, fitMin, fitMax, 4);
    fitExp->SetParNames("Const", "#lambda", "R (fm)", "#epsilon");
    fitExp->SetParameters(1.0, 0.5, 5.0, 0.0);
    fitExp->SetLineColor(gStyle->GetColorPalette(50));
    fitExp->SetLineWidth(2);

    TF1 *fitGauss = new TF1("fitGauss", func2_gauss, fitMin, fitMax, 4);
    fitGauss->SetParNames("Const", "#lambda", "R (fm)", "#epsilon");
    fitGauss->SetParameters(1.0, 0.5, 5.0, 0.0);
    fitGauss->SetLineColor(gStyle->GetColorPalette(150));
    fitGauss->SetLineWidth(2);

    TF1 *fitLevy = new TF1("fitLevy", func3_levy, fitMin, fitMax, 5);
    fitLevy->SetParNames("Const", "#lambda", "R (fm)", "#epsilon", "#alpha");
    fitLevy->SetParameters(1.0, 0.5, 5.0, 0.0, 1.2);
    fitLevy->SetParLimits(4, 0.5, 2.0);
    fitLevy->SetLineColor(gStyle->GetColorPalette(220));
    fitLevy->SetLineWidth(2);

    // 4. Perform the fits
    TFitResultPtr resExp   = histToFit->Fit(fitExp,   "SREM");
    TFitResultPtr resGauss = histToFit->Fit(fitGauss, "SREM+");
    TFitResultPtr resLevy  = histToFit->Fit(fitLevy,  "SREM+");

    // 5. Draw everything on the canvas
    histToFit->Draw("E1 P");
    TLine *line = new TLine(0.0, 1.0, 0.6, 1.0);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line->Draw("SAME");

    // 6. Create a legend
    TLegend *legend = new TLegend(0.55, 0.40, 0.88, 0.88);
    legend->SetTextFont(42);
    legend->SetTextSize(0.028);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetHeader(Form("Centrality: %s", centralityLabel), "C");
    legend->AddEntry(histToFit, "Data", "lep");

    legend->AddEntry(fitExp, "Exponential Fit", "l");
    legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitExp->GetParameter(1), fitExp->GetParError(1)), "");
    legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitExp->GetParameter(2), fitExp->GetParError(2)), "");
    legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resExp->Chi2(), resExp->Ndf()), "");
    
    legend->AddEntry(fitGauss, "Gaussian Fit", "l");
    legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitGauss->GetParameter(1), fitGauss->GetParError(1)), "");
    legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitGauss->GetParameter(2), fitGauss->GetParError(2)), "");
    legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resGauss->Chi2(), resGauss->Ndf()), "");
    
    legend->AddEntry(fitLevy, "L#acute{e}vy Fit", "l");
    legend->AddEntry((TObject*)0, Form(" #lambda = %.2f #pm %.2f", fitLevy->GetParameter(1), fitLevy->GetParError(1)), "");
    legend->AddEntry((TObject*)0, Form(" R = %.2f #pm %.2f fm", fitLevy->GetParameter(2), fitLevy->GetParError(2)), "");
    legend->AddEntry((TObject*)0, Form(" #alpha = %.2f #pm %.2f", fitLevy->GetParameter(4), fitLevy->GetParError(4)), "");
    legend->AddEntry((TObject*)0, Form(" #chi^{2}/NDF = %.1f / %d", resLevy->Chi2(), resLevy->Ndf()), "");
    
    legend->Draw();

    // 7. Save the outputs using your helper functions and predefined paths
    
    // Define standard output paths
    const char* imagePath     = "./imgs/test/fit_correlation/";
    const char* histoPath     = "./data/fit_correlation/";

    // Prepare object arrays for the helper functions
    TCanvas *canvases[] = {c1};
    TH1D *histograms[] = {histToFit};

    // Save histogram to a .root file
    save_histograms(histograms, 1, histoPath, outputPrefix);

    // Save canvas images
    save_canvas_images(canvases, 1, imagePath, outputPrefix, "png");
    save_canvas_images(canvases, 1, imagePath, outputPrefix, "pdf");

    // 8. Clean up memory using your helper function
    close_program(canvases, 1, histograms, 1, inputFile);
}

void execute_fitting() {
    // --- Configuration ---
    // 1. Input file and histogram name
    const char* inputFile = "data/50_70_qinv_normqinv_sr.root";
    const char* histName = "sr_cor";

    // 2. Output prefix (this is the base name for all output files)
    const char* outputPrefix  = "fit_sr_50-70";

    // 3. Legend label
    const char* centralityLabel = "50-70%";

    // 4. Fit range
    double fitMin = 0.0;
    double fitMax = 1.0;

    // --- Run the Analysis ---
    fitSingleRatio(inputFile, histName, outputPrefix,
                              centralityLabel, fitMin, fitMax);

    std::cout << "\nAnalysis finished for " << histName << "." << std::endl;
}