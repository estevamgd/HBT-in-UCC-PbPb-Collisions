#ifndef RATIOS_H
#define RATIOS_H

#include <iostream>
#include <string>
#include <filesystem>
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"  
#include "TString.h" 
#include "TSystem.h" 
#include "../include/my_func.h"
#include "../include/normalizer.h"

void saveRatio(TH1D* hRatio,
                    const char* outputPrefix,
                    double plotXMin = 0.0, double plotXMax = 10.0,
                    double plotYMin = 0.9, double plotYMax = 2.1,
                    const char* outputHistName = "sr_cor",
                    const char* plotHeaderText = "PbPb 2.76 TeV | Single Ratio",
                    const char* histTittle = "; q_{inv} [GeV]; C(q_{inv}) = SS/OS",
                    const char* legend = "Single Ratio (SS/OS)") {
    TH1D* hRatioClone = (TH1D*)(hRatio->Clone("hRatioClone"));

    gStyle->SetOptStat(0);     
    hRatioClone->SetTitle(histTittle);
    hRatioClone->SetMarkerStyle(20);
    hRatioClone->SetMarkerSize(0.8);
    hRatioClone->SetMarkerColor(kBlack);
    hRatioClone->SetLineColor(kBlack);
    std::cout << "Ratio histogram '" << outputHistName << "' created successfully." << std::endl;

    TCanvas *cRatio = new TCanvas("cRatio", "Ratio", 1200, 800);
    cRatio->cd();
    cRatio->SetLeftMargin(0.12);
    cRatio->SetBottomMargin(0.12);
    hRatioClone->GetYaxis()->SetTitleOffset(1.2);
    hRatioClone->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
    hRatioClone->GetYaxis()->SetRangeUser(plotYMin, plotYMax);
    hRatioClone->Draw("E1 PLC");
    
    TLegend *legendRatio = new TLegend(0.6, 0.8, 0.88, 0.88);
    legendRatio->AddEntry(hRatioClone, legend, "lep");
    legendRatio->SetFillStyle(0);
    legendRatio->SetBorderSize(0);
    legendRatio->Draw();

    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", plotHeaderText);

    const char* histoPath = "./data/correlation_ratios/";
    TH1D *histogramsToSave[] = {hRatioClone};

    save_histograms(histogramsToSave, 1, histoPath, outputPrefix);
    std::cout << "Output histograms saved to a new ROOT file in " << histoPath << std::endl;

    const char* imagePath = "./imgs/test/correlation_ratios/";
    TCanvas *canvasesToSave[] = { cRatio };
    save_canvas_images(canvasesToSave, 1, imagePath, outputPrefix, "png");
    save_canvas_images(canvasesToSave, 1, imagePath, outputPrefix, "pdf");
    std::cout << "Output plots saved to " << imagePath << std::endl;
    
    delete hRatioClone;
    delete legendRatio;
    delete cRatio;
}

TH1D* histhistRatio(TH1D* hist1, TH1D* hist2,
                    Double_t q1, Double_t q2,
                    const char* outputHistName) {
    if (!hist1 || !hist2) {
        std::cerr << "histhistRatio ERROR: null input histogram" << std::endl;
        return nullptr;
    }

    if (hist1->GetNbinsX() != hist2->GetNbinsX()) {
        std::cerr << "histhistRatio ERROR: different binning" << std::endl;
        return nullptr;
    }

    TH1D* hist1_norm = (TH1D*)(hist1->Clone("hist1_normalized"));
    TH1D* hist2_norm = (TH1D*)(hist2->Clone("hist2_normalized"));

    TH1D* toNorm[] = {hist1_norm, hist2_norm};
    normalizer(toNorm, 2, q1, q2, 1.0);

    TH1D* hRatio = (TH1D*)(hist1_norm->Clone(outputHistName));

    hRatio->Divide(hist2_norm);

    delete hist1_norm;
    delete hist2_norm;

    return hRatio;
}


TH1D* histfuncRatio(TH1D* histogram, TF1* function,
                    const char* outputHistName) {
    if (!histogram || !function) {
        std::cerr << "histfuncRatio ERROR: null input" << std::endl;
        return nullptr;
    }

    TH1D* hRatio = dynamic_cast<TH1D*>(
        histogram->Clone(outputHistName)
    );
    hRatio->SetDirectory(nullptr);

    hRatio->Divide(function);

    return hRatio;
}

#endif