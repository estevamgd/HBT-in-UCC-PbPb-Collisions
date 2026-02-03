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
#include "TLine.h"  
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

void saveRatio(TH2D* hRatio,
                    const char* outputPrefix,
                    double plotXMin = 0.0, double plotXMax = 10.0,
                    double plotYMin = 0.9, double plotYMax = 2.1,
                    double plotZMin = 0.6, double plotZMax = 1.4,
                    const char* outputHistName = "sr_cor",
                    const char* plotHeaderText = "PbPb 2.76 TeV | Single Ratio",
                    const char* histTittle = "; q_{inv} [GeV]; C(q_{inv}) = SS/OS",
                    const char* legend = "Single Ratio (SS/OS)") {
    TH2D* hRatioClone = (TH2D*)(hRatio->Clone("hRatioClone"));
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);     
    hRatioClone->SetTitle(histTittle);
    hRatioClone->SetMarkerStyle(20);
    hRatioClone->SetMarkerSize(0.8);
    hRatioClone->SetMarkerColor(kBlack);
    hRatioClone->SetLineColor(kBlack);
    hRatioClone->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
    hRatioClone->GetYaxis()->SetRangeUser(plotYMin, plotYMax);
    hRatioClone->GetZaxis()->SetRangeUser(plotZMin, plotZMax);
    std::cout << "Ratio histogram '" << outputHistName << "' created successfully." << std::endl;

    TCanvas *cRatio = new TCanvas("cRatio", "Ratio", 1200, 800);
    cRatio->cd();
    //cRatio->SetLeftMargin(0.12);
    //cRatio->SetBottomMargin(0.12);
    //hRatioClone->GetYaxis()->SetTitleOffset(1.2);
    hRatioClone->Draw("COLZ");
    
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", plotHeaderText);

    const char* histoPath = "./data/correlation_ratios/";
    TH2D *histogramsToSave[] = {hRatioClone};

    save_histograms(histogramsToSave, 1, histoPath, outputPrefix);
    std::cout << "Output histograms saved to a new ROOT file in " << histoPath << std::endl;

    const char* imagePath = "./imgs/test/correlation_ratios/";
    TCanvas *canvasesToSave[] = { cRatio };
    save_canvas_images(canvasesToSave, 1, imagePath, outputPrefix, "png");
    save_canvas_images(canvasesToSave, 1, imagePath, outputPrefix, "pdf");
    std::cout << "Output plots saved to " << imagePath << std::endl;
    
    delete hRatioClone;
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
    
    hist1_norm->SetDirectory(0);
    hist2_norm->SetDirectory(0);

    TH1D* toNorm[] = {hist1_norm, hist2_norm};
    std::vector<TH1D> normed = normalizer(toNorm, 2, q1, q2, 1.0);
    
    TH1D* hRatio = (TH1D*)(normed[0].Clone(outputHistName));

    hRatio->Divide(&normed[1]);

    delete hist1_norm;
    delete hist2_norm;

    return hRatio;
}

TH1D* histhistRatioWithComp(TH1D* hist1, TH1D* hist2,
                            Double_t q1, Double_t q2,
                            const char* outputHistName,
                            const char* saveImagePrefix = nullptr)
{   
    if (!hist1 || !hist2) {
        std::cerr << "histhistRatioWithComp ERROR: null input histogram" << std::endl;
        return nullptr;
    }

    if (hist1->GetNbinsX() != hist2->GetNbinsX()) {
        std::cerr << "histhistRatioWithComp ERROR: different binning" << std::endl;
        return nullptr;
    }

    // 1. Prepare Clones for Normalization
    TH1D* hist1_norm = (TH1D*)(hist1->Clone("hist1_normalized"));
    TH1D* hist2_norm = (TH1D*)(hist2->Clone("hist2_normalized"));
    
    hist1_norm->SetDirectory(0);
    hist2_norm->SetDirectory(0);

    // 2. Perform Normalization
    TH1D* toNorm[] = {hist1_norm, hist2_norm};
    std::vector<TH1D> normed = normalizer(toNorm, 2, q1, q2, 1.0);
    
    // 3. Create and Save Comparison Plot (Before vs After)
    if (saveImagePrefix != nullptr) {
        TCanvas* cComp = new TCanvas("cComp", "Normalization Comparison", 1200, 600);
        cComp->Divide(2, 1);

        // --- Left Pad: BEFORE (Raw) ---
        cComp->cd(1);
        gStyle->SetOptStat(0);     
        gPad->SetGrid();
        gPad->SetLeftMargin(0.15);
        
        hist1->SetLineColor(kRed);
        hist1->SetTitle("Before Normalization (Raw)");
        hist1->Draw("HIST");
        
        hist2->SetLineColor(kBlue);
        hist2->Draw("HIST SAME");

        TLegend* leg1 = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg1->AddEntry(hist1, "Signal (Raw)", "l");
        leg1->AddEntry(hist2, "Background (Raw)", "l");
        leg1->Draw();

        // --- Right Pad: AFTER (Normalized) ---
        cComp->cd(2);
        gStyle->SetOptStat(0);     
        gPad->SetGrid();
        gPad->SetLeftMargin(0.15);

        // Get temporary pointers to the normalized results for drawing
        TH1D* hNormSig = (TH1D*)normed[0].Clone("draw_norm_sig");
        TH1D* hNormBkg = (TH1D*)normed[1].Clone("draw_norm_bkg");

        hNormSig->SetLineColor(kRed);
        hNormSig->SetTitle("After Normalization");
        hNormSig->Draw("HIST");

        hNormBkg->SetLineColor(kBlue);
        hNormBkg->Draw("HIST SAME");

        // Draw visual guides for the normalization range
        TLine* l1 = new TLine(q1, 0, q1, hNormSig->GetMaximum());
        TLine* l2 = new TLine(q2, 0, q2, hNormSig->GetMaximum());
        l1->SetLineStyle(2); l2->SetLineStyle(2);
        l1->Draw(); l2->Draw();

        TLegend* leg2 = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg2->AddEntry(hNormSig, "Signal (Norm)", "l");
        leg2->AddEntry(hNormBkg, "Background (Norm)", "l");
        leg2->Draw();

        // --- NEW: CHECK AND CREATE DIRECTORY ---
        const char* checkPath = "./imgs/test/norm_checks/"; 
        
        // gSystem->AccessPathName returns true if the path DOES NOT exist
        if (gSystem->AccessPathName(checkPath)) {
            // The 'true' argument makes it recursive (mkdir -p), creating parents if needed
            gSystem->mkdir(checkPath, true);
            std::cout << "Created directory: " << checkPath << std::endl;
        }

        TCanvas* canvArr[] = {cComp};
        save_canvas_images(canvArr, 1, checkPath, saveImagePrefix, "png");
        
        // Clean up visualization objects
        delete cComp; 
        delete hNormSig;
        delete hNormBkg;
    }

    // 4. Calculate Ratio
    TH1D* hRatio = (TH1D*)(normed[0].Clone(outputHistName));
    hRatio->Divide(&normed[1]);

    // 5. Cleanup
    delete hist1_norm;
    delete hist2_norm;

    return hRatio;
}

TH2D* histhistRatio2d(TH2D* hist1, TH2D* hist2,
                    Double_t q1x, Double_t q2x,
                    Double_t q1y, Double_t q2y,
                    const char* outputHistName) {
    if (!hist1 || !hist2) {
        std::cerr << "histhistRatio ERROR: null input histogram" << std::endl;
        return nullptr;
    }

    if (hist1->GetNbinsX() != hist2->GetNbinsX()) {
        std::cerr << "histhistRatio ERROR: different binning" << std::endl;
        return nullptr;
    }

    TH2D* hist1_norm = (TH2D*)(hist1->Clone("hist1_normalized"));
    TH2D* hist2_norm = (TH2D*)(hist2->Clone("hist2_normalized"));

    TH2D* toNorm[] = {hist1_norm, hist2_norm};
    normalizer(toNorm, 2, q1x, q2x, q1y, q2y, 1.0);

    TH2D* hRatio = (TH2D*)(hist1_norm->Clone(outputHistName));

    hRatio->Divide(hist2_norm);

    delete hist1_norm;
    delete hist2_norm;

    return hRatio;
}

TH2D* histhistRatioDeltaEtaDeltaPhi(TH2D* hist1, TH2D* hist2,
                    Double_t etaCut,
                    const char* outputHistName) {
    if (!hist1 || !hist2) {
        std::cerr << "histhistRatio ERROR: null input histogram" << std::endl;
        return nullptr;
    }

    if (hist1->GetNbinsX() != hist2->GetNbinsX()) {
        std::cerr << "histhistRatio ERROR: different binning" << std::endl;
        return nullptr;
    }

    TH2D* hist1_norm = (TH2D*)(hist1->Clone("hist1_normalized"));
    TH2D* hist2_norm = (TH2D*)(hist2->Clone("hist2_normalized"));

    double normA = hist1_norm->Integral(hist1_norm->GetXaxis()->FindBin(0.04), hist1_norm->GetNbinsX(), 0., hist1_norm->GetNbinsY());
    double normB = hist2_norm->Integral(hist2_norm->GetXaxis()->FindBin(0.04), hist2_norm->GetNbinsX(), 0., hist2_norm->GetNbinsY());

    TH2D* hRatio = (TH2D*)(hist1_norm->Clone(outputHistName));

    hRatio->Divide(hist2_norm);
    hRatio->Scale(normB / normA);

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