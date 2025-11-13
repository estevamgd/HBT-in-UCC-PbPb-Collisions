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


void getSingleRatio(const char* ss_FileName, const char* ss_HistName,
                    const char* os_FileName, const char* os_HistName,
                    const char* outputPrefix, const char* outputHistName = "sr_cor",
                    const char* plotHeaderText = "PbPb 2.76 TeV | Single Ratio") {
    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    TFile *ssFile = TFile::Open(ss_FileName, "READ");
    if (!ssFile || ssFile->IsZombie()) {
        std::cerr << "Error: Could not open Same-Sign file '" << ss_FileName << "'!" << std::endl;
        return;
    }
    TH1D *hSS = (TH1D*)ssFile->Get(ss_HistName);
    if (!hSS) {
        std::cerr << "Error: Same-Sign histogram '" << ss_HistName << "' not found in '" << ss_FileName << "'!" << std::endl;
        ssFile->Close();
        delete ssFile;
        return;
    }

    TFile *osFile = TFile::Open(os_FileName, "READ");
    if (!osFile || osFile->IsZombie()) {
        std::cerr << "Error: Could not open Opposite-Sign file '" << os_FileName << "'!" << std::endl;
        ssFile->Close(); delete ssFile;
        return;
    }
    TH1D *hOS = (TH1D*)osFile->Get(os_HistName);
    if (!hOS) {
        std::cerr << "Error: Opposite-Sign histogram '" << os_HistName << "' not found in '" << os_FileName << "'!" << std::endl;
        osFile->Close(); delete osFile;
        ssFile->Close(); delete ssFile;
        return;
    }

    TH1D *hSS_norm = (TH1D*)hSS->Clone("hSS_normalized");
    hSS_norm->SetTitle("Same-Sign Pairs (Normalized)");

    TH1D *hOS_norm = (TH1D*)hOS->Clone("hOS_normalized");
    hOS_norm->SetTitle("Opposite-Sign Pairs (Normalized)");

    std::cout << "Successfully loaded Same-Sign and Opposite-Sign histograms." << std::endl;

    TH1D *tonormhist[] = {hSS_norm, hOS_norm};
    int numTonorm = 2;
    Double_t scale = 1;
    normalizer(tonormhist, numTonorm, scale);

    TH1D *hRatio = (TH1D*)hSS_norm->Clone(outputHistName);
    hRatio->Divide(hOS_norm);

    hRatio->SetTitle("; q_{inv} (GeV); C(q_{inv}) = SS/OS");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.8);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);

    std::cout << "Single Ratio histogram '" << outputHistName << "' created successfully." << std::endl;

    TCanvas *cRatio = new TCanvas("cRatio", "Single Ratio", 1200, 800);
    cRatio->cd();
    cRatio->SetLeftMargin(0.12);
    cRatio->SetBottomMargin(0.12);
    hRatio->GetYaxis()->SetTitleOffset(1.2);

    hRatio->Draw("E1 PLC");
    
    TLegend *legendRatio = new TLegend(0.6, 0.8, 0.88, 0.88);
    legendRatio->AddEntry(hRatio, "Single Ratio (SS/OS)", "lep");
    legendRatio->SetFillStyle(0);
    legendRatio->SetBorderSize(0);
    legendRatio->Draw();

    drawCMSHeaders("#bf{CMS} #it{Preliminary}", plotHeaderText);

    const char* histoPath = "./data/correlation_ratios/";
    TH1D *histogramsToSave[] = {hRatio, hSS_norm, hOS_norm};

    save_histograms(histogramsToSave, 3, histoPath, outputPrefix);
    std::cout << "Output histograms saved to a new ROOT file in " << histoPath << std::endl;

    const char* imagePath = "./imgs/correlation_ratios/";
    TCanvas *canvasesToSave[] = { cRatio };
    save_canvas_images(canvasesToSave, 1, imagePath, outputPrefix, "png");
    save_canvas_images(canvasesToSave, 1, imagePath, outputPrefix, "pdf");
    std::cout << "Output plots saved to " << imagePath << std::endl;

    TH1D *histogramsToClose[] = { hSS_norm, hOS_norm, hRatio };
    TCanvas *canvasesToClose[] = { cRatio };
    TLegend *legendsToClose[] = { legendRatio };

    close_program(canvasesToClose, 1, histogramsToClose, 3, legendsToClose, 1, ssFile);
    
    osFile->Close();
    delete osFile;
    
    std::cout << "Cleanup complete." << std::endl;
}

int main() {
    // To compile use:
    // g++ -std=c++17 -pthread getSingleRatio.cpp -o getSingleRatio `root-config --cflags --libs`
    // To run use:
    // ./getSingleRatio

    // --- Configuration ---
    // Define the analysis parameters
    ControlVar selectedControlVar = ControlVar::CENTHF;
    double bin_low = 3200.0;
    double bin_high = 3300.0;
    
    const char* selectionVarName = getSelVarName(selectedControlVar);

    TString searchPattern = TString::Format(
        "./data/signal_mix/sig_double_loop*%s_%f-%f*.root", 
        selectionVarName, bin_low, bin_high
    );

    // Use findFile (from my_func.h) to get the newest file matching the pattern
    TString dataFile = findFile(searchPattern);

    if (dataFile.IsNull()) {
        std::cerr << "Error: No data file found matching pattern: " << searchPattern << std::endl;
        std::cerr << "Please run the 'sig_double_loop' or 'sig_double_loop_parallel' analysis first." << std::endl;
        return 1;
    }
    std::cout << "Found data file: " << dataFile << std::endl;

    // --- Analysis Set 1: Corrected Ratio ---
    std::cout << "\n--- Running Analysis 1: Corrected Ratio ---" << std::endl;
    
    const char* ssHist_cor = "h_qinvSSCor_signal_2l"; 
    const char* osHist_cor = "h_qinvOSCor_signal_2l"; 

    TString outputPrefix_cor = TString::Format(
        "SingleRatio_Cor_%s_%.0f-%.0f", 
        selectionVarName, bin_low, bin_high
    );
    const char* outputHist_cor = "sr_cor"; 
    
    TString plotHeader_cor = TString::Format(
        "PbPb 2.76 TeV | %s %.0f-%.0f (Cor)",
        selectionVarName, bin_low, bin_high
    );

    getSingleRatio(dataFile.Data(), ssHist_cor,
                   dataFile.Data(), osHist_cor, 
                   outputPrefix_cor.Data(), outputHist_cor,
                   plotHeader_cor.Data()); 

    // --- Analysis Set 2: Uncorrected Ratio  ---
    std::cout << "\n--- Running Analysis 2: Uncorrected Ratio ---" << std::endl;

    const char* ssHist_uncor = "h_qinvSS_signal_2l"; 
    const char* osHist_uncor = "h_qinvOS_signal_2l"; 

    TString outputPrefix_uncor = TString::Format(
        "SingleRatio_Uncor_%s_%.0f-%.0f", 
        selectionVarName, bin_low, bin_high
    );

    const char* outputHist_uncor = "sr_uncor";
    
    TString plotHeader_uncor = TString::Format(
        "PbPb 2.76 TeV | %s %.0f-%.0f (Uncor)",
        selectionVarName, bin_low, bin_high
    );

    getSingleRatio(dataFile.Data(), ssHist_uncor,
                   dataFile.Data(), osHist_uncor,
                   outputPrefix_uncor.Data(), outputHist_uncor,
                   plotHeader_uncor.Data());

    std::cout << "\nAll Single Ratio creation processes finished." << std::endl;
    return 0;
}