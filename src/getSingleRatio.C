#include <iostream>
#include <string>
#include <filesystem>

#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../include/my_func.h" 

/**
 * @brief Creates the single ratio (SR) correlation histogram by dividing Same-Sign pairs by Opposite-Sign pairs.
 *
 * This function implements the standard single ratio method. It reads a Same-Sign (SS) histogram
 * and an Opposite-Sign (OS) histogram from ROOT files, normalizes them to have an integral of 1,
 * and then computes the ratio SR = SS / OS. The resulting histogram and its normalized components
 * are saved into a new ROOT file, ready for fitting.
 *
 * @param ss_FileName     Path to the ROOT file containing the Same-Sign (SS) histogram.
 * @param ss_HistName     Name of the SS histogram in the file.
 * @param os_FileName     Path to the ROOT file for the Opposite-Sign (OS) histogram.
 * @param os_HistName     Name of the OS histogram in the file.
 * @param outputPrefix    Base name for the output ROOT file.
 * @param outputHistName  Name for the final SR histogram inside the new ROOT file.
 */
void getSingleRatio(const char* ss_FileName, const char* ss_HistName,
                         const char* os_FileName, const char* os_HistName,
                         const char* outputPrefix, const char* outputHistName = "sr_cor") {

    ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);

    // 1. Open the Same-Sign file and get the histogram
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
    // Clone to work with a local copy
    TH1D *hSS_norm = (TH1D*)hSS->Clone("hSS_normalized");
    hSS_norm->SetTitle("Same-Sign Pairs (Normalized)");


    // 2. Open the Opposite-Sign file and get the histogram
    TFile *osFile = TFile::Open(os_FileName, "READ");
    if (!osFile || osFile->IsZombie()) {
        std::cerr << "Error: Could not open Opposite-Sign file '" << os_FileName << "'!" << std::endl;
        ssFile->Close(); delete ssFile;
        delete hSS_norm;
        return;
    }
    TH1D *hOS = (TH1D*)osFile->Get(os_HistName);
    if (!hOS) {
        std::cerr << "Error: Opposite-Sign histogram '" << os_HistName << "' not found in '" << os_FileName << "'!" << std::endl;
        osFile->Close(); delete osFile;
        ssFile->Close(); delete ssFile;
        delete hSS_norm;
        return;
    }
    // Clone to work with a local copy
    TH1D *hOS_norm = (TH1D*)hOS->Clone("hOS_normalized");
    hOS_norm->SetTitle("Opposite-Sign Pairs (Normalized)");


    std::cout << "Successfully loaded Same-Sign and Opposite-Sign histograms." << std::endl;

    // 3. Normalize both histograms to have an integral of 1.0
    if (hSS_norm->Integral() > 0) {
        hSS_norm->Scale(1.0 / hSS_norm->Integral());
    }
    if (hOS_norm->Integral() > 0) {
        hOS_norm->Scale(1.0 / hOS_norm->Integral());
    }

    // 4. Create the ratio histogram by dividing SS / OS
    TH1D *hRatio = (TH1D*)hSS_norm->Clone(outputHistName);
    hRatio->Divide(hOS_norm);

    // 5. Style the final ratio histogram
    hRatio->SetTitle("Single Ratio (SS/OS); q_{inv} (GeV/c); C(q_{inv}) = SS/OS");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);

    std::cout << "Single Ratio histogram '" << outputHistName << "' created successfully." << std::endl;

    // 6. Save the results into a new ROOT file using your helper function
    const char* histoPath = "./data/correlation_ratios/";
    TH1D *histogramsToSave[] = {
        hRatio,     // The final result: SS / OS
        hSS_norm,   // The normalized Same-Sign histogram
        hOS_norm    // The normalized Opposite-Sign histogram
    };
    save_histograms(histogramsToSave, 3, histoPath, outputPrefix);

    std::cout << "Output saved to a new ROOT file in " << histoPath << std::endl;

    // 7. Clean up memory
    delete hRatio;
    delete hSS_norm;
    delete hOS_norm;
    ssFile->Close();
    osFile->Close();
    delete ssFile;
    delete osFile;
}

void execute_division() {
    // --- Configuration ---

    // 1. Define the Same-Sign (SS) file and histogram name. This is the numerator.
    const char* ssFile = "./data/sigal_mix/signal_mix_CENTHF_3200.000000-3300.000000.root";
    const char* ssHist = "h_qinvSSCor_signal_1l"; // Corrected Same-Sign pairs

    // 2. Define the Opposite-Sign (OS) file and histogram name. This is the denominator.
    //    This can be the same file as the SS file if both histograms are in one place.
    const char* osFile = "./data/sigal_mix/signal_mix_CENTHF_3200.000000-3300.000000.root";
    const char* osHist = "h_qinvOSCor_signal_1l"; // Corrected Opposite-Sign pairs

    // 3. Define the prefix for the output file and the name for the new ratio histogram.
    const char* outputPrefix = "SingleRatio_CENTHF_3200-3300";
    const char* outputHist = "sr_cor"; // "sr_cor" stands for "single ratio corrected"

    // --- Run the Macro ---
    getSingleRatio(ssFile, ssHist,
                        osFile, osHist,
                        outputPrefix, outputHist);

    std::cout << "\nSingle Ratio creation process finished." << std::endl;
}