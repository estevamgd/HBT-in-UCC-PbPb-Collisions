#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TFitResultPtr.h"

#include "../include/ratios.h"
#include "../include/fits.h"
#include "../include/data_func.h"
#include "../include/my_func.h"


struct FitRangeScanConfig {
    ControlVar selectedControlVar;
    qMode mode;

    double bin_low;
    double bin_high;

    double norm_q1;
    double norm_q2;

    FitFunctionType fitType;
    FitInit fitInit;

    double fitMinLow;
    double fitMinHigh;
    double fitMinStep;

    double fitMaxLow;
    double fitMaxHigh;
    double fitMaxStep;

    const char* outDir;
    const char* tag;
};


void findBestFitRangeScan(const FitRangeScanConfig& cfg) {

    // ================================
    // Basic labels
    // ================================
    const char* qmodename =
        (cfg.mode == qMode::QINV) ? "qinv" : "qlcms";

    const char* selVarName =
        getSelVarName(cfg.selectedControlVar);

    // ================================
    // Load histograms
    // ================================
    TString pattern = TString::Format(
        "./data/signal_mix/sig_%s_double_loop*%s_%f-%f*.root",
        qmodename, selVarName, cfg.bin_low, cfg.bin_high
    );

    TH1D* hSS = getHistogram(
        pattern,
        Form("h_%sSSCor_signal_2l", qmodename)
    );

    TH1D* hOS = getHistogram(
        pattern,
        Form("h_%sOSCor_signal_2l", qmodename)
    );

    if (!hSS || !hOS) {
        std::cerr << "ERROR: Failed to load input histograms\n";
        return;
    }

    TH1D* singleRatio =
        histhistRatio(hSS, hOS,
                      cfg.norm_q1,
                      cfg.norm_q2,
                      "sr_fitRangeScan");

    // ================================
    // Scan binning
    // ================================
    int nBinsMin =
        std::ceil((cfg.fitMinHigh - cfg.fitMinLow) / cfg.fitMinStep);

    int nBinsMax =
        std::ceil((cfg.fitMaxHigh - cfg.fitMaxLow) / cfg.fitMaxStep);

    // ================================
    // Output histograms
    // ================================
    TH2D* hCL = new TH2D(
        "hFitRangeCL",
        ";Fit min [GeV];Fit max [GeV];p-value",
        nBinsMin, cfg.fitMinLow, cfg.fitMinHigh,
        nBinsMax, cfg.fitMaxLow, cfg.fitMaxHigh
    );

    TH2D* hChi2NDF = new TH2D(
        "hFitRangeChi2NDF",
        ";Fit min [GeV];Fit max [GeV];#chi^{2}/NDF",
        nBinsMin, cfg.fitMinLow, cfg.fitMinHigh,
        nBinsMax, cfg.fitMaxLow, cfg.fitMaxHigh
    );

    // ================================
    // Scan fits
    // ================================
    for (double fitMin = cfg.fitMinLow;
         fitMin <= cfg.fitMinHigh;
         fitMin += cfg.fitMinStep) {

        for (double fitMax = cfg.fitMaxLow;
             fitMax <= cfg.fitMaxHigh;
             fitMax += cfg.fitMaxStep) {

            if (fitMax <= fitMin) continue;

            TFitResultPtr res;
            TF1* f = nullptr;

            try {
                f = fitHistogram(
                    singleRatio,
                    &res,
                    cfg.fitType,
                    cfg.fitInit,
                    fitMin,
                    fitMax
                );
            } catch (...) {
                continue;
            }

            if (res.Get() && res->IsValid()) {
                double chi2ndf =
                    (res->Ndf() > 0)
                        ? res->Chi2() / res->Ndf()
                        : 0.0;

                hCL->Fill(fitMin, fitMax, res->Prob());
                hChi2NDF->Fill(fitMin, fitMax, chi2ndf);
            }

            delete f;
        }
    }

    // ================================
    // Save ROOT output
    // ================================
    TFile outFile(
        Form("%s/fitRangeScan_%s_%s_%.0f-%.0f.root",
             cfg.outDir,
             qmodename,
             selVarName,
             cfg.bin_low,
             cfg.bin_high),
        "RECREATE"
    );

    hCL->Write();
    hChi2NDF->Write();

    outFile.Close();

    // ================================
    // Plotting
    // ================================
    gStyle->SetOptStat(0);

    const FitModelConfig fitCfg =
        getFitModelConfig(cfg.fitType);

    // ---- p-value canvas
    TCanvas* cCL =
        new TCanvas("cFitRangeCL",
                     "Fit range scan (p-value)",
                     1100, 900);

    cCL->SetRightMargin(0.15);
    hCL->Draw("COLZ");

    drawCMSHeaders(
        "#bf{CMS} #it{Work in Progress}",
        Form("%s | %s %.0f-%.0f",
             fitCfg.displayName,
             selVarName,
             cfg.bin_low,
             cfg.bin_high)
    );

    // ---- chi2/NDF canvas
    TCanvas* cChi2 =
        new TCanvas("cFitRangeChi2NDF",
                     "Fit range scan (chi2/NDF)",
                     1100, 900);

    cChi2->SetRightMargin(0.15);
    hChi2NDF->Draw("COLZ");

    drawCMSHeaders(
        "#bf{CMS} #it{Work in Progress}",
        Form("%s | %s %.0f-%.0f",
             fitCfg.displayName,
             selVarName,
             cfg.bin_low,
             cfg.bin_high)
    );

    // ================================
    // Save images
    // ================================
    TCanvas* canvasesCL[]   = { cCL };
    TCanvas* canvasesChi2[] = { cChi2 };

    save_canvas_images(
        canvasesCL, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_pvalue", cfg.tag),
        "png"
    );

    save_canvas_images(
        canvasesCL, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_pvalue", cfg.tag),
        "pdf"
    );

    save_canvas_images(
        canvasesChi2, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_chi2ndf", cfg.tag),
        "png"
    );

    save_canvas_images(
        canvasesChi2, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_chi2ndf", cfg.tag),
        "pdf"
);

    // ================================
    // Cleanup
    // ================================
    delete cCL;
    delete cChi2;
    delete hCL;
    delete hChi2NDF;
    delete singleRatio;
    delete hSS;
    delete hOS;
}


int main() {
    // ============================================================
    // Compilation:
    // cd src
    // g++ -std=c++17 findBestFitRange.cpp -o findBestFitRange `root-config --cflags --libs`
    //
    // Execution:
    // ./findBestFitRange
    // ============================================================

    FitRangeScanConfig cfg;

    cfg.selectedControlVar = ControlVar::CENTHF;
    cfg.mode = qMode::QINV;

    cfg.bin_low  = 3200.0;
    cfg.bin_high = 3300.0;

    cfg.norm_q1 = 4.82;
    cfg.norm_q2 = 6.40;

    cfg.fitType = FitFunctionType::LEVY;
    cfg.fitInit = {{1.0, 0.5, 5.0, 0.0, 1.2}};

    cfg.fitMinLow  = 0.05;
    cfg.fitMinHigh = 0.10;
    cfg.fitMinStep = 0.005;

    cfg.fitMaxLow  = 0.10;
    cfg.fitMaxHigh = 0.30;
    cfg.fitMaxStep = 0.005;

    cfg.outDir = "./data/fit_range_scan/";
    cfg.tag    = "CL_scan";

    findBestFitRangeScan(cfg);

    return 0;
}
