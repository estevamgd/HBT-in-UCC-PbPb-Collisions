#ifndef FITRANGESCAN_H
#define FITRANGESCAN_H

#include <cmath>
#include <iostream>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "../include/fits.h"
#include "../include/my_func.h"

struct FitRangeScanConfig {
    FitFunctionType fitType;

    double fitMinLow;
    double fitMinHigh;
    double fitMinStep;

    double fitMaxLow;
    double fitMaxHigh;
    double fitMaxStep;

    const char* outDir;
};


inline void fitRangeScanHistogram(
    TH1D* hist,
    const FitRangeScanConfig& cfg,
    const FitInit& fitInit,
    const char* subtitle
) {
    if (!hist) {
        std::cerr << "ERROR: Null histogram passed to fitRangeScanHistogram\n";
        return;
    }

    const FitModelConfig fitCfg = getFitModelConfig(cfg.fitType);

    int nBinsMin =
        std::ceil((cfg.fitMinHigh - cfg.fitMinLow) / cfg.fitMinStep);

    int nBinsMax =
        std::ceil((cfg.fitMaxHigh - cfg.fitMaxLow) / cfg.fitMaxStep);

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
                    hist,
                    &res,
                    cfg.fitType,
                    fitInit,
                    fitMin,
                    fitMax
                );
            } catch (...) {
                continue;
            }

            if (res.Get() && res->IsValid() && res->Ndf() > 0) {
                hCL->Fill(fitMin, fitMax, res->Prob());
                hChi2NDF->Fill(
                    fitMin,
                    fitMax,
                    res->Chi2() / res->Ndf()
                );
            }

            delete f;
        }
    }

    // ===== Save ROOT =====
    TFile outFile(
        Form("%s/fitRangeScan_%s.root",
             cfg.outDir,
             fitCfg.name),
        "RECREATE"
    );
    hCL->Write();
    hChi2NDF->Write();
    outFile.Close();

    // ===== Plotting =====
    gStyle->SetOptStat(0);

    TCanvas* cCL =
        new TCanvas("cFitRangeCL", "Fit range scan (p-value)", 1100, 900);
    cCL->SetRightMargin(0.15);
    hCL->Draw("COLZ");

    drawCMSHeaders(
        "#bf{CMS} #it{Work in Progress}",
        Form("%s | %s", fitCfg.displayName, subtitle)
    );

    TCanvas* cChi2 =
        new TCanvas("cFitRangeChi2NDF", "Fit range scan (#chi^{2}/NDF)", 1100, 900);
    cChi2->SetRightMargin(0.15);
    hChi2NDF->Draw("COLZ");

    drawCMSHeaders(
        "#bf{CMS} #it{Work in Progress}",
        Form("%s | %s", fitCfg.displayName, subtitle)
    );

    TCanvas* canvCL[]   = { cCL };
    TCanvas* canvChi2[] = { cChi2 };

    save_canvas_images(
        canvCL, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_CL", fitCfg.name),
        "png"
    );
    save_canvas_images(
        canvCL, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_CL", fitCfg.name),
        "pdf"
    );

    save_canvas_images(
        canvChi2, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_Chi2NDF", fitCfg.name),
        "png"
    );
    save_canvas_images(
        canvChi2, 1,
        "./imgs/test/fit_range_scan/",
        Form("%s_Chi2NDF", fitCfg.name),
        "pdf"
    );

    delete cCL;
    delete cChi2;
    delete hCL;
    delete hChi2NDF;
}

#endif
