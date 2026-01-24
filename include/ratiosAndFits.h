#ifndef RATIOSANDFITS_H
#define RATIOSANDFITS_H

#include "../include/ratios.h"
#include "../include/data_func.h"
#include "../include/my_func.h"
#include "../include/fits.h"
#include "../include/fitRangeScan.h"

void doubleRatioFit(ControlVar selectedControlVar,
                    std::vector<std::pair<FitFunctionType, FitInit>> models,
                    double bin_low, double bin_high,
                    Double_t q1, Double_t q2,
                    qMode mode,
                    double fitMin = 0.0, double fitMax = 10.0, double fitMinBG = 0.15,
                    double plotXMin = 0.0, double plotXMax = 10.0, 
                    double plotYMin = 0.0, double plotYMax = 2.1,
                    const FitRangeScanConfig* scanCfgSR = nullptr,
                    const FitRangeScanConfig* scanCfgDR = nullptr,
                    const char* fileName = nullptr)
{   
    const char* qmodename = (mode == qMode::QINV) ? "qinv" : "qlcms";
    const char* selectionVarName = getSelVarName(selectedControlVar);
    TString searchPattern;

    if (!fileName) {
        searchPattern = TString::Format(
            "./data/signal_mix/sig_%s_double_loop*%s_%f-%f*.root",
            qmodename, selectionVarName, bin_low, bin_high
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/sig_%s_double_loop*%s_%f-%f_%s.root",
            qmodename, selectionVarName, bin_low, bin_high, fileName
        );
    }
    const char* name_ssHist = Form("h_%sSSCor_signal_2l", qmodename);
    const char* name_osHist = Form("h_%sOSCor_signal_2l", qmodename);
    
    TH1D* q_ssHist_cor = getHistogram(searchPattern.Data(), name_ssHist);
    TH1D* q_osHist_cor = getHistogram(searchPattern.Data(), name_osHist);
    
    // ===== Single ratio ===== //
    TH1D* singleRatio = histhistRatio(
        q_ssHist_cor, q_osHist_cor, q1, q2, "sr_cor"
    );

    saveRatio(singleRatio, 
        Form("sr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high), 
        plotXMin, plotXMax, plotYMin, plotYMax,
        "sr_cor", "PbPb 2.76 TeV | Single Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C(q_{LCMS}) = SS/OS"
            : "; q_{inv} [GeV]; C(q_{inv}) = SS/OS",
        "Single Ratio (SS/OS)"
    );

    // ===== Single Ratio Fits ===== //
    std::vector<std::pair<FitFunctionType, FitInit>> models1 = models;
    
    auto fits1 = fitHistogramMultiple(
        singleRatio,
        models1,
        fitMin, fitMax
    );

    // ===== Scan on Single Ratio =====
    const char* subtitleSR =  "Single Ratio";
    if (scanCfgSR) {
        for (const auto& m : models1) {
            if (m.first == scanCfgSR->fitType) {
                fitRangeScanHistogram(
                    singleRatio,
                    *scanCfgSR,
                    m.second,
                    subtitleSR
                );
                break;
            }
        }
    }

    // ===== Draw & Save once ===== //
    saveFits(
        singleRatio,
        fits1,
        "c_single_ratio_fits",
        Form("PbPb 2.76 TeV | Single Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C(q_{LCMS}) = SS/OS"
            : "; q_{inv} [GeV]; C(q_{inv}) = SS/OS",
        mode, plotXMin, plotXMax, plotYMin, plotYMax);
    
    // ===== Background Fit ===== //
    FitInit bgInit{{1.0, 0.1, 2.0, 0.1, 2.0}};
    TFitResultPtr bgRes;

    TF1* bgFit = fitHistogram(
        singleRatio,
        &bgRes,
        FitFunctionType::BACKGROUND,
        bgInit,
        fitMinBG
    );

    // ===== Double Ratio ===== //
    TH1D* doubleRatio = histfuncRatio(singleRatio, bgFit, "dr_cor");

    saveRatio(doubleRatio,
        Form("dr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax,
        "dr_cor", "PbPb 2.76 TeV | Double Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C(q_{inv}) = Data/Fit",
        "Double Ratio (Data/Fit)"
    );

    // ===== Final Fits ===== //
    std::vector<std::pair<FitFunctionType, FitInit>> models2 = models;
    
    auto fits = fitHistogramMultiple(
        doubleRatio,
        models2,
        fitMin, fitMax
    );
    
    // ===== Draw & Save once ===== //
    saveFits(
        doubleRatio,
        fits,
        "c_double_ratio_fits",
        Form("PbPb 2.76 TeV | Double Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C(q_{inv}) = Data/Fit",
        mode, plotXMin, plotXMax, plotYMin, plotYMax);
    
    // ===== Scan on Double Ratio =====
    const char* subtitleDR =  "Double Ratio";
    if (scanCfgDR) {
        for (const auto& m : models2) {
            if (m.first == scanCfgDR->fitType) {
                fitRangeScanHistogram(
                    doubleRatio,
                    *scanCfgDR,
                    m.second,
                    subtitleDR
                );
                break;
            }
        }
    }
    

    // Cleanup memory
    for (auto& f : fits)
    delete f.function;
    delete bgFit;
    delete q_ssHist_cor;
    delete q_osHist_cor;
    delete singleRatio;
    delete doubleRatio;
    
}

#endif