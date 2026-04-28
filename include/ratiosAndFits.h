#ifndef RATIOSANDFITS_H
#define RATIOSANDFITS_H

#include "../include/ratios.h"
#include "../include/data_func.h"
#include "../include/my_func.h"
#include "../include/fits.h"
#include "../include/fitRangeScan.h"
#include <stdexcept>


struct MultipleFitResults { 
    std::vector<TFitResultPtr> results;
    std::vector<std::string> modelNames;
};

struct AlphaSeedFitOutput {
    FitOutput fixedAlphaFit;
    FitOutput freeAlphaFit;
    FitInit freeAlphaInit;
};

FitInit findModelInitOrDefault(
    const std::vector<std::pair<FitFunctionType, FitInit>>& models,
    FitFunctionType type,
    const FitInit& fallback)
{
    for (const auto& model : models) {
        if (model.first == type) {
            return model.second;
        }
    }

    return fallback;
}

AlphaSeedFitOutput fitDoubleLevyAlphaOneThenFreeDetailed(
    TH1D* hist,
    const FitInit& init,
    double fitMin = 0.0,
    double fitMax = 10.0)
{
    auto cfg = getFitModelConfig(FitFunctionType::DOUBLE_LEVY);
    if ((int)init.values.size() != cfg.nPar) {
        throw std::runtime_error("Wrong number of initial parameters for DOUBLE_LEVY");
    }

    static int alphaOneSeedCounter = 0;
    TString fixedFitName = TString::Format("fitLevyDR_alpha1_fixed_%d", alphaOneSeedCounter++);
    TF1* fixedAlphaFit = new TF1(fixedFitName, cfg.func, fitMin, fitMax, cfg.nPar);

    for (int i = 0; i < cfg.nPar; ++i) {
        fixedAlphaFit->SetParameter(i, init.values[i]);
        fixedAlphaFit->SetParName(i, cfg.parNames[i].c_str());

        if (i == 2) {
            continue;
        }

        if (i < (int)cfg.parLimits.size() &&
            cfg.parLimits[i].first != cfg.parLimits[i].second) {
            fixedAlphaFit->SetParLimits(
                i,
                cfg.parLimits[i].first,
                cfg.parLimits[i].second
            );
        }
    }

    fixedAlphaFit->FixParameter(2, 1.0);
    TFitResultPtr fixedAlphaResult = hist->Fit(fixedAlphaFit, "S R E M");

    std::vector<double> seededValues = init.values;
    for (int i = 0; i < cfg.nPar; ++i) {
        seededValues[i] = fixedAlphaFit->GetParameter(i);
    }
    seededValues[2] = 1.0;

    TFitResultPtr freeAlphaResult;
    TF1* finalFit = fitHistogram(
        hist,
        &freeAlphaResult,
        FitFunctionType::DOUBLE_LEVY,
        FitInit{seededValues},
        fitMin,
        fitMax
    );

    return {
        {
            FitFunctionType::DOUBLE_LEVY,
            fixedAlphaFit,
            fixedAlphaResult,
            cfg.legendParams,
            "Levy Fit (#alpha fixed to 1)"
        },
        {
            FitFunctionType::DOUBLE_LEVY,
            finalFit,
            freeAlphaResult,
            cfg.legendParams,
            "Levy Fit (#alpha free, #alpha=1 seed)"
        },
        FitInit{seededValues}
    };
}

FitOutput fitDoubleLevyAlphaOneThenFree(
    TH1D* hist,
    const FitInit& init,
    double fitMin = 0.0,
    double fitMax = 10.0)
{
    AlphaSeedFitOutput detailedFit = fitDoubleLevyAlphaOneThenFreeDetailed(
        hist,
        init,
        fitMin,
        fitMax
    );

    delete detailedFit.fixedAlphaFit.function;

    return {
        detailedFit.freeAlphaFit.type,
        detailedFit.freeAlphaFit.function,
        detailedFit.freeAlphaFit.result,
        detailedFit.freeAlphaFit.params,
        detailedFit.freeAlphaFit.displayName
    };
}

void doubleRatioFit(ControlVar selectedControlVar,
                    std::vector<std::pair<FitFunctionType, FitInit>> models,
                    double bin_low, double bin_high,
                    Double_t etaMin, Double_t etaMax,
                    Double_t ptMin,
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
            "./data/signal_mix/"
            "%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d*.root",
            qmodename,etaMin, (int)bin_low, (int)bin_high, (int)bin_low, (int)bin_high 
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/"
            "%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d-%s.root",
            qmodename,etaMin,(int)bin_low, (int)bin_high, (int)bin_low, (int)bin_high, fileName 
        );
    }

    const char* name_ssHist = Form("h_%sSSCor_signal_2l", qmodename);
    const char* name_osHist = Form("h_%sOSCor_signal_2l", qmodename);
    
    TH1D* q_ssHist_cor = getHistogram(searchPattern.Data(), name_ssHist);
    TH1D* q_osHist_cor = getHistogram(searchPattern.Data(), name_osHist);
    
    // ===== Single ratio ===== //
    TH1D* singleRatio = histhistRatioWithComp(
        q_ssHist_cor, q_osHist_cor, q1, q2, "sr_cor", "testOS"
    );

    saveRatio(singleRatio, 
        Form("sr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high), 
        plotXMin, plotXMax, plotYMin, plotYMax,
        "sr_cor", "PbPb 2.76 TeV | Single Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = SS/OS"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = SS/OS",
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
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = SS/OS"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = SS/OS",
        mode, plotXMin, plotXMax, plotYMin, plotYMax);
    
    // ===== Background Fit ===== //
    FitInit bgInit{{1.0, 0.1, 2.0, 0.1, 2.0}};
    TFitResultPtr bgRes;
    std::vector<std::pair<FitFunctionType, FitInit>> bgModel = {
        {FitFunctionType::BACKGROUND, bgInit}
    };
    
    auto fitBg = fitHistogramMultiple(
        singleRatio,
        bgModel,
        fitMinBG
    );

    TF1* bgFit = fitHistogram(
        singleRatio,
        &bgRes,
        FitFunctionType::BACKGROUND,
        bgInit,
        fitMinBG
    );

    // ===== Draw & Save once ===== //
    saveFits(
        singleRatio,
        fitBg,
        "c_single_ratio_bg_fit",
        Form("PbPb 2.76 TeV | Single Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = SS/OS"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = SS/OS",
        mode, plotXMin, plotXMax, plotYMin, plotYMax);

    // ===== Double Ratio ===== //
    TH1D* doubleRatio = histfuncRatio(singleRatio, bgFit, "dr_cor");

    saveRatio(doubleRatio,
        Form("dr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax,
        "dr_cor", "PbPb 2.76 TeV | Double Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Data/Fit",
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
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Data/Fit",
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

MultipleFitResults doubleRatioMixFit(ControlVar selectedControlVar,
                    std::vector<std::pair<FitFunctionType, FitInit>> models,
                    double bin_low, double bin_high,
                    Double_t etaMin, Double_t etaMax,
                    Double_t ptMin,
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
            "./data/signal_mix/"
            "%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d*.root",
            qmodename,etaMin, (int)bin_low, (int)bin_high, (int)bin_low, (int)bin_high 
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/"
            "%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d-%s.root",
            qmodename,etaMin,(int)bin_low, (int)bin_high, (int)bin_low, (int)bin_high, fileName 
        );
    }

    const char* name_sigHist = "hSigSS";
    const char* name_mixHist = "hMixSS";
    
    TH1D* q_sigHist_cor = getHistogram(searchPattern.Data(), name_sigHist);
    TH1D* q_mixHist_cor = getHistogram(searchPattern.Data(), name_mixHist);
    
    // ===== Single ratio ===== //
    TH1D* singleRatio = histhistRatioWithComp(
        q_sigHist_cor, q_mixHist_cor, q1, q2, "sr_cor", "testMix"
    );

    saveRatio(singleRatio, 
        Form("sr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high), 
        plotXMin, plotXMax, plotYMin, plotYMax,
        "sr_cor", "PbPb 2.76 TeV | Single Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Sig/Mix",
        "Single Ratio (Sig/Mix)"
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
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Sig/Mix",
        mode, plotXMin, plotXMax, plotYMin, plotYMax);
    
    // ===== Background Fit ===== //
    FitInit bgInit{{1.0, 0.1, 2.0, 0.1, 2.0}};
    TFitResultPtr bgRes;
    std::vector<std::pair<FitFunctionType, FitInit>> bgModel = {
        {FitFunctionType::BACKGROUND, bgInit}
    };
    
    auto fitBg = fitHistogramMultiple(
        singleRatio,
        bgModel,
        fitMinBG
    );

    TF1* bgFit = fitHistogram(
        singleRatio,
        &bgRes,
        FitFunctionType::BACKGROUND,
        bgInit,
        fitMinBG
    );

    // ===== Draw & Save once ===== //
    saveFits(
        singleRatio,
        fitBg,
        "c_single_ratio_bg_fit",
        Form("PbPb 2.76 TeV | Single Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Sig/Mix",
        mode, plotXMin, plotXMax, plotYMin, plotYMax);

    // ===== Double Ratio ===== //
    TH1D* doubleRatio = histfuncRatio(singleRatio, bgFit, "dr_cor");

    saveRatio(doubleRatio,
        Form("dr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax,
        "dr_cor", "PbPb 2.76 TeV | Double Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Data/Fit",
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
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Data/Fit",
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
    std::string singleRatioBaseName = "sr";
    std::string doubleRatioBaseName = "dr";
    MultipleFitResults multiResults;
    
    for (const auto& f1 : fits1){
        std::string fitName = Form("%s_%s", singleRatioBaseName.c_str(), f1.displayName.c_str());
        multiResults.results.push_back(f1.result);
        multiResults.modelNames.push_back(fitName);
    }

    multiResults.results.push_back(bgRes);
    multiResults.modelNames.push_back("Background Fit");

    for (const auto& f : fits){
        std::string fitName = Form("%s_%s", doubleRatioBaseName.c_str(), f.displayName.c_str());
        multiResults.results.push_back(f.result);
        multiResults.modelNames.push_back(fitName);
    }
    // Cleanup memory
    for (auto& f : fits)
    delete f.function;
    delete bgFit;
    delete q_sigHist_cor;
    delete q_mixHist_cor;
    delete singleRatio;
    delete doubleRatio;

    return multiResults;
    
}

MultipleFitResults doubleRatioMixFitLevyAlphaOneSeed(
                    ControlVar selectedControlVar,
                    std::vector<std::pair<FitFunctionType, FitInit>> models,
                    double bin_low, double bin_high,
                    Double_t etaMin, Double_t etaMax,
                    Double_t ptMin,
                    Double_t q1, Double_t q2,
                    qMode mode,
                    double fitMin = 0.0, double fitMax = 10.0, double fitMinBG = 0.15,
                    double plotXMin = 0.0, double plotXMax = 10.0,
                    double plotYMin = 0.0, double plotYMax = 2.1,
                    const char* fileName = nullptr)
{
    const char* qmodename = (mode == qMode::QINV) ? "qinv" : "qlcms";
    const char* selectionVarName = getSelVarName(selectedControlVar);
    TString searchPattern;

    if (!fileName) {
        searchPattern = TString::Format(
            "./data/signal_mix/"
            "%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d*.root",
            qmodename, etaMin, (int)bin_low, (int)bin_high, (int)bin_low, (int)bin_high
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/"
            "%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d-%s.root",
            qmodename, etaMin, (int)bin_low, (int)bin_high, (int)bin_low, (int)bin_high, fileName
        );
    }

    const char* name_sigHist = "hSigSS";
    const char* name_mixHist = "hMixSS";

    TH1D* q_sigHist_cor = getHistogram(searchPattern.Data(), name_sigHist);
    TH1D* q_mixHist_cor = getHistogram(searchPattern.Data(), name_mixHist);

    TH1D* singleRatio = histhistRatioWithComp(
        q_sigHist_cor, q_mixHist_cor, q1, q2, "sr_cor", "testMix"
    );

    saveRatio(
        singleRatio,
        Form("sr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax,
        "sr_cor", "PbPb 2.76 TeV | Single Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Sig/Mix",
        "Single Ratio (Sig/Mix)"
    );

    FitInit levyInit = findModelInitOrDefault(
        models,
        FitFunctionType::DOUBLE_LEVY,
        {{0.6, 4.0, 1.0, 0.0, 1.0}}
    );

    FitOutput singleRatioLevyFit = fitDoubleLevyAlphaOneThenFree(
        singleRatio,
        levyInit,
        fitMin,
        fitMax
    );

    saveFits(
        singleRatio,
        {singleRatioLevyFit},
        "c_single_ratio_levy_alpha1_seed_fit",
        Form("PbPb 2.76 TeV | Single Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_alpha1seed_sr_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Sig/Mix",
        mode, plotXMin, plotXMax, plotYMin, plotYMax
    );

    FitInit bgInit{{1.0, 0.1, 2.0, 0.1, 2.0}};
    TFitResultPtr bgRes;
    TF1* bgFit = fitHistogram(
        singleRatio,
        &bgRes,
        FitFunctionType::BACKGROUND,
        bgInit,
        fitMinBG
    );

    TH1D* doubleRatio = histfuncRatio(singleRatio, bgFit, "dr_cor");

    saveRatio(
        doubleRatio,
        Form("dr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax,
        "dr_cor", "PbPb 2.76 TeV | Double Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Data/Fit",
        "Double Ratio (Data/Fit)"
    );

    FitOutput doubleRatioLevyFit = fitDoubleLevyAlphaOneThenFree(
        doubleRatio,
        levyInit,
        fitMin,
        fitMax
    );

    saveFits(
        doubleRatio,
        {doubleRatioLevyFit},
        "c_double_ratio_levy_alpha1_seed_fit",
        Form("PbPb 2.76 TeV | Double Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_alpha1seed_dr_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C_{2}(q_{inv}) = Data/Fit",
        mode, plotXMin, plotXMax, plotYMin, plotYMax
    );

    MultipleFitResults multiResults;
    multiResults.results.push_back(singleRatioLevyFit.result);
    multiResults.modelNames.push_back("sr_Levy Fit (#alpha free, #alpha=1 seed)");
    multiResults.results.push_back(doubleRatioLevyFit.result);
    multiResults.modelNames.push_back("dr_Levy Fit (#alpha free, #alpha=1 seed)");

    delete singleRatioLevyFit.function;
    delete doubleRatioLevyFit.function;
    delete bgFit;
    delete q_sigHist_cor;
    delete q_mixHist_cor;
    delete singleRatio;
    delete doubleRatio;

    return multiResults;
}

void DeltaPhiDeltaEtaRatio(ControlVar selectedControlVar,
                    double bin_low, double bin_high,
                    double plotXMin, double plotXMax, 
                    double plotYMin, double plotYMax,
                    double plotZMin, double plotZMax,
                    Double_t etaCut,
                    const char* fileName = nullptr)
{   
    const char* selectionVarName = getSelVarName(selectedControlVar);
    TString searchPattern;

    if (!fileName) {
        searchPattern = TString::Format(
            "./data/signal_mix/DeltaEtaDeltaPhi_eta-abs0.95_pT-from-0.5_%s_%.0f-%.0f_cent%.0fto%.0f-*.root",
            selectionVarName, bin_low, bin_high, bin_low, bin_high
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/DeltaEtaDeltaPhi_eta-abs0.95_pT-from-0.5_%s_%.0f-%.0f_cent%.0fto%.0f-%s.root",
            selectionVarName, bin_low, bin_high, bin_low, bin_high, fileName
        );
    }

    const char* name_sigHist = Form("hSigSS_DeltaEtaDeltaPhi");
    const char* name_mixHist = Form("hMixSS_DeltaEtaDeltaPhi");
    
    TH2D* q_sigHist_cor = getHistogram2d(searchPattern.Data(), name_sigHist);
    TH2D* q_mixHist_cor = getHistogram2d(searchPattern.Data(), name_mixHist);
    
    // ===== Single ratio ===== //
    TH2D* singleRatio = histhistRatioDeltaEtaDeltaPhi(
        q_sigHist_cor, q_mixHist_cor, etaCut, "sr_cor"
    );

    saveRatio(singleRatio, 
        Form("sr_DeltaEtaDeltaPhi_%s_%.0f-%.0f", selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax, plotZMin, plotZMax,
        "sr_DeltaEtaDeltaPhi", "PbPb 2.76 TeV | Single Ratio",
        "; #Delta#eta; #Delta#varphi; C(#Delta#eta, #Delta#varphi) = Sig/Mix",
        "Single Ratio (Sig/Mix)"
    );

    // Cleanup memory
    delete q_sigHist_cor;
    delete q_mixHist_cor;
    delete singleRatio;
}

void QtQzQ0QCorrelationRatios(ControlVar selectedControlVar,
                    double bin_low, double bin_high,
                    double plotXMin, double plotXMax,
                    double plotYMin, double plotYMax,
                    double plotZMin, double plotZMax,
                    Double_t q1x, Double_t q2x,
                    Double_t q1y, Double_t q2y,
                    const char* fileName = nullptr)
{
    const char* selectionVarName = getSelVarName(selectedControlVar);
    TString searchPattern;

    if (!fileName) {
        searchPattern = TString::Format(
            "./data/signal_mix/*qtqzq0q*_%s_%.0f-%.0f*.root",
            selectionVarName, bin_low, bin_high
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/*qtqzq0q*_%s_%.0f-%.0f*%s*.root",
            selectionVarName, bin_low, bin_high, fileName
        );
    }

    // Load once and support both naming schemes:
    // 1) h_qtqz_SS_cor / h_qtqz_OS_cor / h_q0q_SS_cor / h_q0q_OS_cor
    // 2) hSigSS_qtqz / hMixSS_qtqz / hSigSS_q0q / hMixSS_q0q from buildQtqzq0q()
    TString dataFile = findFile(searchPattern);
    if (dataFile.IsNull()) {
        std::cerr << "Error: No data file found matching pattern: "
                  << searchPattern << std::endl;
        return;
    } else {
        std::cout << "Found data file: " << dataFile << std::endl;
    }

    TFile* file = TFile::Open(dataFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file: " << dataFile << std::endl;
        return;
    }

    auto clone2d = [&](const char* name) -> TH2D* {
        TH2D* h = dynamic_cast<TH2D*>(file->Get(name));
        if (!h) return nullptr;
        TH2D* hc = dynamic_cast<TH2D*>(h->Clone(Form("%s_clone", name)));
        hc->SetDirectory(nullptr);
        return hc;
    };

    TH2D* h_q0q_ss  = clone2d("h_q0q_SS_cor");
    TH2D* h_q0q_os  = clone2d("h_q0q_OS_cor");
    TH2D* h_q0q_mix  = clone2d("hMixSS_q0q");

    if (!h_q0q_ss)  h_q0q_ss  = clone2d("hSigSS_q0q");

    file->Close();
    delete file;

    if (!h_q0q_ss) {
        std::cerr << "Error: Required qtqz/q0q histograms not found in file." << std::endl;
        delete h_q0q_ss;
        delete h_q0q_os;
        delete h_q0q_mix;
        return;
    }

    bool has_os = (h_q0q_os);
    bool has_mix = (h_q0q_mix);

    TH2D* ratio_q0q = nullptr;
    const char* q0q_prefix = nullptr;
    const char* q0q_legend = nullptr;
    const char* q0q_title = nullptr;

    if (has_os) {
        ratio_q0q = histhistRatio2d(
            h_q0q_ss, h_q0q_os,
            q1x, q2x, q1y, q2y,
            "sr_q0q_cor"
        );
        q0q_prefix = "sr_q0q_cor";
        q0q_legend = "C(q_{0}^{2}, |q|^{2}) (SS/OS)";
        q0q_title = "; q_{0}^{2} [GeV^{2}]; |#vec{q}|^{2} [GeV^{2}]; C(q_{0}^{2}, |q|^{2}) = SS/OS";
    } else if (has_mix) {
        ratio_q0q = histhistRatio2d(
            h_q0q_ss, h_q0q_mix,
            q1x, q2x, q1y, q2y,
            "sr_q0q_cor"
        );
        q0q_prefix = "sr_q0q_cor";
        q0q_legend = "C(q_{0}^{2}, |q|^{2}) (Sig/Mix)";
        q0q_title = "; q_{0}^{2} [GeV^{2}]; |#vec{q}|^{2} [GeV^{2}]; C(q_{0}^{2}, |q|^{2}) = Sig/Mix";
    } else {
        ratio_q0q  = dynamic_cast<TH2D*>(h_q0q_ss->Clone("sig_q0q_cor"));
        q0q_prefix = "sig_q0q_cor";
        q0q_legend = "C(q_{0}^{2}, |q|^{2}) (Signal SS)";
        q0q_title = "; q_{0}^{2} [GeV^{2}]; |#vec{q}|^{2} [GeV^{2}]; C(q_{0}^{2}, |q|^{2})";
    }

    saveRatio(
        ratio_q0q,
        Form("%s_%s_%.0f-%.0f", q0q_prefix, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax, plotZMin, plotZMax,
           q0q_prefix, "PbPb 2.76 TeV | C(q_{0}^{2}, |q|^{2})",
           q0q_title, q0q_legend
    );

    delete h_q0q_ss;
    delete h_q0q_os;
    delete h_q0q_mix;
    delete ratio_q0q;
}
/*
void QtQzQ0QCorrelationRatios(ControlVar selectedControlVar,
                    double bin_low, double bin_high,
                    double plotXMin, double plotXMax,
                    double plotYMin, double plotYMax,
                    double plotZMin, double plotZMax,
                    Double_t q1x, Double_t q2x,
                    Double_t q1y, Double_t q2y,
                    const char* fileName = nullptr)
{
    const char* selectionVarName = getSelVarName(selectedControlVar);
    TString searchPattern;

    if (!fileName) {
        searchPattern = TString::Format(
            "./data/signal_mix/*qtqzq0q*_%s_%.0f-%.0f*.root",
            selectionVarName, bin_low, bin_high
        );
    } else {
        searchPattern = TString::Format(
            "./data/signal_mix/*qtqzq0q*_%s_%.0f-%.0f*%s*.root",
            selectionVarName, bin_low, bin_high, fileName
        );
    }

    // Load once and support both naming schemes:
    // 1) h_qtqz_SS_cor / h_qtqz_OS_cor / h_q0q_SS_cor / h_q0q_OS_cor
    // 2) hSigSS_qtqz / hMixSS_qtqz / hSigSS_q0q / hMixSS_q0q from buildQtqzq0q()
    TString dataFile = findFile(searchPattern);
    if (dataFile.IsNull()) {
        std::cerr << "Error: No data file found matching pattern: "
                  << searchPattern << std::endl;
        return;
    }

    TFile* file = TFile::Open(dataFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file: " << dataFile << std::endl;
        return;
    }

    auto clone2d = [&](const char* name) -> TH2D* {
        TH2D* h = dynamic_cast<TH2D*>(file->Get(name));
        if (!h) return nullptr;
        TH2D* hc = dynamic_cast<TH2D*>(h->Clone(Form("%s_clone", name)));
        hc->SetDirectory(nullptr);
        return hc;
    };

    TH2D* h_qtqz_ss = clone2d("h_qtqz_SS_cor");
    TH2D* h_qtqz_os = clone2d("h_qtqz_OS_cor");
    TH2D* h_q0q_ss  = clone2d("h_q0q_SS_cor");
    TH2D* h_q0q_os  = clone2d("h_q0q_OS_cor");
    TH2D* h_qtqz_mix = clone2d("hMixSS_qtqz");
    TH2D* h_q0q_mix  = clone2d("hMixSS_q0q");

    if (!h_qtqz_ss) h_qtqz_ss = clone2d("hSigSS_qtqz");
    if (!h_q0q_ss)  h_q0q_ss  = clone2d("hSigSS_q0q");

    file->Close();
    delete file;

    if (!h_qtqz_ss || !h_q0q_ss) {
        std::cerr << "Error: Required qtqz/q0q histograms not found in file." << std::endl;
        delete h_qtqz_ss;
        delete h_qtqz_os;
        delete h_qtqz_mix;
        delete h_q0q_ss;
        delete h_q0q_os;
        delete h_q0q_mix;
        return;
    }

    bool has_os = (h_qtqz_os && h_q0q_os);
    bool has_mix = (h_qtqz_mix && h_q0q_mix);

    TH2D* ratio_qtqz = nullptr;
    TH2D* ratio_q0q = nullptr;
    const char* qtqz_prefix = nullptr;
    const char* q0q_prefix = nullptr;
    const char* qtqz_legend = nullptr;
    const char* q0q_legend = nullptr;
    const char* qtqz_title = nullptr;
    const char* q0q_title = nullptr;

    if (has_os) {
        ratio_qtqz = histhistRatio2d(
            h_qtqz_ss, h_qtqz_os,
            q1x, q2x, q1y, q2y,
            "sr_qtqz_cor"
        );
        ratio_q0q = histhistRatio2d(
            h_q0q_ss, h_q0q_os,
            q1x, q2x, q1y, q2y,
            "sr_q0q_cor"
        );
        qtqz_prefix = "sr_qtqz_cor";
        q0q_prefix = "sr_q0q_cor";
        qtqz_legend = "C(q_{z}^{2}, q_{t}^{2}) (SS/OS)";
        q0q_legend = "C(q_{0}^{2}, |q|^{2}) (SS/OS)";
        qtqz_title = "; q_{t}^{2} [GeV^{2}]; q_{z}^{2} [GeV^{2}]; C(q_{z}^{2}, q_{t}^{2}) = SS/OS";
        q0q_title = "; q_{0}^{2} [GeV^{2}]; |#vec{q}|^{2} [GeV^{2}]; C(q_{0}^{2}, |q|^{2}) = SS/OS";
    } else if (has_mix) {
        ratio_qtqz = histhistRatio2d(
            h_qtqz_ss, h_qtqz_mix,
            q1x, q2x, q1y, q2y,
            "sr_qtqz_cor"
        );
        ratio_q0q = histhistRatio2d(
            h_q0q_ss, h_q0q_mix,
            q1x, q2x, q1y, q2y,
            "sr_q0q_cor"
        );
        qtqz_prefix = "sr_qtqz_cor";
        q0q_prefix = "sr_q0q_cor";
        qtqz_legend = "C(q_{z}^{2}, q_{t}^{2}) (Sig/Mix)";
        q0q_legend = "C(q_{0}^{2}, |q|^{2}) (Sig/Mix)";
        qtqz_title = "; q_{t}^{2} [GeV^{2}]; q_{z}^{2} [GeV^{2}]; C(q_{z}^{2}, q_{t}^{2}) = Sig/Mix";
        q0q_title = "; q_{0}^{2} [GeV^{2}]; |#vec{q}|^{2} [GeV^{2}]; C(q_{0}^{2}, |q|^{2}) = Sig/Mix";
    } else {
        ratio_qtqz = dynamic_cast<TH2D*>(h_qtqz_ss->Clone("sig_qtqz_cor"));
        ratio_q0q  = dynamic_cast<TH2D*>(h_q0q_ss->Clone("sig_q0q_cor"));
        qtqz_prefix = "sig_qtqz_cor";
        q0q_prefix = "sig_q0q_cor";
        qtqz_legend = "C(q_{z}^{2}, q_{t}^{2}) (Signal SS)";
        q0q_legend = "C(q_{0}^{2}, |q|^{2}) (Signal SS)";
        qtqz_title = "; q_{t}^{2} [GeV^{2}]; q_{z}^{2} [GeV^{2}]; C(q_{z}^{2}, q_{t}^{2})";
        q0q_title = "; q_{0}^{2} [GeV^{2}]; |#vec{q}|^{2} [GeV^{2}]; C(q_{0}^{2}, |q|^{2})";
    }

    saveRatio(
        ratio_qtqz,
        Form("%s_%s_%.0f-%.0f", qtqz_prefix, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax, plotZMin, plotZMax,
           qtqz_prefix, "PbPb 2.76 TeV | C(q_{z}^{2}, q_{t}^{2})",
           qtqz_title, qtqz_legend
    );

    saveRatio(
        ratio_q0q,
        Form("%s_%s_%.0f-%.0f", q0q_prefix, selectionVarName, bin_low, bin_high),
        plotXMin, plotXMax, plotYMin, plotYMax, plotZMin, plotZMax,
           q0q_prefix, "PbPb 2.76 TeV | C(q_{0}^{2}, |q|^{2})",
           q0q_title, q0q_legend
    );

    delete h_qtqz_ss;
    delete h_qtqz_os;
    delete h_qtqz_mix;
    delete h_q0q_ss;
    delete h_q0q_os;
    delete h_q0q_mix;
    delete ratio_qtqz;
    delete ratio_q0q;
}
*/
#endif
