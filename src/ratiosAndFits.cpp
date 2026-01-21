#include "../include/ratios.h"
#include "../include/data_func.h"
#include "../include/my_func.h"
#include "../include/fits.h"

void doubleRatioFit(ControlVar selectedControlVar,
                    double bin_low, double bin_high,
                    Double_t q1, Double_t q2,
                    qMode mode,
                    double fitMin = 0.0, double fitMax = 10.0, double fitMinBG = 0.15,
                    double plotXMin = 0.0, double plotXMax = 10.0,
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

    saveRatio(singleRatio, "sr_cor", plotXMin, plotXMax, 
        "sr_cor", "PbPb 2.76 TeV | Single Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C(q_{LCMS}) = SS/OS"
            : "; q_{inv} [GeV]; C(q_{inv}) = SS/OS",
        "Single Ratio (SS/OS)"
    );

    // ===== Single Ratio Fits ===== //
    std::vector<std::pair<FitFunctionType, FitInit>> models1 = {
        {FitFunctionType::EXPONENTIAL, {{1.0, 0.5, 5.0, 0.0}}},
        {FitFunctionType::GAUSSIAN, {{1.0, 0.5, 5.0, 0.0}}},
        {FitFunctionType::LEVY, {{1.0, 0.5, 5.0, 0.0, 1.2}}},
        //{FitFunctionType::LEVY2, {{0.6, 4.0, 1.5}}},
        {FitFunctionType::DOUBLE_LEVY, {{0.6, 4.0, 1.5, 0.0, 1.0}}}
    };
    
    auto fits1 = fitHistogramMultiple(
        singleRatio,
        models1,
        fitMin, fitMax
    );

    // ===== Draw & Save once ===== //
    drawAndSaveFits(
        singleRatio,
        fits1,
        "c_single_ratio_fits",
        Form("PbPb 2.76 TeV | Single Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        mode, plotXMin, plotXMax);
    
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

    saveRatio(doubleRatio, "dr_cor", plotXMin, plotXMax, 
        "dr_cor", "PbPb 2.76 TeV | Double Ratio",
        (mode == qMode::QLCMS)
            ? "; q_{LCMS} [GeV]; C(q_{LCMS}) = Data/Fit"
            : "; q_{inv} [GeV]; C(q_{inv}) = Data/Fit",
        "Double Ratio (Data/Fit)"
    );

    // ===== Final Fits ===== //
    std::vector<std::pair<FitFunctionType, FitInit>> models = {
        {FitFunctionType::EXPONENTIAL, {{1.0, 0.5, 5.0, 0.0}}},
        {FitFunctionType::GAUSSIAN, {{1.0, 0.5, 5.0, 0.0}}},
        {FitFunctionType::LEVY, {{1.0, 0.5, 5.0, 0.0, 1.2}}},
        //{FitFunctionType::LEVY2, {{0.6, 4.0, 1.5}}},
        {FitFunctionType::DOUBLE_LEVY, {{0.6, 4.0, 1.5, 0.0, 1.0}}}
    };
    
    auto fits = fitHistogramMultiple(
        doubleRatio,
        models,
        fitMin, fitMax
    );
    
    // ===== Draw & Save once ===== //
    drawAndSaveFits(
        doubleRatio,
        fits,
        "c_double_ratio_fits",
        Form("PbPb 2.76 TeV | Double Ratio | %s: %.0f-%.0f", selectionVarName, bin_low, bin_high),
        Form("fit_%s_%s_%.0f-%.0f", qmodename, selectionVarName, bin_low, bin_high),
        mode, plotXMin, plotXMax);

    // Cleanup memory
    for (auto& f : fits)
    delete f.function;
    delete bgFit;
    delete q_ssHist_cor;
    delete q_osHist_cor;
    delete singleRatio;
    delete doubleRatio;
    
}


int main() {
    // to compile use:
    // g++ -std=c++17 -pthread ratiosAndFits.cpp -o ratiosAndFits `root-config --cflags --libs`
    // to run use:
    // ./ratiosAndFits

    // --- Configuration ---
    // Define the analysis parameters
    ControlVar selectedControlVar = ControlVar::CENTHF;
    qMode modeLCMS = qMode::QLCMS;
    qMode modeQinv = qMode::QINV;
    
    double bin_low = 3200.0;
    double bin_high = 3300.0;
    
    double plotXMin = 0.0;
    double plotXMax = 0.25;
    
    double fitMin = 0.05;
    double fitMax = 0.25;
    double fitMinBg = 0.2;

    // Normalization qinv range
    Double_t q1 = 4.82;
    Double_t q2 = 6.4;

    // --- Analysis
    doubleRatioFit(selectedControlVar, bin_low, bin_high, q1, q2, modeLCMS, fitMin, fitMax, fitMinBg, plotXMin, plotXMax);
    //doubleRatioFit(selectedControlVar, bin_low, bin_high, q1, q2, modeQinv, fitMin, fitMax, fitMinBg, plotXMin, plotXMax);

    return 0;
}