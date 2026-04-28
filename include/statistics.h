#include "TF1.h"
#include "TString.h"
#include "../include/data_func.h"
#include "../include/my_func.h" 
#include "../include/fits.h" 
#include <algorithm>
#include <stdexcept>

inline int findAlphaParameterIndex(const FitModelConfig& cfg)
{
    for (int i = 0; i < cfg.nPar; ++i) {
        if (cfg.parNames[i] == "#alpha") {
            return i;
        }
    }

    throw std::runtime_error("FitModelConfig has no #alpha parameter");
}

void f_test(TH1D* data, 
    FitModelConfig full_config, FitModelConfig reduced_config,
    const FitInit& full_init, const FitInit& reduced_init){
    // --- Configuration ---
    // Define the analysis parameters
    ControlVar selectedControlVar = ControlVar::CENTHF;
    qMode modeLCMS = qMode::QLCMS;
    qMode modeQinv = qMode::QINV;
    
    Double_t etaMin = 0.95;
    Double_t etaMax = -0.95;
    Double_t ptMin  = 0.50;

    std::vector<std::pair<FitFunctionType, FitInit>> models = {
        {FitFunctionType::EXPONENTIAL, {{0.6, 4.0, 0.0, 0.0, 1.0}}},
        {FitFunctionType::GAUSSIAN, {{0.6, 4.0, 0.0, 0.0, 1.0}}},
        {FitFunctionType::DOUBLE_LEVY, {{0.6, 4.0, 1.5, 0.0, 1.0}}}
    };

    double bin_low = 3700.0;
    double bin_high = 3800.0;
    
    double plotXMin = 0.0;
    double plotXMax = 0.5;

    double plotYMin = 0.9;
    double plotYMax = 2.1;
    
    double fitMin = 0.04;
    double fitMax = 0.2;
    double fitMinBg = 0.2;

    // Normalization qinv range
    Double_t q1 = 6.82;
    Double_t q2 = 8.4;

    TFitResultPtr fitResultFull, fitResultReduced;
    FitOutput full_fit_output, reduced_fit_output;

    if ((int)reduced_init.values.size() != reduced_config.nPar) {
        throw std::runtime_error("Wrong number of initial parameters for reduced fit");
    }
    if ((int)full_init.values.size() != full_config.nPar) {
        throw std::runtime_error("Wrong number of initial parameters for full fit");
    }

    int reducedAlphaIndex = findAlphaParameterIndex(reduced_config);
    int fullAlphaIndex = findAlphaParameterIndex(full_config);

    static int fixedAlphaFitCounter = 0;
    TString fixedAlphaName = TString::Format("fit_alpha1_seed_%d", fixedAlphaFitCounter++);
    TF1* fixedAlphaFit = new TF1(
        fixedAlphaName,
        reduced_config.func,
        fitMin,
        fitMax,
        reduced_config.nPar
    );

    for (int i = 0; i < reduced_config.nPar; ++i) {
        fixedAlphaFit->SetParameter(i, reduced_init.values[i]);
        fixedAlphaFit->SetParName(i, reduced_config.parNames[i].c_str());

        if (i == reducedAlphaIndex) {
            continue;
        }

        if (i < (int)reduced_config.parLimits.size() &&
            reduced_config.parLimits[i].first != reduced_config.parLimits[i].second) {
            fixedAlphaFit->SetParLimits(
                i,
                reduced_config.parLimits[i].first,
                reduced_config.parLimits[i].second
            );
        }
    }

    fixedAlphaFit->FixParameter(reducedAlphaIndex, 1.0);
    fitResultReduced = data->Fit(fixedAlphaFit, "S R E M");

    std::vector<double> seededFullValues = full_init.values;
    int nSeededParams = std::min(full_config.nPar, reduced_config.nPar);
    for (int i = 0; i < nSeededParams; ++i) {
        seededFullValues[i] = fixedAlphaFit->GetParameter(i);
    }
    seededFullValues[fullAlphaIndex] = 1.0;

    reduced_fit_output = {
        FitFunctionType::UNKNOWN,
        fixedAlphaFit,
        fitResultReduced,
        reduced_config.legendParams,
        reduced_config.displayName
    };

    full_fit_output = fitHistogramCustom(
        data,
        &fitResultFull,
        FitInit{seededFullValues},
        full_config,
        fitMin,
        fitMax
    );
    
}
