#include "../include/ratiosAndFits.h"

void whatRatio(const char* rationame) {
    // -- Configuration --
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
    
    double fitMin = 0.02;
    double fitMax = 0.2;
    double fitMinBg = 0.2;

    // Normalization qinv range
    Double_t q1 = 6.82;
    Double_t q2 = 8.4;

    // ===== Fit-range scan configuration (Single Ratio) =====
    FitRangeScanConfig scanCfgSR;
    scanCfgSR.fitType = FitFunctionType::LEVY;

    scanCfgSR.fitMinLow  = 0.05;
    scanCfgSR.fitMinHigh = 0.10;
    scanCfgSR.fitMinStep = 0.005;

    scanCfgSR.fitMaxLow  = 0.10;
    scanCfgSR.fitMaxHigh = 0.30;
    scanCfgSR.fitMaxStep = 0.005;

    scanCfgSR.outDir = "./data/fit_range_scan/";

    // ===== Fit-range scan configuration (Double Ratio) =====
    FitRangeScanConfig scanCfgDR = scanCfgSR;

    // --- 2D configuration ---
    // --- Configuration ---
    // Define the analysis parameters
    ControlVar selectedControlVar2 = ControlVar::CENTHF;
    
    double bin_low2 = 3200.0;
    double bin_high2 = 3300.0;

    double plotXMin2d = 0.0;
    double plotXMax2d = 0.015;

    double plotYMin2d = 0.0;
    double plotYMax2d = 0.015;

    double plotZMin2d = 0.6;
    double plotZMax2d = 1.4;

    // Normalization qinv range
    Double_t etaCut2d = 0.04;

    Double_t q1x2d = 0.0;
    Double_t q2x2d = 0.015;
    Double_t q1y2d = 0.0;
    Double_t q2y2d = 0.015;

    std::cout << "Requested ratio: " << rationame << std::endl;
    if (strcmp(rationame, "singleRatio") == 0) {
        std::cout << "Performing Single Ratio fit..." << std::endl;
        // Call the function for single ratio fit
    } else if (strcmp(rationame, "doubleRatio") == 0) {
        std::cout << "Performing Double Ratio fit..." << std::endl;
        doubleRatioFit(
            selectedControlVar,
            models,
            bin_low, bin_high,
            etaMin, etaMax, ptMin,
            q1, q2,
            modeLCMS,
            fitMin, fitMax, fitMinBg,
            plotXMin, plotXMax, plotYMin, plotYMax,
            &scanCfgSR,
            &scanCfgDR
        );
    } else if (strcmp(rationame, "doubleRatioMix") == 0) {
        std::cout << "Performing Double Ratio fit..." << std::endl;
        doubleRatioMixFit(
            selectedControlVar,
            models,
            bin_low, bin_high,
            etaMin, etaMax, ptMin,
            q1, q2,
            modeLCMS,
            fitMin, fitMax, fitMinBg,
            plotXMin, plotXMax, plotYMin, plotYMax
        );
    } else if (strcmp(rationame, "DeltaEtaDeltaPhi") == 0) {
        std::cout << "Performing Delta Eta Delta Phi fit..." << std::endl;
        DeltaPhiDeltaEtaRatio(
            selectedControlVar2,
            bin_low2, bin_high2,
            plotXMin2d, plotXMax2d, plotYMin2d, plotYMax2d, plotZMin2d, plotZMax2d,
            etaCut2d
        );
    } else if (strcmp(rationame, "QtQzQ0Q") == 0) {
        std::cout << "Performing Qt Qz Q0 Q fit..." << std::endl;
           QtQzQ0QCorrelationRatios(
            selectedControlVar2,
            bin_low2, bin_high2,
            plotXMin2d, plotXMax2d, plotYMin2d, plotYMax2d, plotZMin2d, plotZMax2d,
            q1x2d, q2x2d, q1y2d, q2y2d
        );
    } else {
        std::cerr << "Unknown ratio name: " << rationame << std::endl;
    }
}

int main() {
    // to compile use:
    // g++ -std=c++17 -pthread ratiosAndFits.cpp -o ratiosAndFits `root-config --cflags --libs`
    // to run use:
    // ./ratiosAndFits
    // Possible arguments: singleRatio, doubleRatio, doubleRatioMix, DeltaEtaDeltaPhi, QtQzQ0Q
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    const char* rationames[] = {"singleRatio", "doubleRatioMix"};
    
    for (const char* rationame : rationames) {
        whatRatio(rationame);
    }
    return 0;
}
