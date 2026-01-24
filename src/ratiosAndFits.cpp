#include "../include/ratiosAndFits.h"


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
    
    std::vector<std::pair<FitFunctionType, FitInit>> models = {
        {FitFunctionType::EXPONENTIAL, {{1.0, 0.5, 5.0, 0.0}}},
        {FitFunctionType::GAUSSIAN, {{1.0, 0.5, 5.0, 0.0}}},
        {FitFunctionType::LEVY2, {{0.6, 4.0, 1.5}}},
        {FitFunctionType::DOUBLE_LEVY, {{0.6, 4.0, 1.5, 0.0, 1.0}}}
    };

    double bin_low = 3200.0;
    double bin_high = 3300.0;
    
    double plotXMin = 0.0;
    double plotXMax = 0.3;

    double plotYMin = 0.9;
    double plotYMax = 2.1;
    
    double fitMin = 0.04;
    double fitMax = 0.25;
    double fitMinBg = 0.2;

    // Normalization qinv range
    Double_t q1 = 4.82;
    Double_t q2 = 6.4;

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

    // ===== Analysis =====
    doubleRatioFit(
        selectedControlVar,
        models,
        bin_low, bin_high,
        q1, q2,
        modeLCMS,
        fitMin, fitMax, fitMinBg,
        plotXMin, plotXMax, plotYMin, plotYMax
    );
    /*
    doubleRatioFit(
        selectedControlVar,
        models,
        bin_low, bin_high,
        q1, q2,
        modeQinv,
        fitMin, fitMax, fitMinBg,
        plotXMin, plotXMax, plotYMin, plotYMax
    );
    */
    return 0;
}