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
    
    double bin_low = 3200.0;
    double bin_high = 3300.0;
    
    double plotXMin = 0.0;
    double plotXMax = 0.3;
    
    double fitMin = 0.05;
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
        bin_low, bin_high,
        q1, q2,
        modeLCMS,
        fitMin, fitMax, fitMinBg,
        plotXMin, plotXMax,
        nullptr,   
        &scanCfgDR       
    );

    return 0;
}