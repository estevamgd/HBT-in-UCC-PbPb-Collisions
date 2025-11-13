#include "../include/fitSingleRatio.h"
#include <iostream>
#include <string>

int main() {
    // To compile use:
    // g++ -std=c++17 -pthread fitSingleRatio.cpp -o fitSingleRatio `root-config --cflags --libs`
    // To run use:
    // ./fitSingleRatio
    
    // --- Configuration ---
    const char* controlVarName = "CENTHF"; // "CENT", "MULT", or "CENTHF"
    double selVarMin = 3200.0;
    double selVarMax = 3300.0;
    double fitMin = 0.0;
    double fitMax = 0.1; 
    double plotXMin = 0.0; 
    double plotXMax = 0.1; 
    const std::string searchPath = "./data/correlation_ratios/";

    processRatio(searchPath, controlVarName, selVarMin, selVarMax, "Cor", "sr_cor", fitMin, fitMax, plotXMin, plotXMax);
    processRatio(searchPath, controlVarName, selVarMin, selVarMax, "Uncor", "sr_uncor", fitMin, fitMax, plotXMin, plotXMax);

    std::cout << "\n\nAll processing finished." << std::endl;
    return 0;
}