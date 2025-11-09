#include <iostream>
#include <vector>
#include <string>
#include <filesystem>     
#include "../include/validation_func.h" 


int main() {
    // To compile use:
    // g++ -std=c++17 -pthread runValidation.cpp -o runValidation `root-config --cflags --libs`
    // To run use:
    // ./runValidation
    std::string runTimestamp = getTimestamp();

    // ---COMMON CONFIGURATION
    ControlVar methodSelectionVar = ControlVar::CENTHF;
    std::vector<double> centralityBins = {3200, 3300}; 

    // --- HISTOGRAM RATIO CONFIGURATION
    const std::string histDataPath = "./data/signal_mix/";
    const std::string baseHistOutputPath = "./imgs/test/validation/";
    const std::string histOutputPath = baseHistOutputPath + runTimestamp + "/";
    const std::vector<std::string> hists_to_compare = {
        "h_qinvSS_signal_2l", "h_qinvSSCor_signal_2l",
        "h_qinvOS_signal_2l", "h_qinvOSCor_signal_2l",
        "h_qinvDiv_signal_2l", "h_qinvDivCor_signal_2l"
    };
    const std::string histFilePrefixRegular = "sig_double_loop";
    const std::string histFilePrefixNew = "sig_double_loop_parallel";
    const std::string histRatioYAxisTitle = "Serialized / Parallel";

    // --- TIMING BENCHMARK CONFIGURATION
    const std::string benchmarkPath = "benchmarks";
    const std::string baseTimingOutputPath = "./imgs/test/validation/";
    const std::string timingOutputPath = baseTimingOutputPath + runTimestamp + "/";
    const std::vector<std::string> method_names = {"Serialized", "Parallel"};
    const std::vector<std::string> method_prefixes = {"sig_double_loop", "sig_double_loop_parallel"};
    
    const std::vector<std::string> time_labels_to_plot = {
        "Total Time", 
        "Signal Time", 
        "Mix Time"
    }; 
    const std::string timing_output_filename_base = "timing_comparison_TGraph";
    
    std::cout << "Saving all outputs to subfolder: " << runTimestamp << std::endl;
    
    runFullValidation(
        methodSelectionVar,
        centralityBins,
        histDataPath,
        histOutputPath,
        hists_to_compare,
        histFilePrefixRegular,
        histFilePrefixNew,
        histRatioYAxisTitle,
        benchmarkPath,
        method_names,
        method_prefixes,
        time_labels_to_plot,
        timingOutputPath,
        timing_output_filename_base
    );
    
    return 0;
}

