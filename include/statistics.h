#ifndef STATISTICS_H
#define STATISTICS_H

#include "TF1.h"
#include "TString.h"
#include "Math/ProbFunc.h"
#include "../include/data_func.h"
#include "../include/my_func.h" 
#include "../include/fits.h" 
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"


struct FTestSummary {
    std::string binLabel;
    std::string ratioName;
    double chi2Reduced;
    int ndfReduced;
    double chi2Full;
    int ndfFull;
    int numeratorDof;
    int denominatorDof;
    double msr;
    double mse;
    double fStatistic;
    double pValue;
};

inline FTestSummary fTestAlphaFixedVsFree(const std::string& binLabel,
                                          const std::string& ratioName,
                                          const FitOutput& reducedFit,
                                          const FitOutput& fullFit)
{
    FTestSummary summary;
    summary.binLabel = binLabel;
    summary.ratioName = ratioName;
    summary.chi2Reduced = reducedFit.result->Chi2();
    summary.ndfReduced = reducedFit.result->Ndf();
    summary.chi2Full = fullFit.result->Chi2();
    summary.ndfFull = fullFit.result->Ndf();
    summary.numeratorDof = summary.ndfReduced - summary.ndfFull;
    summary.denominatorDof = summary.ndfFull;
    summary.msr = std::numeric_limits<double>::quiet_NaN();
    summary.mse = std::numeric_limits<double>::quiet_NaN();
    summary.fStatistic = std::numeric_limits<double>::quiet_NaN();
    summary.pValue = std::numeric_limits<double>::quiet_NaN();

    if (summary.numeratorDof <= 0 || summary.denominatorDof <= 0) {
        std::cerr << "Cannot compute F-test for " << ratioName
                  << " in bin " << binLabel
                  << ": invalid degrees of freedom." << std::endl;
        return summary;
    }

    summary.msr =
        (summary.chi2Reduced - summary.chi2Full) / summary.numeratorDof;
    summary.mse = summary.chi2Full / summary.denominatorDof;

    if (summary.mse > 0.0) {
        summary.fStatistic = summary.msr / summary.mse;
    }

    if (std::isfinite(summary.fStatistic) && summary.fStatistic >= 0.0) {
        summary.pValue = ROOT::Math::fdistribution_cdf_c(
            summary.fStatistic,
            summary.numeratorDof,
            summary.denominatorDof
        );
    }

    std::cout << ratioName << " | F-test alpha fixed vs free"
              << " | MSR = " << summary.msr
              << ", MSE = " << summary.mse
              << ", F* = " << summary.fStatistic
              << " with dof = (" << summary.numeratorDof
              << ", " << summary.denominatorDof << ")"
              << ", p = " << summary.pValue << std::endl;

    return summary;
}

inline int findAlphaParameterIndex(const FitModelConfig& cfg)
{
    for (int i = 0; i < cfg.nPar; ++i) {
        if (cfg.parNames[i] == "#alpha") {
            return i;
        }
    }

    throw std::runtime_error("FitModelConfig has no #alpha parameter");
}

inline void f_test(TH1D* data,
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
    
    double fitMin = 0.02;
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

void plot_ftest() {
    // Define input and output paths based on your WSL environment
    std::string input_csv = "/home/estevamgd/projects/HBT/build/data/statistics/hfsumet_f_test_summary.csv";
    std::string output_dir = "/home/estevamgd/projects/HBT/build/imgs/test/statistics";
    
    // 1. Open the CSV file
    std::ifstream file(input_csv);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << input_csv << std::endl;
        std::cerr << "Make sure the path is correct and the file exists." << std::endl;
        return;
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    // Data containers
    std::vector<double> x_sr, F_sr, p_sr;
    std::vector<double> x_dr, F_dr, p_dr;

    // 2. Read and parse the CSV
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> cols;
        
        while (std::getline(ss, cell, ',')) {
            cols.push_back(cell);
        }

        if (cols.size() < 12) continue; // Safety check

        // Parse HFsumET_bin (e.g., "3200-3300") to get the bin center for the X-axis
        std::string bin_str = cols[0];
        size_t dash_pos = bin_str.find('-');
        double x_val = 0;
        if (dash_pos != std::string::npos) {
            double min_val = std::stod(bin_str.substr(0, dash_pos));
            double max_val = std::stod(bin_str.substr(dash_pos + 1));
            x_val = (min_val + max_val) / 2.0;
        }

        std::string ratio = cols[1];
        double f_stat = std::stod(cols[10]);
        double p_val  = std::stod(cols[11]);

        // Separate data into Single Ratio (sr) and Double Ratio (dr)
        if (ratio.find("sr") != std::string::npos) {
            x_sr.push_back(x_val);
            F_sr.push_back(f_stat);
            p_sr.push_back(p_val);
        } else if (ratio.find("dr") != std::string::npos) {
            x_dr.push_back(x_val);
            F_dr.push_back(f_stat);
            p_dr.push_back(p_val);
        }
    }
    file.close();

    // 3. Create TGraphs
    TGraph *gr_F_sr = new TGraph(x_sr.size(), &x_sr[0], &F_sr[0]);
    TGraph *gr_F_dr = new TGraph(x_dr.size(), &x_dr[0], &F_dr[0]);
    TGraph *gr_p_sr = new TGraph(x_sr.size(), &x_sr[0], &p_sr[0]);
    TGraph *gr_p_dr = new TGraph(x_dr.size(), &x_dr[0], &p_dr[0]);

    // 4. Formatting Graph Styles
    
    // F_statistic (SR): red, marker 24, line style 3
    gr_F_sr->SetMarkerStyle(24);
    gr_F_sr->SetMarkerColor(kRed);
    gr_F_sr->SetLineColor(kRed);
    gr_F_sr->SetLineStyle(3);
    gr_F_sr->SetLineWidth(2);

    // F_statistic (DR): red, marker 26, line style 7
    gr_F_dr->SetMarkerStyle(26);
    gr_F_dr->SetMarkerColor(kRed);
    gr_F_dr->SetLineColor(kRed);
    gr_F_dr->SetLineStyle(7);
    gr_F_dr->SetLineWidth(2);

    // p_value (SR): blue, marker 24, line style 3
    gr_p_sr->SetMarkerStyle(24);
    gr_p_sr->SetMarkerColor(kBlue);
    gr_p_sr->SetLineColor(kBlue);
    gr_p_sr->SetLineStyle(3);
    gr_p_sr->SetLineWidth(2);

    // p_value (DR): blue, marker 26, line style 7
    gr_p_dr->SetMarkerStyle(26);
    gr_p_dr->SetMarkerColor(kBlue);
    gr_p_dr->SetLineColor(kBlue);
    gr_p_dr->SetLineStyle(7);
    gr_p_dr->SetLineWidth(2);

    // 5. Drawing with TMultiGraph
    TCanvas *c1 = new TCanvas("c1", "F-Test Summary", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.08); // Ensure space for the CMS header

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr_F_sr, "LP"); // L = Line, P = Marker
    mg->Add(gr_F_dr, "LP");
    mg->Add(gr_p_sr, "LP");
    mg->Add(gr_p_dr, "LP");

    // Draw axes first before setting titles
    mg->Draw("A"); 
    mg->GetXaxis()->SetTitle("HFsumET bin");
    mg->GetYaxis()->SetTitle("F-statistic / p-value");
    
    mg->SetMinimum(0.0);

    // 6. Add a Legend
    TLegend *leg = new TLegend(0.65, 0.70, 0.90, 0.88); // Top Right
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(gr_F_sr, "F_{stat} (Single Ratio)", "lp");
    leg->AddEntry(gr_F_dr, "F_{stat} (Double Ratio)", "lp");
    leg->AddEntry(gr_p_sr, "p-value (Single Ratio)", "lp");
    leg->AddEntry(gr_p_dr, "p-value (Double Ratio)", "lp");
    leg->Draw();

    // 7. Add the Work in Progress Header (Top Left)
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "");

    // 8. Update and save
    c1->Modified();
    c1->Update();
    
    // Create the directory if it doesn't exist (kTRUE makes it act like 'mkdir -p')
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    
    // Save to the exact paths
    std::string pdf_path = output_dir + "/F_test_summary_graph.pdf";
    std::string png_path = output_dir + "/F_test_summary_graph.png";
    
    c1->SaveAs(pdf_path.c_str());
    c1->SaveAs(png_path.c_str());
    
    std::cout << "Successfully saved graphs to: \n" 
              << pdf_path << "\n" 
              << png_path << std::endl;
}

#endif
