#ifndef VALIDATION_FUNC_H
#define VALIDATION_FUNC_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <ctime>
#include <fstream>      
#include <sstream>      
#include <map>          
#include <cfloat>       
#include <numeric>      
#include <algorithm>
#include <TSystem.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include <TString.h>
#include <TColor.h>
#include <TGraph.h>
#include <TText.h>      
#include <TLatex.h>
#include "../include/my_func.h" 


std::map<std::string, double> compareHists(
    const TString& file1_path,
    const TString& file2_path,
    const std::vector<std::string>& hist_names1,
    const std::vector<std::string>& hist_names2,
    const TString& canvasTitle = "Histogram Comparison",
    const TString& outputBasePath = "./",
    const TString& outputPrefix = "comparison",
    const TString& yAxisTitle = "Ratio" 
) {
    std::map<std::string, double> average_ratios;
    
    // --- Input Validation ---
    if (hist_names1.size() != hist_names2.size()) {
        std::cerr << "Error: The two histogram name vectors must have the same size." << std::endl;
        return average_ratios;
    }
    if (hist_names1.empty()) {
        std::cout << "Info: Histogram name vectors are empty. Nothing to plot." << std::endl;
        return average_ratios;
    }

    TFile* file1 = TFile::Open(file1_path, "READ");
    TFile* file2 = TFile::Open(file2_path, "READ");

    // --- File Opening Validation ---
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error: Could not open file: " << file1_path << std::endl;
        if(file2) file2->Close();
        return average_ratios;
    }
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Could not open file: " << file2_path << std::endl;
        if(file1) file1->Close();
        return average_ratios;
    }

    // --- Canvas Setup ---
    int nHists = hist_names1.size();
    int nCols = static_cast<int>(ceil(sqrt(nHists)));
    int nRows = static_cast<int>(ceil(static_cast<double>(nHists) / nCols));

    TCanvas* canvas = new TCanvas("canvas", canvasTitle, 400 * nCols, 400 * nRows);
    canvas->Divide(nCols, nRows);
    gStyle->SetOptStat(0);

    for (int i = 0; i < nHists; ++i) {
        canvas->cd(i + 1);
        gPad->SetLeftMargin(0.15); 
        gPad->SetTopMargin(0.08);

        // --- Get Histograms ---
        TH1D* h1 = (TH1D*)file1->Get(hist_names1[i].c_str());
        TH1D* h2 = (TH1D*)file2->Get(hist_names2[i].c_str());

        // --- Histogram Validation ---
        if (!h1) {
            std::cerr << "Warning: Histogram '" << hist_names1[i] << "' not found in " << file1_path << ". Skipping." << std::endl;
            average_ratios[hist_names1[i]] = 0.0; // Record failure
            continue;
        }
        if (!h2) {
            std::cerr << "Warning: Histogram '" << hist_names2[i] << "' not found in " << file2_path << ". Skipping." << std::endl;
            average_ratios[hist_names1[i]] = 0.0; // Record failure
            continue;
        }

        // --- Ratio Calculation ---
        TString ratioName = TString::Format("ratio_%d", i);
        TH1D* h_ratio = (TH1D*)h1->Clone(ratioName);
        h_ratio->SetDirectory(0); 
        h_ratio->Sumw2();
        if(h2->GetSumw2N() == 0) h2->Sumw2();
        h_ratio->Divide(h2);

        // --- Auto-adjust Y-axis & Calculate Average Ratio ---
        double max_dev = 0.0;
        double total_ratio = 0.0;
        double num_valid_bins = 0.0;

        // This loop calculates the average ratio per bin (per q_inv)
        for (int bin = 1; bin <= h_ratio->GetNbinsX(); ++bin) {
            // Only consider bins where at least one histogram has entries
            if (h1->GetBinContent(bin) == 0 && h2->GetBinContent(bin) == 0) continue;
            
            double content = h_ratio->GetBinContent(bin);
            if (content == 0) continue; // Skip bins with no ratio

            double dev = std::abs(content - 1.0);
            if (dev > max_dev) {
                max_dev = dev;
            }
            total_ratio += content;
            num_valid_bins++;
        }
        
        // Store the average ratio for this histogram
        average_ratios[hist_names1[i]] = (num_valid_bins > 0) ? total_ratio / num_valid_bins : 0.0;

        if (max_dev < 0.05) { max_dev = 0.1; } // Set minimum zoom
        max_dev *= 1.1; // Add 10% margin
        h_ratio->GetYaxis()->SetRangeUser(1.0 - max_dev, 1.0 + max_dev);

        h_ratio->SetTitle("");
        h_ratio->GetYaxis()->SetTitle(yAxisTitle); 
        gStyle->SetPalette(kLake);
        h_ratio->Draw("HIST PLC");

        TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0); 
        line->SetLineColor(kRed); line->SetLineStyle(kDashed); line->Draw("SAME"); 

        TLegend* legend = new TLegend(0.60, 0.78, 0.88, 0.88); 
        legend->AddEntry(h_ratio, yAxisTitle, "pl"); 
        legend->AddEntry(line, "Ratio = 1.0", "l"); 
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->Draw("SAME"); 

        // Draw CMS Headers for each pad
        TString rightHeaderText = TString::Format("Ratio: %s", hist_names1[i].c_str());
        double sizeFactor = 1.0;
        // Reduce size if header text is long 
        if (rightHeaderText.Length() > 30) { 
            sizeFactor = 0.85; 
        } else if (rightHeaderText.Length() > 25) {
             sizeFactor = 0.9;
        }
        drawCMSHeaders("#bf{CMS} #it{Preliminary}", rightHeaderText.Data(), sizeFactor);
    }

    TCanvas* canvases[] = {canvas};
    save_canvas_images(canvases, 1, outputBasePath.Data(), outputPrefix.Data(), "png");
    save_canvas_images(canvases, 1, outputBasePath.Data(), outputPrefix.Data(), "pdf");

    file1->Close(); file2->Close();
    delete file1; delete file2; delete canvas;
    
    return average_ratios;
}


std::map<std::string, double> parseBenchmarkData(const std::string& filename) {
    std::map<std::string, double> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open benchmark file: " << filename << std::endl;
        return data;
    }
    std::string line;
    while (std::getline(file, line)) {
        size_t colon_pos = line.find(":");
        if (colon_pos == std::string::npos) continue;
        std::string key = line.substr(0, colon_pos);
        std::string value_str = line.substr(colon_pos + 1);
        std::stringstream ss(value_str);
        double value;
        if (ss >> value) {
            data[key] = value;
        }
    }
    file.close();
    return data;
}


void plot_timing_comparison(
    const std::vector<std::string>& method_names_expanded,
    const std::vector<std::vector<double>>& all_normalized_timings, 
    const std::vector<std::string>& centrality_labels,
    const std::vector<std::string>& time_labels_to_plot, 
    const char* output_filename_base = "timing_comparison_TGraph",
    ControlVar selectionVarType = ControlVar::CENT,
    const std::string& outputPath = "./imgs/test/timing_comparison/"
) {
    // --- Input Validation ---
    if (method_names_expanded.size() != all_normalized_timings.size()) { 
        std::cerr << "Erro: Tamanhos de vetores de entrada inconsistentes." << std::endl; return;
    }
    for (const auto& timings : all_normalized_timings) { 
        if (timings.size() != centrality_labels.size()) { 
            std::cerr << "Erro: Tamanho dos dados de tempo não corresponde ao número de bins." << std::endl; return;
        }
    }

    const int n_graphs_total = method_names_expanded.size();
    const int n_bins = centrality_labels.size();
    if (n_bins == 0) { std::cerr << "Erro: Não há dados para plotar." << std::endl; return; }

    double max_normalized_time = 0.0;
    double min_normalized_time = DBL_MAX; 
    std::vector<TGraph*> graphs;
    std::vector<double> x_values(n_bins);
    std::iota(x_values.begin(), x_values.end(), 0.5); 

    for (int i = 0; i < n_graphs_total; ++i) {
        const std::vector<double>& normalized_times = all_normalized_timings[i]; // Get pre-normalized data
        
        for (int j = 0; j < n_bins; ++j) {
            if (normalized_times[j] > max_normalized_time) max_normalized_time = normalized_times[j];
            if (normalized_times[j] < min_normalized_time && normalized_times[j] > 0) min_normalized_time = normalized_times[j];
        }
        graphs.push_back(new TGraph(n_bins, x_values.data(), normalized_times.data()));
    }

    TCanvas *c = new TCanvas("c_timing_comparison", "Timing", 1000, 700);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08); c->SetRightMargin(0.05);
    gStyle->SetOptStat(0); gStyle->SetPalette(kLake);

    const char * frameTitle; 
    TString rightHeaderText = "Timing"; 
    if (selectionVarType == ControlVar::CENTHF) frameTitle = ";Centrality (HFsumET);Time (s) / Event";
    else if (selectionVarType == ControlVar::MULT) frameTitle = ";Multiplicity (#Tracks);Time (s) / Event";
    else frameTitle = ";Centrality (%);Time (s) / Event";

    TH1F *frame = new TH1F("frame", frameTitle, n_bins, 0, n_bins);
    frame->SetStats(0);
    
    // Auto-scale Y-axis
    double y_min_range = min_normalized_time * 0.9;
    double y_max_range = max_normalized_time * 1.1;
    if (y_min_range == y_max_range) { y_min_range *= 0.9; y_max_range *= 1.1; }
    if (y_min_range == 0 && y_max_range == 0) y_max_range = 1.0;
    
    frame->GetYaxis()->SetRangeUser(y_min_range, y_max_range); 
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->SetRangeUser(0, n_bins);

    // Set X-axis labels
    for (int j = 0; j < n_bins; ++j) {
        frame->GetXaxis()->SetBinLabel(j + 1, centrality_labels[j].c_str());
    }

    TLegend *legend = new TLegend(0.65, 0.92 - n_graphs_total * 0.05, 0.93, 0.9);
    legend->SetBorderSize(0); legend->SetFillStyle(0); 

    std::vector<int> base_colors = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kCyan + 1};
    std::vector<int> line_styles = {1, 2, 7}; // 1=Solid, 2=Dashed, 7=LongDash
    std::vector<int> marker_styles = {20, 21, 22}; // 20=Circle, 21=Square, 22=Triangle

    int n_metrics = time_labels_to_plot.size();
    if (n_metrics == 0) n_metrics = 1;

    for (int i = 0; i < n_graphs_total; ++i) {
        int method_idx = i / n_metrics; // 0 for Serialized, 1 for Parallel
        int metric_idx = i % n_metrics; // 0 for Total, 1 for Signal, 2 for Mix
        graphs[i]->SetMarkerSize(1.2); graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerStyle(marker_styles[metric_idx % marker_styles.size()]);
        graphs[i]->SetLineStyle(line_styles[metric_idx % line_styles.size()]);
        graphs[i]->SetMarkerColor(base_colors[method_idx % base_colors.size()]);
        graphs[i]->SetLineColor(base_colors[method_idx % base_colors.size()]);
        legend->AddEntry(graphs[i], method_names_expanded[i].c_str(), "lp");
    }

    frame->Draw("AXIS"); 
    for (int i = 0; i < n_graphs_total; ++i) graphs[i]->Draw("LP SAME");
    legend->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", rightHeaderText.Data());

    TCanvas* canvases[1] = {c};
    save_canvas_images(canvases, 1, outputPath.c_str(), output_filename_base, "png");
    save_canvas_images(canvases, 1, outputPath.c_str(), output_filename_base, "pdf");

    delete c; delete frame; delete legend;
    for (auto g : graphs) delete g;
}


void plotTimingRatios(
    const std::vector<std::string>& ratio_metric_names, 
    const std::vector<std::vector<double>>& all_timing_ratios,
    const std::vector<std::string>& centrality_labels,
    const std::string& yAxisTitle, 
    const char* output_filename_base,
    ControlVar selectionVarType,
    const std::string& outputPath
) {
    const int n_graphs_total = ratio_metric_names.size();
    const int n_bins = centrality_labels.size();
    if (n_bins == 0) { std::cerr << "Erro: Não há dados para plotar ratios." << std::endl; return; }

    // --- Data & TGraph Creation ---
    double min_ratio = DBL_MAX;
    double max_ratio = -DBL_MAX;
    std::vector<TGraph*> graphs;
    std::vector<double> x_values(n_bins);

    for (int i = 0; i < n_graphs_total; ++i) {
        std::vector<double> ratio_values;
        for (int j = 0; j < n_bins; ++j) {
            if (i == 0) x_values[j] = j + 0.5;
            double ratio = all_timing_ratios[i][j];
            ratio_values.push_back(ratio);
            if (ratio == 0) continue; // Don't let 0 skew the y-axis range
            if (ratio < min_ratio) min_ratio = ratio;
            if (ratio > max_ratio) max_ratio = ratio;
        }
        graphs.push_back(new TGraph(n_bins, x_values.data(), ratio_values.data()));
    }

    TCanvas *c = new TCanvas("c_timing_ratios", "Timing Ratio", 1000, 700);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08); c->SetRightMargin(0.05);
    gStyle->SetOptStat(0); gStyle->SetPalette(kLake);

    TString frameTitle; // Title for axis labels only
    TString rightHeaderText = "Timing Ratio"; // Text for top-right header
    if (selectionVarType == ControlVar::CENTHF) frameTitle = TString::Format(";Centrality (HFsumET);%s", yAxisTitle.c_str());
    else if (selectionVarType == ControlVar::MULT) frameTitle = TString::Format(";Multiplicity (#Tracks);%s", yAxisTitle.c_str());
    else frameTitle = TString::Format(";Centrality (%%);%s", yAxisTitle.c_str());

    TH1F *frame = new TH1F("frame", frameTitle.Data(), n_bins, 0, n_bins);
    frame->SetStats(0);
    
    // Auto-scale Y-axis with 10% margin
    double y_min_range = min_ratio * 0.9;
    double y_max_range = max_ratio * 1.1;
    if (std::abs(y_min_range - y_max_range) < 1e-6) { 
         y_min_range *= 0.9; 
         y_max_range *= 1.1;
    }
    if (y_min_range == 0 && y_max_range == 0) y_max_range = 1.0;
    frame->GetYaxis()->SetRangeUser(y_min_range, y_max_range); 
    
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->SetRangeUser(0, n_bins);

    for (int j = 0; j < n_bins; ++j) {
        frame->GetXaxis()->SetBinLabel(j + 1, centrality_labels[j].c_str());
    }

    TLegend *legend = new TLegend(0.65, 0.92 - (n_graphs_total) * 0.05, 0.93, 0.9);
    legend->SetBorderSize(0); legend->SetFillStyle(0); 
    legend->SetHeader(TString::Format("Ratio: %s", yAxisTitle.c_str()), "C"); 

    std::vector<int> base_colors = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kCyan + 1};
    std::vector<int> marker_styles = {20, 21, 22};

    for (int i = 0; i < n_graphs_total; ++i) {
        graphs[i]->SetMarkerSize(1.2); graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerStyle(marker_styles[i % marker_styles.size()]);
        graphs[i]->SetLineStyle(1);
        graphs[i]->SetMarkerColor(base_colors[i % base_colors.size()]);
        graphs[i]->SetLineColor(base_colors[i % base_colors.size()]);
        legend->AddEntry(graphs[i], ratio_metric_names[i].c_str(), "lp");
    }
    
    frame->Draw("AXIS");
    for (int i = 0; i < n_graphs_total; ++i) {
        graphs[i]->Draw("LP SAME");
    }
    
    // Add TText labels for each point
    TText *txt = new TText();
    txt->SetTextSize(0.025);
    txt->SetTextAlign(12); 
    for (int i = 0; i < n_graphs_total; ++i) {
        int n_points = graphs[i]->GetN();
        double *x = graphs[i]->GetX();
        double *y = graphs[i]->GetY();
        for (int p = 0; p < n_points; ++p) {
            txt->SetTextColor(graphs[i]->GetLineColor());
            txt->DrawText(x[p] + 0.03, y[p], TString::Format("%.3f", y[p]));
        }
    }

    legend->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", rightHeaderText.Data());

    TCanvas* canvases[1] = {c};
    save_canvas_images(canvases, 1, outputPath.c_str(), output_filename_base, "png");
    save_canvas_images(canvases, 1, outputPath.c_str(), output_filename_base, "pdf");

    delete c; delete frame; delete legend; delete txt; 
    for (auto g : graphs) delete g;
}


void plotAverageRatios(
    const std::map<std::string, std::vector<double>>& all_avg_hist_ratios,
    const std::vector<std::string>& hist_names_in_order,
    const std::string& yAxisTitle,
    const char* output_filename_base,
    const std::string& outputPath
) {
    int n_bins = hist_names_in_order.size();
    if (n_bins == 0) { std::cerr << "Erro: Não há dados para plotar ratios de histograma." << std::endl; return; }

    TCanvas *c = new TCanvas("c_avg_hist_ratios", "Data Ratio", 1200, 700);
    c->SetLeftMargin(0.1); c->SetBottomMargin(0.25);
    c->SetTopMargin(0.08); c->SetRightMargin(0.05);
    gStyle->SetOptStat(0); gStyle->SetPalette(kLake);

    TString frameTitle = TString::Format(";Hist.;%s", yAxisTitle.c_str()); // Changed X-axis title
    TString rightHeaderText = "Data Ratio"; // Text for top-right header
    
    TH1F *frame = new TH1F("h_avg_hist_frame", frameTitle, n_bins, 0, n_bins);
    TGraph *g_summary = new TGraph();
    
    double max_dev = 0.0;
    
    // --- Data & TGraph Creation ---
    for (int i = 0; i < n_bins; ++i) {
        const std::string& hist_name = hist_names_in_order[i];
        frame->GetXaxis()->SetBinLabel(i + 1, hist_name.c_str());
        
        if (all_avg_hist_ratios.find(hist_name) == all_avg_hist_ratios.end()) {
            std::cerr << "Warning: No ratio data found for histogram " << hist_name << std::endl;
            g_summary->SetPoint(i, i + 0.5, 0.0);
            continue;
        }

        const std::vector<double>& ratios = all_avg_hist_ratios.at(hist_name);
        double sum = std::accumulate(ratios.begin(), ratios.end(), 0.0);
        double avg = (ratios.size() > 0) ? sum / ratios.size() : 0.0;
        
        g_summary->SetPoint(i, i + 0.5, avg);
        
        if (avg == 0) continue;
        double dev = std::abs(avg - 1.0);
        if (dev > max_dev) max_dev = dev;
    }

    // Auto-scale Y-axis
    if (max_dev < 0.05) max_dev = 0.1;
    max_dev *= 1.1;
    frame->GetYaxis()->SetRangeUser(1.0 - max_dev, 1.0 + max_dev);
    
    // X-Axis Labels
    frame->GetXaxis()->SetLabelSize(0.03);
    frame->GetXaxis()->LabelsOption("v");
    frame->GetXaxis()->CenterTitle(false);
    frame->GetXaxis()->SetTitleOffset(1.2);
    
    frame->Draw("AXIS");

    TLine* line = new TLine(0, 1.0, n_bins, 1.0); 
    line->SetLineColor(kRed); line->SetLineStyle(kDashed); line->SetLineWidth(2);
    line->Draw("SAME");

    g_summary->SetMarkerStyle(20);
    g_summary->SetMarkerSize(1.5);
    g_summary->SetMarkerColor(kBlue + 1);
    g_summary->SetLineColor(kBlue + 1);
    g_summary->Draw("P SAME");

    TLegend *legend = new TLegend(0.7, 0.8, 0.93, 0.9);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(g_summary, "Avg. Ratio (All Bins)", "p");
    legend->AddEntry(line, "Ratio = 1.0", "l");
    legend->Draw("SAME");
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", rightHeaderText.Data());

    TCanvas* canvases[1] = {c};
    save_canvas_images(canvases, 1, outputPath.c_str(), output_filename_base, "png");
    save_canvas_images(canvases, 1, outputPath.c_str(), output_filename_base, "pdf");

    delete c; delete frame; delete g_summary; delete legend; delete line;
}


void runFullValidation(
    ControlVar methodSelectionVar,
    const std::vector<double>& centralityBins,
    const std::string& histDataPath,
    const std::string& histOutputPath,
    const std::vector<std::string>& hists_to_compare,
    const std::string& histFilePrefixRegular, 
    const std::string& histFilePrefixNew,     
    const std::string& histRatioYAxisTitle,
    const std::string& benchmarkPath,
    const std::vector<std::string>& method_names,       
    const std::vector<std::string>& method_prefixes,    
    const std::vector<std::string>& time_labels_to_plot, 
    const std::string& timingOutputPath,
    const std::string& timing_output_filename_base
) {
    std::cout << "\n=========================================" << std::endl;
    std::cout << "STARTING: Full Validation Run" << std::endl;
    std::cout << "=========================================\n" << std::endl;

    const char* selVarName = getSelVarName(methodSelectionVar);

    // --- Data structures for timing plot ---
    int n_metrics = time_labels_to_plot.size();
    int n_methods = method_names.size();
    std::vector<std::vector<double>> all_timings_expanded(n_methods * n_metrics);
    std::vector<std::string> method_names_expanded;
    std::vector<double> n_events_sig; 
    std::vector<double> n_events_mix; 
    std::vector<std::string> centrality_labels;

    // --- Data structure for Avg. Hist. Ratio plot ---
    std::map<std::string, std::vector<double>> all_avg_hist_ratios;


    // --- Main Loop over Centrality Bins ---
    for (size_t i = 0; i < centralityBins.size() - 1; ++i) {
        double bin_low = centralityBins[i];
        double bin_high = centralityBins[i+1];
        
        TString bin_label = TString::Format("%.0f-%.0f", bin_low, bin_high);
        if (methodSelectionVar == ControlVar::CENTHF) bin_label += " HF";
        
        std::cout << "\n--- Processing Bin: " << bin_label.Data() << " ---" << std::endl;

        // --- HISTOGRAM RATIO VALIDATION ---
        std::cout << "  (1) Running Histogram Comparison..." << std::endl;
        
        TString pattern1 = TString::Format("%s/%s_%s_%f-%f*.root",
                                          histDataPath.c_str(), histFilePrefixRegular.c_str(), 
                                          selVarName, bin_low, bin_high);
        TString pattern2 = TString::Format("%s/%s_%s_%f-%f*.root",
                                          histDataPath.c_str(), histFilePrefixNew.c_str(), 
                                          selVarName, bin_low, bin_high);

        TString file1 = findFile(pattern1, "");
        TString file2 = findFile(pattern2, "");

        if (file1.IsNull() || file2.IsNull()) {
            std::cerr << "  Error: Could not find histogram files for this bin. Skipping ratio plot." << std::endl;
            if (file1.IsNull()) std::cerr << "  -> Pattern not found: " << pattern1 << std::endl;
            if (file2.IsNull()) std::cerr << "  -> Pattern not found: " << pattern2 << std::endl;
        } else {
            std::cout << "  -> Comparing: " << file1 << std::endl;
            std::cout << "  -> With:      " << file2 << std::endl;
            
            TString commonTitle = TString::Format("Comparison (%s %.0f-%.0f)", selVarName, bin_low, bin_high);
            TString outPrefix = TString::Format("%s_vs_%s_%s_%.0f-%.0f", 
                                               histFilePrefixRegular.c_str(), histFilePrefixNew.c_str(), 
                                               selVarName, bin_low, bin_high);
            
            // Store the returned average ratios
            std::map<std::string, double> bin_avg_ratios = compareHists(
                file1, file2, hists_to_compare, hists_to_compare, 
                commonTitle, histOutputPath.c_str(), outPrefix, histRatioYAxisTitle.c_str()
            );

            // Populate the master map
            for(const auto& pair : bin_avg_ratios) {
                all_avg_hist_ratios[pair.first].push_back(pair.second);
            }
        }

        // --- TIMING BENCHMARK DATA COLLECTION ---
        std::cout << "  (2) Collecting Timing Data..." << std::endl;
        centrality_labels.push_back(bin_label.Data());
        bool n_events_fetched = false;

        for (size_t j = 0; j < method_prefixes.size(); ++j) {
            TString pattern = TString::Format("%s/%s_%s_%f-%f_benchmark-*.txt",
                                              benchmarkPath.c_str(),
                                              method_prefixes[j].c_str(),
                                              selVarName,
                                              bin_low,
                                              bin_high);
            TString file = findFile(pattern, ""); 
            if (file.IsNull()) {
                std::cerr << "  -> ERROR: Benchmark file not found for method '" << method_names[j] << "'" << std::endl;
                std::cerr << "     Pattern: " << pattern << std::endl;
                for (size_t k = 0; k < time_labels_to_plot.size(); ++k) {
                    int expanded_idx = j * n_metrics + k;
                    all_timings_expanded[expanded_idx].push_back(0.0);
                    if (i == 0) { 
                         method_names_expanded.push_back(method_names[j] + " - " + time_labels_to_plot[k]);
                    }
                }
                continue;
            }
            
            std::cout << "  -> Found '" << method_names[j] << "' file: " << file << std::endl;
            auto data = parseBenchmarkData(file.Data());

            for (size_t k = 0; k < time_labels_to_plot.size(); ++k) {
                const std::string& label = time_labels_to_plot[k];
                int expanded_idx = j * n_metrics + k;
                if (i == 0) {
                    method_names_expanded.push_back(method_names[j] + " - " + label);
                }
                if (data.find(label) != data.end()) {
                    all_timings_expanded[expanded_idx].push_back(data[label]);
                } else {
                    std::cerr << "     Warning: Metric '" << label << "' not found in " << file << ". Using 0.0." << std::endl;
                    all_timings_expanded[expanded_idx].push_back(0.0);
                }
            }
            
            if (!n_events_fetched) {
                if (data.count("Number of Processed Events(Sig)")) {
                    n_events_sig.push_back(data.at("Number of Processed Events(Sig)"));
                } else {
                    std::cerr << "     Warning: 'Number of Processed Events(Sig)' not found in " << file << ". Using 0." << std::endl;
                    n_events_sig.push_back(0);
                }
                if (data.count("Number of Processed Events(Mix)")) {
                    n_events_mix.push_back(data.at("Number of Processed Events(Mix)"));
                } else {
                    std::cerr << "     Warning: 'Number of Processed Events(Mix)' not found in " << file << ". Using 0." << std::endl;
                    n_events_mix.push_back(0);
                }
                n_events_fetched = true;
            }
        } 
    } 

    
    // --- PLOT TIMING AND RATIO GRAPHS ---
    if (n_events_sig.empty()) { 
        std::cerr << "\nError: No timing data was successfully loaded. Aborting timing plots." << std::endl;
    } else {
        
        // --- PRE-NORMALIZATION OF TIMING DATA ---
        std::cout << "\n--- Normalizing Timing Data ---" << std::endl;
        std::vector<std::vector<double>> all_timings_normalized(n_methods * n_metrics);
        for (size_t j = 0; j < n_methods; ++j) {
            for (size_t k = 0; k < n_metrics; ++k) {
                int expanded_idx = j * n_metrics + k;

                    // ===================================================================
                    // --- Always normalize by number of signal events ---
                    const std::string& label = time_labels_to_plot[k];
                    std::cout << "  -> Normalizing '" << label << "' using nEvents(Sig)" << std::endl;
                    const std::vector<double>& n_events_to_use = n_events_sig;
                    // ===================================================================

                const std::vector<double>& raw_times = all_timings_expanded[expanded_idx];

                if (raw_times.size() != centrality_labels.size()) {
                    std::cerr << "FATAL ERROR: Raw times vector size mismatch for " << method_names_expanded[expanded_idx] << std::endl;
                    continue; 
                }
                if (n_events_to_use.size() != centrality_labels.size()) {
                    std::cerr << "FATAL ERROR: N_events vector size mismatch for " << method_names_expanded[expanded_idx] << std::endl;
                    continue;
                }

                for (size_t bin_idx = 0; bin_idx < centrality_labels.size(); ++bin_idx) {
                    double norm_time = (n_events_to_use[bin_idx] > 0) ? raw_times[bin_idx] / n_events_to_use[bin_idx] : 0;
                    all_timings_normalized[expanded_idx].push_back(norm_time);
                }
            }
        }
        // --- END OF PRE-NORMALIZATION ---


        // --- Plot Absolute Timings ---
        std::cout << "\n--- Plotting Absolute Timing Comparison ---" << std::endl;
        plot_timing_comparison(
            method_names_expanded, 
            all_timings_normalized, // Pass normalized data
            centrality_labels, time_labels_to_plot, 
            timing_output_filename_base.c_str(), 
            methodSelectionVar, timingOutputPath
        );
        std::cout << "Absolute timing plot finished." << std::endl;

        // --- Calculate and Plot Timing Ratios ---
        if (n_methods == 2) { 
            std::cout << "\n--- Plotting Timing Ratios ---" << std::endl;
            std::vector<std::vector<double>> all_timing_ratios(n_metrics);
            std::vector<std::string> ratio_metric_names;
            std::string ratio_y_title = method_names[0] + " / " + method_names[1];

            for (int k = 0; k < n_metrics; ++k) {
                std::vector<double>& regular_data = all_timings_expanded[k]; 
                std::vector<double>& new_data = all_timings_expanded[n_metrics + k]; 
                ratio_metric_names.push_back(time_labels_to_plot[k] + " Ratio");
                
                for (size_t j = 0; j < regular_data.size(); ++j) {
                    if (new_data[j] > 0) {
                        all_timing_ratios[k].push_back(regular_data[j] / new_data[j]);
                    } else {
                        all_timing_ratios[k].push_back(0.0); // Avoid division by zero
                    }
                }
            }
            
            std::string ratio_filename = timing_output_filename_base + "_Ratios";
            plotTimingRatios(
                ratio_metric_names, all_timing_ratios, centrality_labels,
                ratio_y_title, ratio_filename.c_str(),
                methodSelectionVar, timingOutputPath
            );
            std::cout << "Timing ratio plot finished." << std::endl;

        } else {
            std::cout << "\n--- Skipping Timing Ratios Plot (not 2 methods) ---" << std::endl;
        }

        // --- Write Summary TXT File ---
        std::string txt_summary_path = timingOutputPath + "/validation_summary.txt";
        std::ofstream txt_summary(txt_summary_path);
        if (!txt_summary.is_open()) {
            std::cerr << "Error: Could not open summary file for writing: " << txt_summary_path << std::endl;
        } else {
            txt_summary << "# Validation Summary\n";
            txt_summary << "# Timing (s/event) for each method and centrality bin\n";
            txt_summary << "# Centrality bins: " << centrality_labels.size() << "\n";
            txt_summary << "# Methods: ";
            for (const auto& name : method_names_expanded) txt_summary << name << ", ";
            txt_summary << "\n\n";

            // Timing per method/metric/centrality
            txt_summary << "Timing (s/event):\n";
            for (size_t i = 0; i < centrality_labels.size(); ++i) {
                txt_summary << "Centrality: " << centrality_labels[i] << "\n";
                for (size_t idx = 0; idx < method_names_expanded.size(); ++idx) {
                    txt_summary << "  " << method_names_expanded[idx] << ": " << all_timings_normalized[idx][i] << "\n";
                }
                txt_summary << "\n";
            }

            // Timing ratios (if two methods)
            if (n_methods == 2) {
                txt_summary << "Timing Ratios (Serialized/Parallel):\n";
                for (int k = 0; k < n_metrics; ++k) {
                    txt_summary << "  Metric: " << time_labels_to_plot[k] << "\n";
                    for (size_t j = 0; j < centrality_labels.size(); ++j) {
                        double ratio = 0.0;
                        double denom = all_timings_normalized[n_metrics + k][j];
                        if (denom > 0)
                            ratio = all_timings_normalized[k][j] / denom;
                        txt_summary << "    " << centrality_labels[j] << ": " << ratio << "\n";
                    }
                }
                txt_summary << "\n";
            }

            // Average ratios for each histogram
            txt_summary << "Average Ratios for Each Histogram:\n";
            for (const auto& hist_name : hists_to_compare) {
                txt_summary << "  " << hist_name << ": ";
                if (all_avg_hist_ratios.count(hist_name)) {
                    const auto& ratios = all_avg_hist_ratios.at(hist_name);
                    for (size_t i = 0; i < ratios.size(); ++i) {
                        txt_summary << ratios[i];
                        if (i < ratios.size() - 1) txt_summary << ", ";
                    }
                } else {
                    txt_summary << "No data";
                }
                txt_summary << "\n";
            }
            txt_summary.close();
            std::cout << "Summary TXT file written to: " << txt_summary_path << std::endl;
        }

        // --- Write Summary CSV File (tab-delimited, block format) ---
        std::string csv_summary_path = timingOutputPath + "/validation_summary.csv";
        std::ofstream csv_summary(csv_summary_path);
        if (!csv_summary.is_open()) {
            std::cerr << "Error: Could not open CSV file for writing: " << csv_summary_path << std::endl;
        } else {
            // Timing block per centrality
            for (size_t i = 0; i < centrality_labels.size(); ++i) {
                csv_summary << "Centrality\t" << centrality_labels[i] << "\n";
                for (size_t idx = 0; idx < method_names_expanded.size(); ++idx) {
                    csv_summary << method_names_expanded[idx] << "\t" << all_timings_normalized[idx][i] << "\n";
                }
                csv_summary << "\n";
            }

            // Timing ratios block per centrality (if two methods)
            if (n_methods == 2) {
                for (size_t i = 0; i < centrality_labels.size(); ++i) {
                    csv_summary << "Centrality\t" << centrality_labels[i] << "\n";
                    for (int k = 0; k < n_metrics; ++k) {
                        double ratio = 0.0;
                        double denom = all_timings_normalized[n_metrics + k][i];
                        if (denom > 0)
                            ratio = all_timings_normalized[k][i] / denom;
                        csv_summary << "Timing Ratio - " << time_labels_to_plot[k] << "\t" << ratio << "\n";
                    }
                    csv_summary << "\n";
                }
            }

            // Average ratios for each histogram block
            csv_summary << "Average Ratios for Each Histogram\t\n";
            csv_summary << "Histogram";
            for (size_t i = 0; i < centrality_labels.size(); ++i) {
                csv_summary << "\t" << centrality_labels[i];
            }
            csv_summary << "\n";
            for (const auto& hist_name : hists_to_compare) {
                csv_summary << hist_name;
                if (all_avg_hist_ratios.count(hist_name)) {
                    const auto& ratios = all_avg_hist_ratios.at(hist_name);
                    for (size_t i = 0; i < centrality_labels.size(); ++i) {
                        if (i < ratios.size())
                            csv_summary << "\t" << ratios[i];
                        else
                            csv_summary << "\t";
                    }
                } else {
                    for (size_t i = 0; i < centrality_labels.size(); ++i) {
                        csv_summary << "\t";
                    }
                }
                csv_summary << "\n";
            }
            csv_summary.close();
            std::cout << "Summary CSV file written to: " << csv_summary_path << std::endl;
        }
    }
    
    // --- PLOT AVERAGE HISTOGRAM RATIOS ---
    if (all_avg_hist_ratios.empty()) {
        std::cerr << "\nError: No histogram ratio data was collected. Aborting average ratio plot." << std::endl;
    } else {
        std::cout << "\n--- Plotting Average Histogram Ratios ---" << std::endl;
        std::string avg_hist_ratio_filename = "summary_hist_ratios";
        plotAverageRatios(
            all_avg_hist_ratios, hists_to_compare, histRatioYAxisTitle,
            avg_hist_ratio_filename.c_str(), histOutputPath
        );
        std::cout << "Average histogram ratio plot finished." << std::endl;
    }


    std::cout << "\n=========================================" << std::endl;
    std::cout << "Validation run finished." << std::endl;
    std::cout << "=========================================\n" << std::endl;
}

#endif