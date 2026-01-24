#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TText.h>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <filesystem>

#include "../include/my_func.h"

void timing_comparison_plot_TGraph(
    ControlVar methodSelectionVar,
    const std::vector<std::vector<Double_t>>& timings,
    const std::vector<Double_t>& n_events,
    const std::vector<std::string>& centrality_labels,
    const std::vector<std::string>& method_labels,
    const double normalization_factor = 1.0 // e.g., 1.0 for seconds, 60.0 for minutes
) {
    // Normalazing to Timming/event
    std::vector<std::vector<Double_t>> normalized_timings = timings;
    std::vector<std::vector<Double_t>> column_number(timings.size());

    for (size_t i = 0; i < timings.size(); ++i) {
        for (size_t j = 0; j < timings[i].size(); ++j) {
            normalized_timings[i][j] = timings[i][j] / (normalization_factor * n_events[j]);
            column_number[i].push_back(j + 0.5);
        }
    }

    std::vector<TGraph*> graphs;
    for (size_t i = 0; i < timings.size(); ++i) {
        graphs.push_back(new TGraph(timings[i].size(), column_number[i].data(), normalized_timings[i].data()));
    }
    
    const char* frameTitle = "Timing Comparison;Centrality (%);Time/Event (s)";
    if (methodSelectionVar == ControlVar::CENTHF) {
        frameTitle = "Timing Comparison; #Sigma E_{T} [GeV];Time / Event (s)";
    } else if (methodSelectionVar == ControlVar::MULT) {
        frameTitle = "Timing Comparison;Multiplicity (#Tracks);Time / Event (s)";
    }

    TH1F *frame = new TH1F("frame", frameTitle, centrality_labels.size(), 0, centrality_labels.size());
    // Style
    double x_max_range = centrality_labels.size();
    double y_max_range = 0.0;
    for (size_t i = 0; i < normalized_timings.size(); ++i) {
        for (size_t j = 0; j < normalized_timings[i].size(); ++j) {
            if (normalized_timings[i][j] > y_max_range) {
                y_max_range = normalized_timings[i][j];
            }
        }
    }
    y_max_range *= 1.1;
    
    frame->SetStats(0);
    frame->GetYaxis()->SetRangeUser(0.0, y_max_range); // Setting Y axis range
    frame->GetXaxis()->SetRangeUser(0.0, x_max_range); // Setting X axis range
    gStyle->SetPalette(kLake);

    // Set bin labels 
    // Force ticks at bin edges
    frame->GetXaxis()->SetNdivisions(centrality_labels.size(), kTRUE);

    // Replace tick labels with your centrality labels
    for (int i = 0; i < centrality_labels.size(); ++i) {
        frame->GetXaxis()->SetBinLabel(i + 1, centrality_labels[i].c_str());
    }

    // Center the text under the ticks
    frame->GetXaxis()->LabelsOption("C");
    //frame->GetXaxis()->SetLabelOffset(0.01);  // adjust spacing below axis

    TCanvas *c = new TCanvas("c_sig_tgraph", "Timing Comparison", 1000, 600);
    c->SetGrid();
    c->SetLeftMargin(0.10);
    gStyle->SetOptStat(0);
    
    // Create legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.88, 0.88);
    
    for (size_t i = 0; i < method_labels.size(); ++i) {
        legend->AddEntry(graphs[i], method_labels[i].c_str(), "lp");
    }
    legend->SetBorderSize(0);
    
    // More Style
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerSize(1.);
        graphs[i]->SetLineWidth(2.);
        graphs[i]->SetMarkerStyle(20);
    }

    // Get colors from the current palette
    int nColors = graphs.size(); // we have 3 graphs
    std::vector<int> colors(nColors);
    for (int i = 0; i < nColors; ++i) {
        double norm = double(i) / double(nColors - 1); // from 0 to 1
        colors[i] = TColor::GetColorPalette(int(norm * (gStyle->GetNumberOfColors() - 1)));
    }

    // Assign colors to graphs
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
    }
    
    // Draw
    frame->Draw();
    
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->Draw("LP PLC SAME");
    }
    legend->Draw();
    
    // Save the result
    TCanvas* canvases[1] = {c};
    const char *ipath = "./imgs/test/timing_comparison/";
    save_canvas_images(canvases, 1, ipath, "timing_comparison_TGraph", "png");
    save_canvas_images(canvases, 1, ipath, "timing_comparison_TGraph", "pdf");

    // Clean up
    delete c;
    delete frame;
    delete legend;
    for (auto g : graphs) delete g;
}

int main() {
    // To compile use:
    // g++ -std=c++17 -pthread timing_comparison_plot_TGraph.cpp -o timing_comparison_plot_TGraph `root-config --cflags --libs`
    // To run use:
    // ./timing_comparison_plot_TGraph
    
    ControlVar methodSelectionVar = ControlVar::CENT;
    std::vector<std::vector<Double_t>> timings = {
        {4248.484, 6127.988, 18661.691, 16563.013, 22991.744, 5719.090, 680.262}, // Method 1
        {5344.226, 5759.330, 23004.541, 18536.400, 22862.583, 5705.505, 660.923}  // Method 2
    };
    std::vector<Double_t> n_events = {2049, 2884, 13764, 16849, 61984, 63052, 62328};
    std::vector<std::string> centrality_labels = {"0-0.02", "0.5-1", "1-5", "5-10", "10-30", "30-50", "50-70"};
    std::vector<std::string> method_labels = {"Single Loop", "Double Loop"};

    timing_comparison_plot_TGraph(methodSelectionVar, timings, n_events, centrality_labels, method_labels, 1.0);

    /* ---- Metadata ---- */
    BenchmarkMetadata meta;
    meta.controlVariable = "CENTRALITY";
    meta.methodNames = method_labels;
    meta.inputFiles = {
        "PbPb_Run2011A.root",
        "PbPb_Run2011B.root"
    };

    /* ---- Per-bin entries ---- */
    std::vector<BenchmarkCentralityEntry> entries;

    for (size_t i = 0; i < centrality_labels.size(); ++i) {
        BenchmarkCentralityEntry e;
        e.label = centrality_labels[i];
        e.nEvents = n_events[i];

        for (size_t m = 0; m < timings.size(); ++m)
            e.times.push_back(timings[m][i]);

        entries.push_back(e);
    }

    save_benchmark_validation(meta, entries, "./benchmarks/", "timing_comparison", 1.0);

    return 0;
}