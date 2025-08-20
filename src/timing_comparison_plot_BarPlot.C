#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <filesystem>

#include "../include/my_func.h"

void timing_comparison_plot_BarPlot() {
    // Original data in seconds
    //std::vector<int> single_loop_s = {6526, 5460, 17189}; // 0-1, 1-2, 2-10
    //std::vector<int> double_loop_s = {6388, 4790, 25842};
    //std::vector<int> single_loop_utn_s = {6875, 6252, 23754};
    
    //std::vector<int> single_loop_s = {35001, 21416, 32888, 6185, 754}; // 0-10, 10-20, 20-60, 60-100, 100-140
    //std::vector<int> double_loop_s = {28869, 16296, 21845, 5133, 608};
    //std::vector<int> single_loop_utn_s = {22172, 16104, 23445, 4362, 490};
    
    std::vector<int> single_loop_s = {5460, 17189, 21416, 32888, 6185, 754}; // 1-2, 2-10, 10-20, 20-60, 60-100, 100-140
    std::vector<int> double_loop_s = {4790, 25842, 16296, 21845, 5133, 608};
    std::vector<int> single_loop_utn_s = {6252, 23754, 16104, 23445, 4362, 490};
    
    // New data vectors in hours (double for precision)
    std::vector<double> single_loop_h(single_loop_s.size());
    std::vector<double> double_loop_h(double_loop_s.size());
    std::vector<double> upper_triangle_h(single_loop_utn_s.size());

    double seconds = 1.0; //3600.0

    for (size_t i = 0; i < single_loop_s.size(); ++i) {
        single_loop_h[i] = static_cast<double>(single_loop_s[i]) / seconds;
        double_loop_h[i] = static_cast<double>(double_loop_s[i]) / seconds;
        upper_triangle_h[i] = static_cast<double>(single_loop_utn_s[i]) / seconds;
    }

    // New data vectors in hours (double for precision)
    std::vector<double> single_loop_hpe(single_loop_s.size());
    std::vector<double> double_loop_hpe(double_loop_s.size());
    std::vector<double> upper_triangle_hpe(single_loop_utn_s.size());
    
    // Number of entries for each centrality
    std::vector<double>  n_entries = {2884, 13764, 16849, 61984, 63052, 62328};

    for (size_t i = 0; i < single_loop_s.size(); ++i) {
        single_loop_hpe[i] = static_cast<double>(single_loop_h[i]) / n_entries[i];
        double_loop_hpe[i] = static_cast<double>(double_loop_h[i]) / n_entries[i];
        upper_triangle_hpe[i] = static_cast<double>(upper_triangle_h[i]) / n_entries[i];
    }

    //std::vector<std::string> centralities = {"0-0.5", "0.5-1", "1-5"};
    std::vector<std::string> centralities = {"0-0.5", "1-5", "5-10", "10-30", "30-50", "50-70"};
    std::vector<std::string> mode_labels = {"Single Loop", "Double Loop", "Upper Triangle"};
    std::vector<int> colors = {kRed+1, kBlue+1, kGreen+2};

    const int nBins = centralities.size();

    TH1D* h1 = new TH1D("h1", "Timing Comparison;Centrality (%);Time/Event (s)", nBins, 0, nBins);
    TH1D* h2 = new TH1D("h2", "", nBins, 0, nBins);
    TH1D* h3 = new TH1D("h3", "", nBins, 0, nBins);

    // Set bin labels and fill with normalized data
    for (int j = 0; j < nBins; ++j) {
        h1->GetXaxis()->SetBinLabel(j+1, centralities[j].c_str());
        h1->SetBinContent(j+1, single_loop_hpe[j]);
        h2->SetBinContent(j+1, double_loop_hpe[j]);
        h3->SetBinContent(j+1, upper_triangle_hpe[j]);
    }

    // Set styles for bar charts
    h1->SetBarWidth(0.2); h1->SetBarOffset(0.15); h1->SetFillColor(colors[0]);
    h2->SetBarWidth(0.2); h2->SetBarOffset(0.40); h2->SetFillColor(colors[1]);
    h3->SetBarWidth(0.2); h3->SetBarOffset(0.65); h3->SetFillColor(colors[2]);

    //h1->GetYaxis()->SetRangeUser(0.0, 0.09);
    //h2->GetYaxis()->SetRangeUser(0.0, 0.09);
    //h3->GetYaxis()->SetRangeUser(0.0, 0.09);


    TCanvas *c = new TCanvas("c_sig_bar_norm", "Timing Comparison", 1000, 600);
    c->SetGrid();
    c->SetLeftMargin(0.10);
    gStyle->SetOptStat(0);

    h1->Draw("bar2");
    h2->Draw("bar2 same");
    h3->Draw("bar2 same");

    // Create legend
    TLegend *legend = new TLegend(0.15, 0.7, 0.32, 0.88);
    legend->AddEntry(h1, mode_labels[0].c_str(), "f");
    legend->AddEntry(h2, mode_labels[1].c_str(), "f");
    legend->AddEntry(h3, mode_labels[2].c_str(), "f");
    legend->SetBorderSize(0);
    legend->Draw();

    // Style
    gStyle->SetPalette(kLake);

    // Save the result
    TCanvas* canvases[1] = {c};
    const char *ipath = "./imgs/test/timing_comparison/";
    save_canvas_images(canvases, 1, ipath, "timing_comparison_bar", "png");
    save_canvas_images(canvases, 1, ipath, "timing_comparison_bar", "pdf");

    delete h1;
    delete h2;
    delete h3;
}
