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

void timing_comparison_plot_TGraph() {
    const Int_t n = 7;

    // Timmings in seconds
    double single_loop_ts[n] = {4248.484, 6127.988, 18661.691, 16563.013, 22991.744, 5719.090, 680.262}; // 0-0.02, 1-2, 2-10, 10-20, 20-60, 60-100, 100-140
    double double_loop_ts[n] = {5344.226, 5759.330, 23004.541, 18536.400, 22862.583, 5705.505, 660.923};
    double single_loop_utn_ts[n] = {4028.943, 6038.020, 23960.601, 20133.965, 31716.626, 5655.094, 694.461};
    
    // Normalization for hours/mins/seconds
    double seconds = 1.0; //3600.0
    
    // Number of events processed for each centrality range
    double  n_entries[n] = {2049, 2884, 13764, 16849, 61984, 63052, 62328};
    double xns[n];

    // Normalazing to Timming/event
    for (size_t i = 0; i < n; ++i) {
        single_loop_ts[i] = single_loop_ts[i] / (seconds * n_entries[i]);
        double_loop_ts[i] = double_loop_ts[i] / (seconds * n_entries[i]);
        single_loop_utn_ts[i] = single_loop_utn_ts[i] / (seconds * n_entries[i]);
        xns[i] = i;

        std::cout << "i: " << i << "| " << "single_loop: " << single_loop_ts[i] << " double_loop: " << double_loop_ts[i] << " single_loop_utn: " << single_loop_utn_ts[i] << " xns: " << xns[i] << std::endl;
    }

    TGraph *gr_single_loop_hpe = new TGraph(n, xns, single_loop_ts);
    TGraph *gr_double_loop_hpe = new TGraph(n, xns, double_loop_ts);
    TGraph *gr_upper_triangle_hpe = new TGraph(n, xns, single_loop_utn_ts);

    std::vector<std::string> centralities = {"0-0.02", "0.5-1", "1-5", "5-10", "10-30", "30-50", "50-70"}; // 0-0.02, 1-2, 2-10, 10-20, 20-60, 60-100, 100-140
    std::vector<std::string> mode_labels = {"Single Loop", "Double Loop", "Upper Triangle"};

    const int nBins = centralities.size();

    TH1F *frame = new TH1F("frame", "Timing Comparison;Centrality (%);Time/Event (s)", n, 0, n);
    // Style
    frame->SetStats(0);
    frame->GetYaxis()->SetRangeUser(0.0, 2.7); // Setting Y axis range
    frame->GetXaxis()->SetRangeUser(0.0, 6.0); // Setting X axis range
    gStyle->SetPalette(kLake);

    // Set bin labels 
    // Force ticks at bin edges
    frame->GetXaxis()->SetNdivisions(n, kTRUE);

    // Replace tick labels with your centrality labels
    for (int j = 0; j < nBins; ++j) {
        frame->GetXaxis()->ChangeLabel(j+1, -1, -1, -1, -1, -1, centralities[j].c_str());
    }

    // Center the text under the ticks
    frame->GetXaxis()->LabelsOption("C");
    frame->GetXaxis()->SetLabelOffset(0.01);  // adjust spacing below axis

    TCanvas *c = new TCanvas("c_sig_tgraph", "Timing Comparison", 1000, 600);
    c->SetGrid();
    c->SetLeftMargin(0.10);
    gStyle->SetOptStat(0);
    
    // Create legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.88, 0.88);
    legend->AddEntry(gr_single_loop_hpe, mode_labels[0].c_str(), "lp");
    legend->AddEntry(gr_double_loop_hpe, mode_labels[1].c_str(), "lp");
    legend->AddEntry(gr_upper_triangle_hpe, mode_labels[2].c_str(), "lp");
    legend->SetBorderSize(0);
    
    // More Style
    gr_single_loop_hpe->SetMarkerSize(1.);
    gr_single_loop_hpe->SetLineWidth(2.);
    gr_single_loop_hpe->SetMarkerStyle(20);
    gr_double_loop_hpe->SetMarkerSize(1.);
    gr_double_loop_hpe->SetLineWidth(2.);
    gr_double_loop_hpe->SetMarkerStyle(20);
    gr_upper_triangle_hpe->SetMarkerSize(1.);
    gr_upper_triangle_hpe->SetLineWidth(2.);
    gr_upper_triangle_hpe->SetMarkerStyle(20);

    // Get colors from the current palette
    int nColors = 3; // we have 3 graphs
    std::vector<int> colors(nColors);
    for (int i = 0; i < nColors; ++i) {
        double norm = double(i) / double(nColors - 1); // from 0 to 1
        colors[i] = TColor::GetColorPalette(int(norm * (gStyle->GetNumberOfColors() - 1)));
    }

    // Assign colors to graphs
    gr_single_loop_hpe->SetMarkerColor(colors[0]);
    
    gr_double_loop_hpe->SetMarkerColor(colors[1]);

    gr_upper_triangle_hpe->SetMarkerColor(colors[2]);
    
    // Draw
    frame->Draw();
    gr_single_loop_hpe->Draw("LP PLC SAME");
    gr_double_loop_hpe->Draw("LP PLC SAME");
    gr_upper_triangle_hpe->Draw("LP PLC SAME");
    legend->Draw();
    
    // Save the result
    TCanvas* canvases[1] = {c};
    const char *ipath = "./imgs/test/timing_comparison/";
    save_canvas_images(canvases, 1, ipath, "timing_comparison_TGraph", "png");
    save_canvas_images(canvases, 1, ipath, "timing_comparison_TGraph", "pdf");

    delete gr_single_loop_hpe;
    delete gr_double_loop_hpe;
    delete gr_upper_triangle_hpe;
}
