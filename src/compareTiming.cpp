#include "../include/validation_func.h"


void plot_timing_comparison(
    const std::vector<std::string>& method_names,
    const std::vector<std::vector<double>>& all_timings,
    const std::vector<double>& n_events,
    const std::vector<std::string>& centrality_labels,
    const char* output_filename_base = "timing_comparison_TGraph",
    ControlVar selectionVarType = ControlVar::CENT
) {
    if (method_names.size() != all_timings.size() || n_events.size() != centrality_labels.size()) {
        std::cerr << "Erro: Tamanhos de vetores de entrada inconsistentes." << std::endl;
        return;
    }
    for (const auto& timings : all_timings) {
        if (timings.size() != n_events.size()) {
            std::cerr << "Erro: Tamanho dos dados de tempo não corresponde ao número de bins." << std::endl;
            return;
        }
    }

    const int n_methods = method_names.size();
    const int n_bins = centrality_labels.size();
    if (n_bins == 0) {
        std::cerr << "Erro: Não há dados para plotar." << std::endl;
        return;
    }

    double max_normalized_time = 0.0;
    double min_normalized_time = DBL_MAX; 

    std::vector<TGraph*> graphs;
    std::vector<double> x_values(n_bins);

    for (int i = 0; i < n_methods; ++i) {
        std::vector<double> normalized_times(n_bins);
        for (int j = 0; j < n_bins; ++j) {
            if (i == 0) {
                x_values[j] = j + 0.5;
            }
            if (n_events[j] > 0) {
                normalized_times[j] = all_timings[i][j] / n_events[j];
            } else {
                normalized_times[j] = 0; 
            }
            
            if (normalized_times[j] > max_normalized_time) {
                max_normalized_time = normalized_times[j];
            }
            if (normalized_times[j] < min_normalized_time && normalized_times[j] > 0) {
                min_normalized_time = normalized_times[j];
            }
        }
        graphs.push_back(new TGraph(n_bins, x_values.data(), normalized_times.data()));
    }

    TCanvas *c = new TCanvas("c_timing_comparison", "Timing Comparison", 1000, 600);
    c->SetGrid();
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetRightMargin(0.05);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    const char * selTypeName;
    if (selectionVarType == ControlVar::CENTHF){
        selTypeName = "Timing Comparison;Centrality (HFsumET);Time / Event (s)";
    } else if (selectionVarType == ControlVar::MULT)
    {
        selTypeName = "Timing Comparison;Multiplicity (#Tracks);Time / Event (s)";
    } else {
        selTypeName = "Timing Comparison;Centrality (%);Time / Event (s)";
    }

    TH1F *frame = new TH1F("frame", selTypeName, n_bins, 0, n_bins);
    frame->SetStats(0);
    
    double y_min_range = min_normalized_time * 0.9;
    double y_max_range = max_normalized_time * 1.1;
    if (y_min_range == y_max_range) {
        y_min_range *= 0.9;
        y_max_range *= 1.1;
    }
    if (y_min_range == 0 && y_max_range == 0){
        y_max_range = 1.0;
    }
    
    frame->GetYaxis()->SetRangeUser(y_min_range, y_max_range); 
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->SetRangeUser(0, n_bins);

    for (int j = 0; j < n_bins; ++j) {
        frame->GetXaxis()->SetBinLabel(j + 1, centrality_labels[j].c_str());
    }
    frame->GetXaxis()->LabelsOption("C"); 

    TLegend *legend = new TLegend(0.7, 0.75, 0.93, 0.9); 
    legend->SetFillStyle(0); 
    legend->SetBorderSize(0); 

    std::vector<int> colors(n_methods);
    for (int i = 0; i < n_methods; ++i) {
        double norm = (n_methods == 1) ? 0.5 : static_cast<double>(i) / (n_methods - 1);
        colors[i] = TColor::GetColorPalette(static_cast<int>(norm * (gStyle->GetNumberOfColors() - 1)));
    }

    for (int i = 0; i < n_methods; ++i) {
        graphs[i]->SetMarkerSize(1.2);
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        legend->AddEntry(graphs[i], method_names[i].c_str(), "lp");
        
    }

    frame->Draw();
    for (int i = 0; i < n_methods; ++i) {
        graphs[i]->Draw("LP SAME");
    }
    legend->Draw();

    TCanvas* canvases[1] = {c};
    const char *ipath = "./imgs/test/timing_comparison/";
    save_canvas_images(canvases, 1, ipath, output_filename_base, "png");
    save_canvas_images(canvases, 1, ipath, output_filename_base, "pdf");

    delete c;
    delete frame;
    delete legend;
    for (auto g : graphs) {
        delete g;
    }
}


void runTimingAnalysis(
    ControlVar methodSelectionVar,
    const std::vector<std::string>& method_names,
    const std::vector<std::string>& method_prefixes,
    const char* benchmarkPath,
    const std::string& time_label_to_plot,
    const std::vector<double>& centralityBins,
    const char* output_filename_base
) {
    const char* selVarName = getSelVarName(methodSelectionVar);

    // --- Data Loading ---
    std::cout << "Loading benchmark data from: " << benchmarkPath << std::endl;
    std::cout << "Plotting metric: " << time_label_to_plot << std::endl;

    std::vector<std::vector<double>> all_timings(method_names.size());
    std::vector<double> n_events;
    std::vector<std::string> centrality_labels;

    // Loop over all defined centrality bins
    for (size_t i = 0; i < centralityBins.size() - 1; ++i) {
        double bin_low = centralityBins[i];
        double bin_high = centralityBins[i+1];
        
        // Create the label for the X-axis
        TString label = TString::Format("%.0f-%.0f", bin_low, bin_high);
        if (methodSelectionVar == ControlVar::CENTHF) label += " HF";
        centrality_labels.push_back(label.Data());

        std::cout << "\nProcessing bin: " << label.Data() << std::endl;
        bool n_events_fetched = false;

        // Loop over all methods
        for (size_t j = 0; j < method_names.size(); ++j) {
            TString pattern = TString::Format("%s/%s_%s_%f-%f_benchmark-*.txt",
                                              benchmarkPath,
                                              method_prefixes[j].c_str(),
                                              selVarName,
                                              bin_low,
                                              bin_high);

            TString file = findFile(pattern, "");

            if (file.IsNull()) {
                std::cerr << "-> ERROR: Benchmark file not found for method '" << method_names[j] << "'" << std::endl;
                std::cerr << "    Pattern: " << pattern << std::endl;
                all_timings[j].push_back(0.0); 
                continue;
            }

            std::cout << "-> Found '" << method_names[j] << "' file: " << file << std::endl;
            
            auto data = parseBenchmarkData(file.Data());

            // Extract the desired timing data
            if (data.find(time_label_to_plot) != data.end()) {
                all_timings[j].push_back(data[time_label_to_plot]);
            } else {
                std::cerr << "    Warning: Metric '" << time_label_to_plot << "' not found in " << file << ". Using 0.0." << std::endl;
                all_timings[j].push_back(0.0);
            }

            // Extract the event count 
            if (!n_events_fetched) {
                if (data.find("Number of Processed Events") != data.end()) {
                    n_events.push_back(data["Number of Processed Events"]);
                    n_events_fetched = true;
                } else {
                     std::cerr << "    Warning: 'Number of Processed Events' not found in " << file << ". Using 0." << std::endl;
                     n_events.push_back(0); 
                     n_events_fetched = true; 
                }
            }
        } 
    } 

    if (n_events.empty()) {
        std::cerr << "\nError: No data was successfully loaded. Aborting plot." << std::endl;
        return; 
    }

    std::cout << "\nAll data loaded. Generating plot..." << std::endl;
    plot_timing_comparison(method_names, all_timings, n_events, centrality_labels, output_filename_base, methodSelectionVar);
    
    std::cout << "Plotting finished." << std::endl;
}


int main() {
    // To compile use:
    // g++ -std=c++17 -pthread compareTiming.cpp -o compareTiming `root-config --cflags --libs`
    // To run use:
    // ./compareTiming
    ControlVar methodSelectionVar = ControlVar::CENTHF;

    std::vector<std::string> method_names = {"Serialized", "Parallel"};
    std::vector<std::string> method_prefixes = {"sig_double_loop", "sig_double_loop_parallel"};
    
    const char* benchmarkPath = "benchmarks";
    const std::string time_label_to_plot = "Total Time";
    const char* output_filename_base = "timing_comparison_TGraph";

    std::vector<double> centralityBins = {3200, 3300}; 
    
    runTimingAnalysis(
        methodSelectionVar,
        method_names,
        method_prefixes,
        benchmarkPath,
        time_label_to_plot,
        centralityBins,
        output_filename_base
    );
    
    return 0;
}