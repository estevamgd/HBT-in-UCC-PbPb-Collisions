#ifndef SAVEFUNCS_H
#define SAVEFUNCS_H

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <TStyle.h>
#include <fstream>
#include <ctime>
#include <sstream>
#include <filesystem>
#include "my_func.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "TMath.h"
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>
#include <chrono>
#include <string>
#include "data_func.h"
#include <TSystem.h>
#include <TLatex.h>
#include "my_func.h"


struct BenchmarkCentralityEntry {
    std::string label;          // e.g. "0–5%"
    double nEvents;             // number of processed events
    std::vector<double> times;  // one per method (seconds)
};

struct BenchmarkMetadata {
    std::string controlVariable;            // CENTRALITY, MULTIPLICITY, CENTHF, etc.
    std::vector<std::string> inputFiles;    // files being compared
    std::vector<std::string> methodNames;   // e.g. Single Loop, Double Loop
};

// Saving the canvases
void save_canvas_images(TCanvas *canvases[], int numCanvases, const char *path, const char *prefix, const char *file_type) {
    // Cria diretório se necessário
    std::filesystem::create_directories(path);
    
    if (numCanvases > 0) {
        // Get the current time
        time_t ttime = time(NULL);
        struct tm date = *localtime(&ttime);
        
        // Loop through all canvases and save them
        for (int i = 0; i < numCanvases; i++) {
            char canvas_name[200];  // Adjust size to hold path, prefix, and time info
            // Generate a unique file name for each canvas
            sprintf(canvas_name, "%s%s-%d-%02d-%02d-%02d-%02d-%02d-%02d.%s", 
                    path, prefix, 
                    date.tm_year + 1900, 
                    date.tm_mon + 1, 
                    date.tm_mday, 
                    date.tm_hour, 
                    date.tm_min, 
                    date.tm_sec, 
                    i + 1,    // Canvas index
                    file_type);  // File type
        
            // Save the canvas as an image with the specified file type
            canvases[i]->Print(canvas_name);
        }
    }
    
}

// Saving histograms
void save_histograms(TH1D *histograms[], int numHistograms, const char *path, const char *prefix,
                    int centmult_min, int centmult_max) {
    // Cria diretório se necessário
    std::filesystem::create_directories(path);
    
    if (numHistograms > 0) {
        const char *file_type = "root";
        time_t ttime = time(NULL);
        struct tm date = *localtime(&ttime);
        char root_name[200];  

        sprintf(root_name, "%s%s_cent%dto%d-%d-%02d-%02d-%02d-%02d-%02d.%s", 
            path, prefix, 
            centmult_min, centmult_max,
            date.tm_year + 1900, 
            date.tm_mon + 1, 
            date.tm_mday, 
            date.tm_hour, 
            date.tm_min, 
            date.tm_sec, 
            file_type);

        TFile saveFile(root_name, "NEw");

        for (int i = 0; i < numHistograms; i++) {
            histograms[i]->Write();
        }
    }
    
}

void save_histograms(TH1D *histograms[], int numHistograms, const char *path, const char *prefix) {
    // Cria diretório se necessário
    std::filesystem::create_directories(path);
    
    if (numHistograms > 0) {
        const char *file_type = "root";
        time_t ttime = time(NULL);
        struct tm date = *localtime(&ttime);
        char root_name[200];  

        sprintf(root_name, "%s%s-%d-%02d-%02d-%02d-%02d-%02d.%s", 
            path, prefix, 
            date.tm_year + 1900, 
            date.tm_mon + 1, 
            date.tm_mday, 
            date.tm_hour, 
            date.tm_min, 
            date.tm_sec, 
            file_type);

        TFile saveFile(root_name, "NEw");
        std::cout << "Saving histograms to file: " << root_name << std::endl;
        for (int i = 0; i < numHistograms; i++) {
            histograms[i]->Write();
        }
    }
    
}

// Save information
void save_benchmark(TStopwatch* stopWatches[], int numSW, const char* path, const char* prefix, int nProcessedEvents = -1) {
    // Cria diretório se necessário
    std::filesystem::create_directories(path);

    // Pega data/hora atual
    time_t ttime = time(NULL);
    struct tm date = *localtime(&ttime);

    // Monta nome do arquivo
    char fileName[256];
    sprintf(fileName, "%s/%s_benchmark-%d-%02d-%02d-%02d-%02d-%02d.txt", 
            path, prefix, 
            date.tm_year + 1900, 
            date.tm_mon + 1, 
            date.tm_mday, 
            date.tm_hour, 
            date.tm_min, 
            date.tm_sec);

    // Abre arquivo
    std::ofstream outFile(fileName);
    if (!outFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo para salvar benchmark: " << fileName << std::endl;
        return;
    }

    outFile << "===== Benchmark - "
            << (1900 + date.tm_year) << "/"
            << (date.tm_mon + 1) << "/"
            << date.tm_mday << " "
            << date.tm_hour << ":" << date.tm_min << ":" << date.tm_sec
            << " =====\n";

    outFile << std::fixed << std::setprecision(3);

    for (int i = 0; i < numSW; ++i) {
        outFile << "Stopwatch " << i << ": "
                << stopWatches[i]->RealTime() << " s (real), "
                << stopWatches[i]->CpuTime()  << " s (CPU)" << "\n";
    }

    outFile << "=========================================\n";
    outFile << "Number of Processed Events: " << nProcessedEvents << "\n";
    outFile << "=========================================\n";
    outFile.close();

    std::cout << "Benchmark salvo em: " << fileName << std::endl;
}

void save_benchmark_chrono(const std::vector<double>& durations, const std::vector<std::string>& labels, const char* path, const char* prefix, 
    int nProcessedEventsSig = -1, int nProcessedEventsMix = -1) {
    std::filesystem::create_directories(path);

    time_t ttime = time(NULL);
    struct tm date = *localtime(&ttime);

    char fileName[256];
    sprintf(fileName, "%s/%s_benchmark-%d-%02d-%02d-%02d-%02d-%02d.txt",
            path, prefix,
            date.tm_year + 1900,
            date.tm_mon + 1,
            date.tm_mday,
            date.tm_hour,
            date.tm_min,
            date.tm_sec);

    std::ofstream outFile(fileName);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file to save benchmark: " << fileName << std::endl;
        return;
    }

    outFile << "===== Benchmark - "
            << (1900 + date.tm_year) << "/"
            << (date.tm_mon + 1) << "/"
            << date.tm_mday << " "
            << date.tm_hour << ":" << date.tm_min << ":" << date.tm_sec
            << " =====\n";

    outFile << std::fixed << std::setprecision(6);

    if (durations.size() != labels.size()) {
         outFile << "Error: Mismatch between number of durations and labels.\n";
    } else {
        for (size_t i = 0; i < durations.size(); ++i) {
            outFile << labels[i] << ": "
                    << durations[i] << " s\n";
        }
    }

    outFile << "=========================================\n";
    outFile << "Number of Processed Events(Sig): " << nProcessedEventsSig << "\n";
    outFile << "=========================================\n";
    outFile << "=========================================\n";
    outFile << "Number of Processed Events(Mix): " << nProcessedEventsMix << "\n";
    outFile << "=========================================\n";
    outFile.close();

    std::cout << "Benchmark saved to: " << fileName << std::endl;
}

void save_benchmark_validation(
    const BenchmarkMetadata& meta,
    const std::vector<BenchmarkCentralityEntry>& entries,
    const char* path,
    const char* prefix,
    double normalization_factor = 1.0 // 1.0 → seconds, 60 → minutes, etc.
) {
    std::filesystem::create_directories(path);

    time_t ttime = time(nullptr);
    tm date = *localtime(&ttime);

    char fileName[512];
    sprintf(fileName,
            "%s/%s_benchmark-%04d-%02d-%02d-%02d-%02d-%02d.txt",
            path, prefix,
            date.tm_year + 1900,
            date.tm_mon + 1,
            date.tm_mday,
            date.tm_hour,
            date.tm_min,
            date.tm_sec);

    std::ofstream out(fileName);
    if (!out.is_open()) {
        std::cerr << "Error opening benchmark file: " << fileName << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(6);

    /* ===================== HEADER ===================== */
    out << "===== BENCHMARK REPORT =====\n";
    out << "Date: "
        << (1900 + date.tm_year) << "-"
        << (date.tm_mon + 1) << "-"
        << date.tm_mday << " "
        << date.tm_hour << ":"
        << date.tm_min << ":"
        << date.tm_sec << "\n\n";

    out << "Control variable: " << meta.controlVariable << "\n\n";

    /* ===================== FILES ===================== */
    out << "Compared input files:\n";
    for (const auto& f : meta.inputFiles)
        out << "  - " << f << "\n";
    out << "\n";

    /* ===================== METHODS ===================== */
    out << "Compared methods:\n";
    for (size_t i = 0; i < meta.methodNames.size(); ++i)
        out << "  [" << i << "] " << meta.methodNames[i] << "\n";
    out << "\n";

    /* ===================== RESULTS ===================== */
    out << "Per-bin timing results:\n";
    out << "------------------------------------------------------------\n";

    for (const auto& entry : entries) {
        out << "Bin: " << entry.label << "\n";
        out << "Processed events: " << entry.nEvents << "\n";

        for (size_t m = 0; m < entry.times.size(); ++m) {
            double t_total = entry.times[m] / normalization_factor;
            double t_per_evt = t_total / entry.nEvents;

            out << "  Method: " << meta.methodNames[m] << "\n";
            out << "    Total time      : " << t_total << " s\n";
            out << "    Time per event  : " << t_per_evt << " s/event\n";
        }

        out << "------------------------------------------------------------\n";
    }

    out << "===== END OF BENCHMARK =====\n";
    out.close();

    std::cout << "Benchmark saved to: " << fileName << std::endl;
}

// Deletes canvases, histograms, legends, and closes the TFile
void close_program(TCanvas *canvases[] = nullptr, int numCanvases = 0,
    TH1D *histograms[] = nullptr, int numHistograms = 0,
    TLegend *legends[] = nullptr, int numLegends = 0,
    TFile *fr = nullptr,
    TH2D *histograms2d[] = nullptr, int numHistograms2d = 0)
{
    // Delete Canvases
    if (canvases != nullptr && numCanvases > 0) {
        for (int i = 0; i < numCanvases; i++) {
            delete canvases[i];
    }
    }

    // Delete 1D Histograms
    if (histograms != nullptr && numHistograms > 0) {
        for (int i = 0; i < numHistograms; i++) {
            delete histograms[i];
    }
    }

    // --- NEW: Delete Legends ---
    if (legends != nullptr && numLegends > 0) {
        for (int i = 0; i < numLegends; i++) {
            delete legends[i];
    }
    }

    // --- CORRECTED: Delete 2D Histograms (was in a faulty else block) ---
    if (histograms2d != nullptr && numHistograms2d > 0) {
        for (int i = 0; i < numHistograms2d; i++) {
            delete histograms2d[i];
    }
    }

    // Close and delete the TFile
    if (fr != nullptr) {
    fr->Close();
    delete fr;
    }
}


#endif