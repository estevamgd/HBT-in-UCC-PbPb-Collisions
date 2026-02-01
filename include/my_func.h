#ifndef MY_FUNC_H
#define MY_FUNC_H

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
#include "AnalysisLog.h"

using FourVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;

enum class ControlVar { CENT = 0, MULT = 1, CENTHF = 2 };

enum class qMode { QINV = 0, QLCMS = 1 };

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

const char* getSelVarName(ControlVar varType) {
    if (varType == ControlVar::CENT) return "CENT";
    if (varType == ControlVar::MULT) return "MULT";
    if (varType == ControlVar::CENTHF) return "CENTHF";
    return "UNKNOWN";
}

enum LoopMode { SINGLE = 0, BOTH = 1, DOUBLE = 2 }; // Define an enum for modes


static const Int_t colors[9] = {
    kRed-7,
    kYellow-7,
    kGreen-7,
    kBlue-7,
    kPink+10,
    kOrange+10,
    kSpring+10,
    kCyan-7,
    kMagenta-4
};

static const Int_t lineStyles[4] = {
    1, // solid
    2, // dashed
    3, // dotted
    4  // dash-dotted
};

struct EventData {
    std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> tracks;
    std::vector<Int_t> charges;
};

void printProgressBar(Long64_t current, Long64_t total, int width = 50) {
    if (total <= 0) return;
    double progress = (double)current / total;
    int pos = width * progress;

    std::cout << "\r[" ;
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")" << std::flush;
    if (current == total) std::cout << std::endl;
}

TString findFile(const TString& pattern, const TString& version_tag = "") {
    TString command;
    if (version_tag.IsNull() || version_tag.IsWhitespace()) {
        // Default behavior: find the newest file by sorting alphabetically in reverse (works with YYYY-MM-DD filenames).
        command = TString::Format("ls -1r %s 2>/dev/null | head -n 1", pattern.Data());
    } else {
        // Specific version: find the file containing the version tag, return the newest match if multiple.
        command = TString::Format("ls -1r %s 2>/dev/null | grep %s | head -n 1", pattern.Data(), version_tag.Data());
    }
    
    TString result = TString(gSystem->GetFromPipe(command)).Strip(TString::kBoth, '\n');
    return result;
}

TH1D* getHistogram(TString fileSearchPattern, const char* histName) {
    TString dataFile = findFile(fileSearchPattern);

    if (dataFile.IsNull()) {
        std::cerr << "Error: No data file found matching pattern: "
                  << fileSearchPattern << std::endl;
        return nullptr;
    }

    std::cout << "Found data file: " << dataFile << std::endl;

    TFile* file = TFile::Open(dataFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file: " << dataFile << std::endl;
        return nullptr;
    }
    AnalysisLog::instance().addInputFile(dataFile.Data());

    TH1D* hist = dynamic_cast<TH1D*>(file->Get(histName));
    if (!hist) {
        std::cerr << "Error: Histogram '" << histName
                  << "' not found in file." << std::endl;
        file->Close();
        delete file;
        return nullptr;
    }

    TH1D* histClone = dynamic_cast<TH1D*>(
        hist->Clone(Form("%s_clone", histName))
    );

    histClone->SetDirectory(nullptr);
    
    file->Close();
    delete file;

    return histClone;
}

TH2D* getHistogram2d(TString fileSearchPattern, const char* histName) {
    TString dataFile = findFile(fileSearchPattern);

    if (dataFile.IsNull()) {
        std::cerr << "Error: No data file found matching pattern: "
                  << fileSearchPattern << std::endl;
        return nullptr;
    }

    std::cout << "Found data file: " << dataFile << std::endl;

    TFile* file = TFile::Open(dataFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file: " << dataFile << std::endl;
        return nullptr;
    }
    AnalysisLog::instance().addInputFile(dataFile.Data());

    TH2D* hist = dynamic_cast<TH2D*>(file->Get(histName));
    if (!hist) {
        std::cerr << "Error: Histogram '" << histName
                  << "' not found in file." << std::endl;
        file->Close();
        delete file;
        return nullptr;
    }

    TH2D* histClone = dynamic_cast<TH2D*>(
        hist->Clone(Form("%s_clone", histName))
    );

    histClone->SetDirectory(nullptr);
    
    file->Close();
    delete file;

    return histClone;
}

std::string getTimestamp() {
    time_t ttime = time(NULL);
    struct tm date = *localtime(&ttime);
    char buffer[80];
    strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", &date);
    return std::string(buffer);
}

// Input file and tree names get tree and file
bool getFileTree(const char* file_name, const char* tree_name, TFile *&fr, TTree *&t) {
    // Open the input file
    fr = new TFile(file_name, "READ");
    if (!fr || fr->IsZombie()) {
        std::cerr << "Error: File could not be opened!" << std::endl;
        return false;
    }
    AnalysisLog::instance().addInputFile(file_name);

    // Get the tree from the specified path
    fr->GetObject(tree_name, t);
    if (!t) {
        std::cerr << "Error: Tree '" << tree_name << "' not found!" << std::endl;
        fr->Close();
        delete fr;
        return false;
    }

    return true;
}

// Bin size calculator
double numBins(double interval, double length ,double scale) {
    double n = (interval / length) * scale;
    return n;
}

/* Histogram layout

Colors: 
kWhite  = 0, kBlack     = 1, kGray      = 920, kRed     = 632, kGreen   = 416, 
kBlue   = 600, kYellow  = 400, kMagenta = 616, kCyan    = 432, kOrange  = 800, 
kSpring = 820, kTeal    = 840, kAzure   = 860, kViolet  = 880, kPink    = 900

Fill Style:
kFDotted1  = 3001, kFDotted2    = 3002, kFDotted3  = 3003,
kFHatched1 = 3004, kHatched2    = 3005, kFHatched3 = 3006,
kFHatched4 = 3007, kFWicker     = 3008, kFScales   = 3009,
kFBricks   = 3010, kFSnowflakes = 3011, kFCircles  = 3012,
kFTiles    = 3013, kFMondrian   = 3014, kFDiamonds = 3015,
kFWaves1   = 3016, kFDashed1    = 3017, kFDashed2  = 3018,
kFAlhambra = 3019, kFWaves2     = 3020, kFStars1   = 3021,
kFStars2   = 3022, kFPyramids   = 3023, kFFrieze   = 3024,
kFMetopes  = 3025, kFEmpty      = 0   , kFSolid    = 1

Line Style:
kSolid = 1, kDashed = 2, kDotted = 3, kDashDotted = 4
*/
TH1D* cHist(const char* name, const char* xAxisTitle, const char* yAxisTitle, 
                        double xMin, double xMax, 
                        double ninterval = 275., double nlength = 0.02, double nscale = 1., 
                        int fill_color = 920, double fill_alpha = 1, int fill_style = 1001, 
                        int line_color = 1, double line_alpha = 1, int line_style = 1001, 
                        int line_width = 1, double label_size = 0.04) {
    // Create the histogram
    double scale = numBins(ninterval, nlength, nscale);
    TH1D* h = new TH1D(name, "", scale, xMin, xMax);

    // Set visual properties
    h->SetFillColorAlpha(fill_color, fill_alpha); 
    h->SetFillStyle(fill_style); 
    h->SetLineColorAlpha(line_color, line_alpha);
    h->SetLineStyle(line_style);
    h->SetLineWidth(line_width); 

    // Set titles and axis labels
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xAxisTitle);
    h->GetYaxis()->SetTitle(yAxisTitle);
    
    // Set label sizes
    h->GetXaxis()->SetLabelSize(label_size);
    h->GetYaxis()->SetLabelSize(label_size);
    
    return h;
}

TH1F* cHistFloat(const char* name, const char* xAxisTitle, const char* yAxisTitle, 
    double xMin, double xMax, 
    double ninterval = 275., double nlength = 0.02, double nscale = 1., 
    int fill_color = 920, double fill_alpha = 1, int fill_style = 1001, 
    int line_color = 1, double line_alpha = 1, int line_style = 1001, 
    int line_width = 1, double label_size = 0.04) {
// Create the histogram
double scale = numBins(ninterval, nlength, nscale);
TH1F* h = new TH1F(name, "", scale, xMin, xMax);

// Set visual properties
h->SetFillColorAlpha(fill_color, fill_alpha); 
h->SetFillStyle(fill_style); 
h->SetLineColorAlpha(line_color, line_alpha);
h->SetLineStyle(line_style);
h->SetLineWidth(line_width); 

// Set titles and axis labels
h->SetTitle("");
h->GetXaxis()->SetTitle(xAxisTitle);
h->GetYaxis()->SetTitle(yAxisTitle);

// Set label sizes
h->GetXaxis()->SetLabelSize(label_size);
h->GetYaxis()->SetLabelSize(label_size);

return h;
}

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
            AnalysisLog::instance().addSavedObject(
                    LogEntryType::SavedCanvas,
                    "Canvas image saved",
                    path,
                    canvas_name
            );
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
        AnalysisLog::instance().addSavedObject(
                LogEntryType::SavedHistogram,
                "Histogram ROOT file saved",
                path,
                root_name
        );
        
        for (int i = 0; i < numHistograms; i++) {
            histograms[i]->Write();
        }
    }
    
}

void save_histograms2d(TH2D *histograms[], int numHistograms, const char *path, const char *prefix,
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
        AnalysisLog::instance().addSavedObject(
                LogEntryType::SavedHistogram,
                "Histogram ROOT file saved",
                path,
                root_name
        );

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
        AnalysisLog::instance().addSavedObject(
                LogEntryType::SavedHistogram,
                "Histogram ROOT file saved",
                path,
                root_name
        );

        std::cout << "Saving histograms to file: " << root_name << std::endl;
        for (int i = 0; i < numHistograms; i++) {
            histograms[i]->Write();
        }
    }
    
}

void save_histograms(TH2D *histograms[], int numHistograms, const char *path, const char *prefix) {
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

        TFile saveFile(root_name, "NEW");
        AnalysisLog::instance().addSavedObject(
                LogEntryType::SavedHistogram,
                "Histogram ROOT file saved",
                path,
                root_name
        );

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

    AnalysisLog::instance().addSavedObject(
            LogEntryType::SavedBenchmark,
            "Benchmark report saved",
            path,
            fileName
    );

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
    
    AnalysisLog::instance().addSavedObject(
            LogEntryType::SavedBenchmark,
            "Benchmark report saved",
            path,
            fileName
    );

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
    
    AnalysisLog::instance().addSavedObject(
            LogEntryType::SavedBenchmark,
            "Benchmark report saved",
            path,
            fileName
    );

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
    TFile *fr = nullptr)
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

    // Close and delete the TFile
    if (fr != nullptr) {
    fr->Close();
    delete fr;
    }
}

// Deletes canvases, histograms, legends, and closes the TFile
void close_program(TCanvas *canvases[] = nullptr, int numCanvases = 0,
    TH2D *histograms[] = nullptr, int numHistograms = 0,
    TLegend *legends[] = nullptr, int numLegends = 0,
    TFile *fr = nullptr)
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

    // Close and delete the TFile
    if (fr != nullptr) {
    fr->Close();
    delete fr;
    }
}

void drawCMSHeaders(
    const char* leftText  = "#bf{CMS} #it{Work in Progress}",
    const char* rightText = "",
    double rightTextSizeFactor = 1.0
) {
    if (!gPad) return;

    // ---- HARD RESET OF PAD STATE ----
    gPad->SetTitle("");
    gPad->SetName("");
    gPad->SetGrid(0, 0);
    gPad->SetTicks(1, 1);
    gPad->Update();

    const double yPos = 1.0 - gPad->GetTopMargin() + 0.01;

    // ===== Left CMS label =====
    TLatex latexLeft;
    latexLeft.SetNDC();
    latexLeft.SetTextFont(42);
    latexLeft.SetTextSize(0.035);
    latexLeft.SetTextAlign(11);

    latexLeft.DrawLatex(
        gPad->GetLeftMargin(),
        yPos,
        leftText
    );

    // ===== Right header =====
    if (rightText && rightText[0] != '\0') {

        TLatex latexRight;
        latexRight.SetNDC();
        latexRight.SetTextFont(42);
        latexRight.SetTextAlign(31);

        const double xRight = 1.0 - gPad->GetRightMargin();
        const double xLeft  = gPad->GetLeftMargin() + 0.25;
        const double maxWidth = xRight - xLeft;

        double textSize = 0.04 * rightTextSizeFactor;
        latexRight.SetTextSize(textSize);

        latexRight.SetTitle(rightText);
        double textWidth = latexRight.GetXsize();

        while (textWidth > maxWidth && textSize > 0.015) {
            textSize *= 0.95;
            latexRight.SetTextSize(textSize);
            latexRight.SetTitle(rightText);
            textWidth = latexRight.GetXsize();
        }

        latexRight.DrawLatex(xRight, yPos, rightText);
    }

    // ---- FINAL FLUSH ----
    gPad->Modified();
    gPad->Update();
}



void no_statbox(TH1D *histograms[], int numHistograms) {
    for (int i = 0; i < numHistograms; i++) {
        histograms[i]->SetStats(0);
    }
}

/*void setCentralities(std::vector<std::vector<int>>& centralitiesSelected, const std::vector<int>& centralityX) {
    std::vector<int> centrailityAux;
    centrailityAux.reserve(2);

    const size_t numPairs = (centralityX.size() + 1) / 2;
    centralitiesSelected.reserve(centralitiesSelected.size() + numPairs);

    for (size_t i = 0; i < centralityX.size(); ++i) {
        centrailityAux.push_back(centralityX[i]);

        if (i % 2 == 1) {
            centralitiesSelected.push_back(centrailityAux);
            centrailityAux.clear();
        }
    }
}
*/


int setLoopMode(LoopMode mode) {
    int xmode = 2;

    switch (mode) {
        case SINGLE:
            xmode = 0;
            break;
        case BOTH:
            xmode = 1;
            break;
        case DOUBLE:
            xmode = 2;
            break;
        default:
            std::cout << "Choose between: SINGLE, BOTH, or DOUBLE." << std::endl;
    }
    return xmode;
}

//enum class ControlVar { CENT = 0, MULT = 1, VERTEX = 2 }; // Define an enum for types of control var
//constexpr int nControlVar = 3;
//
//struct TrackInfo {
//    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> p4;
//    Int_t charge;    
//    Float_t weight;  
//};



#endif 
