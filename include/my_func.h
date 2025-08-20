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
#include <iomanip>
#include <ctime>
#include <sstream>
#include <filesystem>
#include "my_func.h"


// Input file and tree names get tree and file
bool getFileTree(const char* file_name, const char* tree_name, TFile *&fr, TTree *&t) {
    // Open the input file
    fr = new TFile(file_name, "READ");
    if (!fr || fr->IsZombie()) {
        std::cerr << "Error: File could not be opened!" << std::endl;
        return false;
    }

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

// Delete canvases and histograms
void close_program(TCanvas *canvases[], int numCanvases, TH1D *histograms[], int numHistograms, TFile *fr,
                    TH2D *histograms2d[] = nullptr) {
    if (numCanvases > 0) {
        for (int i = 0; i < numCanvases; i++) {
            delete canvases[i];
        }
    }
    if (histograms != nullptr) {
        if (numHistograms > 0) {
            for (int i = 0; i < numHistograms; i++) {
                delete histograms[i];
            }
        }
    }else {
        if (numHistograms > 0) {
            for (int i = 0; i < numHistograms; i++) {
                delete histograms2d[i];
            }
        }
    }
    

    fr->Close();
    delete fr;
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

enum LoopMode { SINGLE = 0, BOTH = 1, DOUBLE = 2 }; // Define an enum for modes

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
