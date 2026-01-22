#ifndef GETFUNCS_H
#define GETFUNCS_H


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


enum class ControlVar { CENT = 0, MULT = 1, CENTHF = 2 };

const char* getSelVarName(ControlVar varType) {
    if (varType == ControlVar::CENT) return "CENT";
    if (varType == ControlVar::MULT) return "MULT";
    if (varType == ControlVar::CENTHF) return "CENTHF";
    return "UNKNOWN";
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


#endif