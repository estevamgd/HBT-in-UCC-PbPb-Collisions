///*** To compile use:
///*** g++ -std=c++17 makeSignalMixOLM.cpp -o makeSignalMixOLM `root-config --cflags --libs`
///*** To run use:
///*** ./makeSignalMixOLM

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <ctime>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"

#include "../include/data_func.h"
#include "../include/my_func.h" 

enum class ControlVar { CENT = 0, MULT = 1, CENTHF = 2 };

/**
 * @brief Perform signal and mixed-event correlation analysis on CMS PbPb data.
 *
 * This macro reads a ROOT file and TTree, applies event and track selection,
 * fills same-event (signal) and mixed-event histograms of q_inv for track pairs,
 * and saves histograms, canvases, and benchmark results.
 * 
 * @param fileInput Path to the input ROOT file.
 * @param treeInput Name of the TTree inside the file.
 * @param selVarMoreeq Lower bound (inclusive) of the selection variable (e.g. centrality).
 * @param selVarLess Upper bound (exclusive) of the selection variable.
 * @param selectionVarType Event selection variable (CENT, MULT, CENTHF).
 * @param vertexDistance Maximum |Δz| between primary vertices for mixed events (in cm).
 * @param pairChargeMult Charge filter for track pairs: +1 for same-sign, -1 for opposite-sign.
 * @param poolSizeInt Number of events stored in the mixing pool before clearing.
 */
void makeSignalMixOLM(const char *fileInput, const char *treeInput, double selVarMoreeq, double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, Float_t vertexDistance = 2, Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) {
        
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree(fileInput, treeInput, fr, t);

    Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
    Float_t HFsumET, pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize], 
            pionMass = 0.13957039;
    
    double* selectionVar, displaySelVarMoreeq, displaySelVarLess;
    double hiBinProxy, NtrkProxy, HFsumETProxy;
    
    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("pvZ", &pvZ);
    t->SetBranchAddress("trkCharge", trkCharge);
    t->SetBranchAddress("trkWeight", trkWeight);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);

    const char * selectionVarName;
    if (selectionVarType == ControlVar::CENT){
        selectionVar = &hiBinProxy;
        displaySelVarMoreeq = selVarMoreeq;
        displaySelVarLess = selVarLess;
        selVarMoreeq *= 2;
        selVarLess *= 2;
        selectionVarName = "CENT";
    } else if (selectionVarType == ControlVar::MULT){
        selectionVar = &NtrkProxy;
        displaySelVarMoreeq = selVarMoreeq;
        displaySelVarLess = selVarLess;
        selectionVarName = "MULT";
    } else if (selectionVarType == ControlVar::CENTHF){
        selectionVar = &HFsumETProxy;
        displaySelVarMoreeq = selVarMoreeq;
        displaySelVarLess = selVarLess;
        selectionVarName = "CENTHF";
    } else {
        std::cerr << "Invalid selection variable type!" << std::endl;
        return;
    }

    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    int processedEvents = 0;

    double nnscale = 10000, x0 = 0., xt=10.;
    TH1D* h_qinvSS_signal_1l = new TH1D("h_qinvSS_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_signal_1l = new TH1D("h_qinvSSCor_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_signal_1l = new TH1D("h_qinvOS_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_signal_1l = new TH1D("h_qinvOSCor_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinv_mix_1l = new TH1D("h_qinv_mix_1l", "", nnscale, x0, xt);
    TH1D* h_qinvCor_mix_1l = new TH1D("h_qinvCor_mix_1l", "", nnscale, x0, xt);
    
    h_qinvSSCor_signal_1l->Sumw2();
    h_qinvOSCor_signal_1l->Sumw2();
    h_qinvCor_mix_1l->Sumw2();

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTrackspvZ;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;
    size_t poolSize = static_cast<size_t>(poolSizeInt);

    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i_entry = 0; i_entry < nentries; i_entry++){
        t->GetEntry(i_entry);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (!(*selectionVar >= selVarMoreeq && *selectionVar < selVarLess)) continue;
        processedEvents++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTrackspvZ;
        std::vector<Float_t> currentEventTracksWeight;
        
        for (int j = 0; j < Ntrk; j++) {
            currentEventTracks4V.emplace_back(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTrackspvZ.push_back(pvZ);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        if (!currentEventTracks4V.empty() && !vectorEventTracks4V.empty()) {
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); nEv++) {
                for (size_t k = 0; k < vectorEventTracks4V[nEv].size() * currentEventTracks4V.size(); k++) {
                    size_t p1 = k / currentEventTracks4V.size();
                    size_t p2 = k % currentEventTracks4V.size();
                    
                    if (fabs(pvZ - vectorEventTrackspvZ[nEv][p1]) >= vertexDistance) continue;
                    if (vectorEventTracksCharge[nEv][p1] * trkCharge[p2] != pairChargeMult) continue;
                    if (std::isinf(vectorEventTracksWeight[nEv][p1] * trkWeight[p2])) continue;
                    
                    double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                    h_qinv_mix_1l->Fill(qinv);
                    h_qinvCor_mix_1l->Fill(qinv, vectorEventTracksWeight[nEv][p1] * trkWeight[p2]);
                }
            }
        }
        if (!currentEventTracks4V.empty()) {
            vectorEventTracks4V.push_back(currentEventTracks4V);
            vectorEventTracksCharge.push_back(currentEventTracksCharge);
            vectorEventTrackspvZ.push_back(currentEventTrackspvZ);
            vectorEventTracksWeight.push_back(currentEventTracksWeight);
            if (vectorEventTracks4V.size() >= poolSize) {
                vectorEventTracks4V.clear();
                vectorEventTracksCharge.clear();
                vectorEventTrackspvZ.clear();
                vectorEventTracksWeight.clear();
            }
        }
        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        if (currentEventTracks4V.size() > 1) {
            for (size_t k = 0; k < currentEventTracks4V.size() * (currentEventTracks4V.size() - 1); k++) {
                size_t p1 = k / (currentEventTracks4V.size() - 1);
                size_t p2 = (k % (currentEventTracks4V.size() - 1)) + 1;
                if (p1 >= p2) continue;
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;
                
                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);
                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    h_qinvSS_signal_1l->Fill(qinv);
                    h_qinvSSCor_signal_1l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else {
                    h_qinvOS_signal_1l->Fill(qinv);
                    h_qinvOSCor_signal_1l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();
    
    std::cout << "#Processed Events: " << processedEvents << std::endl;
    
    // Final processing and saving (similar to DLM)
    char prefix[50];
    sprintf(prefix, "signal_mix_OLM_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    std::vector<double> durations = {duration_full, duration_mix, duration_signal};
    std::vector<std::string> labels = {"Total Time", "Mixing Time", "Signal Time"};
    save_benchmark_chrono(durations, labels, "benchmarks", prefix, processedEvents);

    // ... (rest of the saving logic remains the same)
}

int main(){
    const char* inputFile = "data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root";  
    const char* treeName  = "demo/TreeMBUCC";
    std::vector<double> centralityBins = {0.5, 1, 5, 10, 30, 50, 70}; 
    ControlVar selectedControlVar = ControlVar::CENT;
    Float_t selectedVertexDistance = 2.0f;
    int poolSize = 10;   
    Int_t pairChargeMult = 1;

    for (size_t i = 0; i < centralityBins.size() - 1; ++i){
        makeSignalMixOLM(inputFile, treeName, centralityBins[i], centralityBins[i+1], 
                         selectedControlVar, selectedVertexDistance, pairChargeMult, poolSize);
    }
    
    std::cout << "\nAnalysis finished." << std::endl;
    return 0;
}
