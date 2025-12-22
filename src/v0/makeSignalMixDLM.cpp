///*** To compile use:
///*** g++ -std=c++17 makeSignalMixDLM.cpp -o makeSignalMixDLM `root-config --cflags --libs`
///*** To run use:
///*** ./makeSignalMixDLM

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

void makeSignalMixDLM(const char *fileInput, const char *treeInput, double selVarMoreeq, double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, Float_t vertexDistance = 2, Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) {
        
    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree(fileInput, treeInput, fr, t);

    // Variables from the tree
    Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
    Float_t HFsumET, pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize], 
            pionMass = 0.13957039; // Pion mass [GeV] from PDG
    
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
        selVarMoreeq = selVarMoreeq*2;
        selVarLess = selVarLess*2;
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

    // === HISTOGRAMS ===
    double nnscale = 10000, x0 = 0., xt=10.;

    TH1D* h_qinvSS_signal_2l = new TH1D("h_qinvSS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_signal_2l = new TH1D("h_qinvSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_signal_2l = new TH1D("h_qinvOS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_signal_2l = new TH1D("h_qinvOSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinv_mix_2l = new TH1D("h_qinv_mix_2l", "", nnscale, x0, xt);
    TH1D* h_qinvCor_mix_2l = new TH1D("h_qinvCor_mix_2l", "", nnscale, x0, xt);
    
    h_qinvSSCor_signal_2l->Sumw2();
    h_qinvOSCor_signal_2l->Sumw2();
    h_qinvCor_mix_2l->Sumw2();
    
    h_qinvSS_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvOS_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvSSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvOSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinv_mix_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinv_mix_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvCor_mix_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvCor_mix_2l->GetYaxis()->SetTitle("#Pairs");

    // Benchmarking
    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTrackspvZ;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;

    size_t poolSize = static_cast<size_t>(poolSizeInt);

    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy   = hiBin;
        NtrkProxy    = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue; 
        processedEvents++;                                                      
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTrackspvZ;
        std::vector<Float_t> currentEventTracksWeight;
        
        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTrackspvZ.push_back(pvZ);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        if (!currentEventTracks4V.empty()) { 
            if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()) { 
                for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); nEv++) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); p1++) {
                        for (size_t p2 = 0; p2 < currentEventTracks4V.size(); p2++) {
                            if (!(fabs(pvZ - vectorEventTrackspvZ[nEv][p1]) < vertexDistance)) continue;
                            if (!(vectorEventTracksCharge[nEv][p1]*trkCharge[p2] == pairChargeMult)) continue;
                            if (std::isinf(vectorEventTracksWeight[nEv][p1]*trkWeight[p2])) continue;
                            
                            double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                            h_qinv_mix_2l->Fill(qinv);
                            h_qinvCor_mix_2l->Fill(qinv, vectorEventTracksWeight[nEv][p1]*trkWeight[p2]);
                        }
                    }
                }
            }
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
        // === SIGNAL ===
        if (currentEventTracks4V.size() > 1) { 
            for (size_t p1 = 0; p1 < currentEventTracks4V.size(); p1++) {
                for (size_t p2 = p1 + 1; p2 < currentEventTracks4V.size(); p2++) {
                    if (std::isinf(trkWeight[p1]*trkWeight[p2])) continue;
                    
                    double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);
                    
                    if (trkCharge[p1]*trkCharge[p2] > 0){
                        h_qinvSS_signal_2l->Fill(qinv);
                        h_qinvSSCor_signal_2l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                    } else {
                        h_qinvOS_signal_2l->Fill(qinv);
                        h_qinvOSCor_signal_2l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                    }
                }
            }
        }
        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events: " << processedEvents << std::endl;
    
    TH1D *h_qinvDivCor_signal_2l = (TH1D *)h_qinvSSCor_signal_2l->Clone("h_qinvDivCor_signal_2l");
    TH1D *h_qinvDiv_signal_2l = (TH1D *)h_qinvSS_signal_2l->Clone("h_qinvDiv_signal_2l");
    h_qinvDivCor_signal_2l->Divide(h_qinvOSCor_signal_2l);
    h_qinvDiv_signal_2l->Divide(h_qinvOS_signal_2l);
    h_qinvDiv_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDiv_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvDivCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDivCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvDivCor_signal_2l->Sumw2();

    // --- Save output ---
    TCanvas *c_sig_SS = new TCanvas("c_sig_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_OS = new TCanvas("c_sig_OS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Corrected Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Corrected Distributions", 1200, 800);
    TCanvas *c_sig_cor_Div = new TCanvas("c_sig_cor_Div", "Signal Corrected Division Distributions", 1200, 800);
    TCanvas *c_sig_Div = new TCanvas("c_sig_Div", "Signal Division Distributions", 1200, 800);
    TCanvas *c_mix = new TCanvas("c_mix", "Mixing Distributions", 1200, 800);
    TCanvas *c_mix_cor = new TCanvas("c_mix_cor", "Mixing Corrected Distributions", 1200, 800);
    
    TCanvas *canvases[] = {c_sig_SS, c_sig_OS, c_sig_cor_OS, c_sig_cor_SS, c_sig_cor_Div, c_sig_Div, c_mix, c_mix_cor};
    int numCanvases = 8;

    h_qinvSS_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");
    h_qinvOS_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");
    h_qinvSSCor_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Corrected");
    h_qinvOSCor_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Corrected");
    h_qinvDivCor_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division Corrected");
    h_qinvDiv_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division");
    h_qinv_mix_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Mixing");
    h_qinvCor_mix_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Mixing Corrected");
    
    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_SS->AddEntry(h_qinvSS_signal_2l, "Signal SS - Double loop", "l");
    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qinvSS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_OS->AddEntry(h_qinvOS_signal_2l, "Signal OS - Double loop", "l");
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qinvOS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();

    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_cor_SS->AddEntry(h_qinvSSCor_signal_2l, "Signal SS Cor - Double loop", "l");
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();

    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_cor_OS->AddEntry(h_qinvOSCor_signal_2l, "Signal OS Cor - Double loop", "l");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();

    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_cor_Div->AddEntry(h_qinvDivCor_signal_2l, "Signal Div Cor - Double loop", "l");
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qinvDivCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();

    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_Div->AddEntry(h_qinvDiv_signal_2l, "Signal Div - Double loop", "l");
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qinvDiv_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();

    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_mix->AddEntry(h_qinv_mix_2l, "Mix - Double loop", "l");
    c_mix->cd(); gStyle->SetOptStat(0); h_qinv_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix->Draw();

    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_mix_cor->AddEntry(h_qinvCor_mix_2l, "Mix Cor - Double loop", "l");
    c_mix_cor->cd(); gStyle->SetOptStat(0); h_qinvCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor->Draw();
    
    Int_t numHistograms = 8;
    TH1D *histograms[] = {h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l, h_qinvDivCor_signal_2l, h_qinvDiv_signal_2l, h_qinv_mix_2l, h_qinvCor_mix_2l};

    TLegend *legends[] = {legend_sig_SS, legend_sig_OS, legend_sig_cor_OS, legend_sig_cor_SS, legend_sig_cor_Div, legend_sig_Div};
    int numLegends = 8;
    
    char prefix[50];
    sprintf(prefix, "signal_mix_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_mix, duration_signal};
    std::vector<std::string> labels = {"Total Time", "Mixing Time", "Signal Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEvents);

    const char *hpath = "./data/sigal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/sigal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);
    
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);      
}

int main(){
    const char* inputFile = "data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root";  
    const char* treeName  = "demo/TreeMBUCC";
    std::vector<double> centralityBins = {0.5, 1, 5, 10, 30, 50, 70}; 
    ControlVar selectedControlVar = ControlVar::CENT;
    Float_t selectedVertexDistance = 2.0f;
    int poolSize = 10;   
    Int_t pairChargeMult = 1;

    std::cout << "Starting analysis for " << centralityBins.size() - 1 << " bin(s)." << std::endl;
    for (size_t i = 0; i < centralityBins.size() - 1; ++i){
        double bin_low = centralityBins[i];
        double bin_high = centralityBins[i+1];
        
        std::cout << "\n--- Running for centrality bin: " << bin_low << "% - " << bin_high << "% ---" << std::endl;
        
        makeSignalMixDLM(inputFile, treeName, bin_low, bin_high, 
                     selectedControlVar, selectedVertexDistance, pairChargeMult, poolSize);
    }
    
    std::cout << "\nAnalysis finished." << std::endl;
    return 0;
}
