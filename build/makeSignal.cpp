#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <thread>
#include <string>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "ROOT/RConfig.hxx"
#include "../include/data_func.h"
#include "../include/process_func.h"


void sig_qinv_double_loop(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, 
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) 
{
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
    t->SetBranchAddress("pvZ", &pvZ); // Z in cm; Use fabs(pvZ) for absolute value of the vertex position in z axis
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

    // Getting how many entries
    Long64_t nentries = t->GetEntries(), trackvec_max_size = 0;
    std::cout << "#Events: " << nentries << std::endl;
    
    // See how many processed events
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMA ===
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;
    //double nnscale = numBins(ninterval, nlength, nscale), x0 = 0., xt=10.;
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

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;

    size_t poolSize = static_cast<size_t>(poolSizeInt);

    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        size_t size_before = vectorEventTracks4V.size();
        if (pvZ < vertexDistance) {
            if (!currentEventTracks4V.empty()) {
                if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()) { 
                    for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); nEv++) {
                        for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); p1++) {
                            for (size_t p2 = 0; p2 < currentEventTracks4V.size(); p2++) {
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
                vectorEventTracksWeight.push_back(currentEventTracksWeight);
                
                if (vectorEventTracks4V.size() >= poolSize) {
                    vectorEventTracks4V.clear();
                    vectorEventTracksCharge.clear();
                    vectorEventTracksWeight.clear();
                }
            }
        }
        size_t size_after = vectorEventTracks4V.size();
        if (size_after == size_before+1) processedEventsMix++;
        if (size_before == poolSize-1 && size_after == 0) processedEventsMix++;
        
        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        if (currentEventTracks4V.size() <= 1) return; // check if event has 2 or more tracks 

        for (size_t p1 = 0; p1 < currentEventTracks4V.size(); p1++) {
            for (size_t p2 = p1+1; p2 < currentEventTracks4V.size(); p2++) {
                // Checks
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue; // Check if weight is infinity

                // Calculate q_inv
                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                // Use a lock_guard to acquire the mutex. It will be automatically released.
                if (trkCharge[p1] * trkCharge[p2] > 0) { // Fills same charge particle pair histogram
                    h_qinvSS_signal_2l->Fill(qinv);
                    h_qinvSSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else { // Fills opposite charge particle pair histogram
                    h_qinvOS_signal_2l->Fill(qinv);
                    h_qinvOSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }

            }
        }

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
        
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;
    
    TH1D *h_qinvDivCor_signal_2l = (TH1D *)h_qinvSSCor_signal_2l->Clone("h_qinvDivCor_signal_2l");
    TH1D *h_qinvDiv_signal_2l = (TH1D *)h_qinvSS_signal_2l->Clone("h_qinvDiv_signal_2l");    

    h_qinvDivCor_signal_2l->Divide(h_qinvOSCor_signal_2l);
    h_qinvDiv_signal_2l->Divide(h_qinvOS_signal_2l);

    h_qinvDiv_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDiv_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvDivCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDivCor_signal_2l->GetYaxis()->SetTitle("#Pairs");

    h_qinvDivCor_signal_2l->Sumw2();

    // Save histograms
    TCanvas *c_sig_SS = new TCanvas("c_sig_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_OS = new TCanvas("c_sig_OS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_Div = new TCanvas("c_sig_cor_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_Div = new TCanvas("c_sig_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_mix = new TCanvas("c_mix", "Mixing Distributions", 1200, 800);
    TCanvas *c_mix_cor = new TCanvas("c_mix_cor", "Mixing Corrected Distributions", 1200, 800);
    
    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_SS->AddEntry(h_qinvSS_signal_2l, "Signal SS - One loop", "l");
    legend_sig_OS->AddEntry(h_qinvOS_signal_2l, "Signal OS - One loop", "l");
    legend_sig_cor_SS->AddEntry(h_qinvSSCor_signal_2l, "Signal SS Cor - One loop", "l");
    legend_sig_cor_OS->AddEntry(h_qinvOSCor_signal_2l, "Signal OS Cor - One loop", "l");
    legend_sig_cor_Div->AddEntry(h_qinvOSCor_signal_2l, "Signal Div Cor - One loop", "l");
    legend_sig_Div->AddEntry(h_qinvOSCor_signal_2l, "Signal Div - One loop", "l");
    legend_mix->AddEntry(h_qinv_mix_2l, "Mix - Double loop", "l");
    legend_mix_cor->AddEntry(h_qinvCor_mix_2l, "Mix Cor - Double loop", "l");
    
    legend_sig_SS->SetFillStyle(0); legend_sig_SS->SetBorderSize(0);
    legend_sig_OS->SetFillStyle(0); legend_sig_OS->SetBorderSize(0);
    legend_sig_cor_SS->SetFillStyle(0); legend_sig_cor_SS->SetBorderSize(0);
    legend_sig_cor_OS->SetFillStyle(0); legend_sig_cor_OS->SetBorderSize(0);
    legend_sig_cor_Div->SetFillStyle(0); legend_sig_cor_Div->SetBorderSize(0);
    legend_sig_Div->SetFillStyle(0); legend_sig_Div->SetBorderSize(0);
    legend_mix->SetFillStyle(0); legend_mix->SetBorderSize(0);
    legend_mix_cor->SetFillStyle(0); legend_mix_cor->SetBorderSize(0);

    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qinvSS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS");
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qinvOS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS");
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS Cor");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS Cor");
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qinvDivCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div Cor");
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qinvDiv_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div");
    c_mix->cd(); gStyle->SetOptStat(0); h_qinv_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing");
    c_mix_cor->cd(); gStyle->SetOptStat(0); h_qinvCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing Corrected");

    TH1D *histograms[] = {h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l, h_qinvDivCor_signal_2l, h_qinvDiv_signal_2l, h_qinv_mix_2l, h_qinvCor_mix_2l};
    Int_t numHistograms = 8;

    TCanvas *canvases[] = {c_sig_SS, c_sig_OS, c_sig_cor_OS, c_sig_cor_SS, c_sig_cor_Div, c_sig_Div, c_mix, c_mix_cor};
    int numCanvases = 8;

    TLegend *legends[] = {legend_sig_SS, legend_sig_OS, legend_sig_cor_OS, legend_sig_cor_SS, legend_sig_cor_Div, legend_sig_Div, legend_mix, legend_mix_cor};
    int numLegends = 8;

    // Prefix
    char prefix[100];
    sprintf(prefix, "sig_qinv_double_loop_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    // Saving Benchmakrs
    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);

    // Saving histograms
    const char *hpath = "./data/signal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    // Saving image
    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);

    AnalysisLog::instance().save("./logs", "makeSignalQinv");
    
    // Closing program
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_qinv_double_loop_parallel(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, 
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) 
{
    
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
    int thread_count = std::thread::hardware_concurrency();

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    
    int processedEventsSig = 0;
    int processedEventsMix = 0;

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

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;

    size_t poolSize = static_cast<size_t>(poolSizeInt);

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        size_t size_before = vectorEventTracks4V.size();
        processMixQinv(
            currentEventTracks4V.size(),
            thread_count,
            poolSize,
            vertexDistance,
            pairChargeMult,
            pvZ,
            vectorEventTracks4V,
            vectorEventTracksCharge,
            vectorEventTracksWeight,
            currentEventTracks4V,
            currentEventTracksCharge,
            currentEventTracksWeight,
            h_qinv_mix_2l,
            h_qinvCor_mix_2l
        );
        size_t size_after = vectorEventTracks4V.size();
        if (size_after == size_before+1) processedEventsMix++;
        if (size_before == poolSize-1 && size_after == 0) processedEventsMix++;

        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        // === SIGNAL (Parallelized) ===
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        processSignalQinv(
            currentEventTracks4V.size(),
            thread_count,
            trkWeight, 
            trkCharge, 
            currentEventTracks4V, 
            h_qinvSS_signal_2l, 
            h_qinvSSCor_signal_2l, 
            h_qinvOS_signal_2l, 
            h_qinvOSCor_signal_2l
        );

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;
    
    TH1D *h_qinvDivCor_signal_2l = (TH1D *)h_qinvSSCor_signal_2l->Clone("h_qinvDivCor_signal_2l");
    TH1D *h_qinvDiv_signal_2l = (TH1D *)h_qinvSS_signal_2l->Clone("h_qinvDiv_signal_2l"); 

    h_qinvDivCor_signal_2l->Divide(h_qinvOSCor_signal_2l);
    h_qinvDiv_signal_2l->Divide(h_qinvOS_signal_2l);

    h_qinvDiv_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDiv_signal_2l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvDivCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDivCor_signal_2l->GetYaxis()->SetTitle("#Pairs");

    h_qinvDivCor_signal_2l->Sumw2();

    // Save histograms
    TCanvas *c_sig_SS = new TCanvas("c_sig_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_OS = new TCanvas("c_sig_OS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_Div = new TCanvas("c_sig_cor_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_Div = new TCanvas("c_sig_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_mix = new TCanvas("c_mix", "Mixing Distributions", 1200, 800);
    TCanvas *c_mix_cor = new TCanvas("c_mix_cor", "Mixing Corrected Distributions", 1200, 800);

    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_SS->AddEntry(h_qinvSS_signal_2l, "Signal SS - One loop", "l");
    legend_sig_OS->AddEntry(h_qinvOS_signal_2l, "Signal OS - One loop", "l");
    legend_sig_cor_SS->AddEntry(h_qinvSSCor_signal_2l, "Signal SS Cor - One loop", "l");
    legend_sig_cor_OS->AddEntry(h_qinvOSCor_signal_2l, "Signal OS Cor - One loop", "l");
    legend_sig_cor_Div->AddEntry(h_qinvOSCor_signal_2l, "Signal Div Cor - One loop", "l");
    legend_sig_Div->AddEntry(h_qinvOSCor_signal_2l, "Signal Div - One loop", "l");
    legend_mix->AddEntry(h_qinv_mix_2l, "Mix - Double loop", "l");
    legend_mix_cor->AddEntry(h_qinvCor_mix_2l, "Mix Cor - Double loop", "l");
    
    legend_sig_SS->SetFillStyle(0); legend_sig_SS->SetBorderSize(0);
    legend_sig_OS->SetFillStyle(0); legend_sig_OS->SetBorderSize(0);
    legend_sig_cor_SS->SetFillStyle(0); legend_sig_cor_SS->SetBorderSize(0);
    legend_sig_cor_OS->SetFillStyle(0); legend_sig_cor_OS->SetBorderSize(0);
    legend_sig_cor_Div->SetFillStyle(0); legend_sig_cor_Div->SetBorderSize(0);
    legend_sig_Div->SetFillStyle(0); legend_sig_Div->SetBorderSize(0);
    legend_mix->SetFillStyle(0); legend_mix->SetBorderSize(0);
    legend_mix_cor->SetFillStyle(0); legend_mix_cor->SetBorderSize(0);

    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qinvSS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS");
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qinvOS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS");
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS Cor");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS Cor");
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qinvDivCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div Cor");
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qinvDiv_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div");
    c_mix->cd(); gStyle->SetOptStat(0); h_qinv_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing");
    c_mix_cor->cd(); gStyle->SetOptStat(0); h_qinvCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing Corrected");

    TH1D *histograms[] = {h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l, h_qinvDivCor_signal_2l, h_qinvDiv_signal_2l, h_qinv_mix_2l, h_qinvCor_mix_2l};
    Int_t numHistograms = 8;

    TCanvas *canvases[] = {c_sig_SS, c_sig_OS, c_sig_cor_OS, c_sig_cor_SS, c_sig_cor_Div, c_sig_Div, c_mix, c_mix_cor};
    int numCanvases = 8;

    TLegend *legends[] = {legend_sig_SS, legend_sig_OS, legend_sig_cor_OS, legend_sig_cor_SS, legend_sig_cor_Div, legend_sig_Div, legend_mix, legend_mix_cor};
    int numLegends = 8;

    char prefix[100];
    sprintf(prefix, "sig_qinv_double_loop_parallel_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);
    
    const char *hpath = "./data/signal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);

    AnalysisLog::instance().save("./logs", "makeSignalQinvParallel");
    
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_qinv_double_loop_parallel_pt(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT,
    Double_t ptMin = 0.5, 
    Double_t ptMax = -1.0,
    Double_t etaMin = -0.95,
    Double_t etaMax = 0.95,  
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    int poolSizeInt = 10) 
{
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
    int thread_count = std::thread::hardware_concurrency();

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMS ===
    double nnscale = 10000, x0 = 0., xt=10.;

    TH1D* h_qinvSSCor_signal_2l = new TH1D("h_qinvSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_signal_2l = new TH1D("h_qinvOSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_mix_2l = new TH1D("h_qinvSSCor_mix_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_mix_2l = new TH1D("h_qinvOSCor_mix_2l", "", nnscale, x0, xt);

    h_qinvSSCor_signal_2l->Sumw2();
    h_qinvOSCor_signal_2l->Sumw2();
    h_qinvSSCor_mix_2l->Sumw2();
    h_qinvOSCor_mix_2l->Sumw2();
    
    
    h_qinvSSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvOSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvSSCor_mix_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_mix_2l->GetYaxis()->SetTitle("#Pairs");
    h_qinvOSCor_mix_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOSCor_mix_2l->GetYaxis()->SetTitle("#Pairs");

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;
    
    size_t poolSize = static_cast<size_t>(poolSizeInt);

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        size_t size_before = vectorEventTracks4V.size();
        processMixQinv(
            currentEventTracks4V.size(),
            thread_count,
            poolSize,
            vertexDistance,
            pvZ,
            vectorEventTracks4V,
            vectorEventTracksCharge,
            vectorEventTracksWeight,
            currentEventTracks4V,
            currentEventTracksCharge,
            currentEventTracksWeight,
            h_qinvSSCor_mix_2l,
            h_qinvOSCor_mix_2l
        );
        size_t size_after = vectorEventTracks4V.size();
        if (size_after == size_before+1) processedEventsMix++;
        if (size_before == poolSize-1 && size_after == 0) processedEventsMix++;

        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        // === SIGNAL (Parallelized) ===
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        processSignalQinv(
            currentEventTracks4V.size(),
            thread_count,
            trkWeight, 
            trkCharge,
            currentEventTracks4V, 
            h_qinvSSCor_signal_2l, 
            h_qinvOSCor_signal_2l
        );

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;

    // Save histograms
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Corrected Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Corrected Distributions", 1200, 800);
    TCanvas *c_mix_cor_SS = new TCanvas("c_mix_cor_SS", "Mixing Corrected Distributions", 1200, 800);
    TCanvas *c_mix_cor_OS = new TCanvas("c_mix_cor_OS", "Mixing Corrected Distributions", 1200, 800);

    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_cor_SS->AddEntry(h_qinvSSCor_signal_2l, "Signal SS Cor - Double loop", "l");
    legend_sig_cor_OS->AddEntry(h_qinvOSCor_signal_2l, "Signal OS Cor - Double loop", "l");
    legend_mix_cor_SS->AddEntry(h_qinvSSCor_mix_2l, "Mix SS Cor- Double loop", "l");
    legend_mix_cor_OS->AddEntry(h_qinvOSCor_mix_2l, "Mix OS Cor - Double loop", "l");

    legend_sig_cor_SS->SetFillStyle(0); legend_sig_cor_SS->SetBorderSize(0);
    legend_sig_cor_OS->SetFillStyle(0); legend_sig_cor_OS->SetBorderSize(0);
    legend_mix_cor_SS->SetFillStyle(0); legend_mix_cor_SS->SetBorderSize(0);
    legend_mix_cor_OS->SetFillStyle(0); legend_mix_cor_OS->SetBorderSize(0);

    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal SS Cor");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal OS Cor");
    c_mix_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Mixing SS Cor");
    c_mix_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Mixing OS Cor");

    TH1D *histograms[] = {h_qinvSSCor_signal_2l, h_qinvOSCor_signal_2l, h_qinvSSCor_mix_2l, h_qinvOSCor_mix_2l};
    Int_t numHistograms = 4;

    TCanvas *canvases[] = {c_sig_cor_OS, c_sig_cor_SS, c_mix_cor_SS, c_mix_cor_OS};
    int numCanvases = 4;

    TLegend *legends[] = {legend_sig_cor_OS, legend_sig_cor_SS, legend_mix_cor_SS, legend_mix_cor_OS};
    int numLegends = 4;

    char prefix[100];
    if (ptMax > 0){
        if (ptMin > 0){
            sprintf(prefix,
                    "sig_qinv_eta-%f-to-%f_pT-from-%f-to-%f_%s_%f-%f",
                    etaMin, etaMax, ptMin, ptMax,
                    selectionVarName,
                    displaySelVarMoreeq, displaySelVarLess);
        } else {
            sprintf(prefix,
                    "sig_qinv_eta-%f-to-%f_pT-to-%f_%s_%f-%f",
                    etaMin, etaMax, ptMax,
                    selectionVarName,
                    displaySelVarMoreeq, displaySelVarLess);
        }
    } else {
        if (ptMin > 0){
            sprintf(prefix,
                    "sig_qinv_eta-%f-to-%f_pT-from-%f_%s_%f-%f",
                    etaMin, etaMax, ptMin,
                    selectionVarName,
                    displaySelVarMoreeq, displaySelVarLess);
        } else {
            sprintf(prefix,
                    "sig_qinv_eta-%f-to-%f_%s_%f-%f",
                    etaMin, etaMax,
                    selectionVarName,
                    displaySelVarMoreeq, displaySelVarLess);
        }
    } 
    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);
    
    const char *hpath = "./data/signal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);

    AnalysisLog::instance().save("./logs", "makeSignalQinvParallelPtEta");
    
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_qlcms_double_loop(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, 
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) 
{
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
    t->SetBranchAddress("pvZ", &pvZ); // Z in cm; Use fabs(pvZ) for absolute value of the vertex position in z axis
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

    // Getting how many entries
    Long64_t nentries = t->GetEntries(), trackvec_max_size = 0;
    std::cout << "#Events: " << nentries << std::endl;
    
    // See how many processed events
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMA ===
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;
    //double nnscale = numBins(ninterval, nlength, nscale), x0 = 0., xt=10.;
    double nnscale = 10000, x0 = 0., xt=10.;

    TH1D* h_qlcmsSS_signal_2l = new TH1D("h_qlcmsSS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsSSCor_signal_2l = new TH1D("h_qlcmsSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsOS_signal_2l = new TH1D("h_qlcmsOS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsOSCor_signal_2l = new TH1D("h_qlcmsOSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcms_mix_2l = new TH1D("h_qlcms_mix_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsCor_mix_2l = new TH1D("h_qlcmsCor_mix_2l", "", nnscale, x0, xt);

    h_qlcmsSSCor_signal_2l->Sumw2();
    h_qlcmsOSCor_signal_2l->Sumw2();
    h_qlcmsCor_mix_2l->Sumw2();

    h_qlcmsSS_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsSS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsOS_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsOS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsSSCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsOSCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsOSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcms_mix_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcms_mix_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsCor_mix_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsCor_mix_2l->GetYaxis()->SetTitle("#Pairs");

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;

    size_t poolSize = static_cast<size_t>(poolSizeInt);

    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        size_t size_before = vectorEventTracks4V.size();
        if (pvZ < vertexDistance) {
            if (!currentEventTracks4V.empty()) {
                if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()) { 
                    for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); nEv++) {
                        for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); p1++) {
                            for (size_t p2 = 0; p2 < currentEventTracks4V.size(); p2++) {
                                if (!(vectorEventTracksCharge[nEv][p1]*trkCharge[p2] == pairChargeMult)) continue;
                                if (std::isinf(vectorEventTracksWeight[nEv][p1]*trkWeight[p2])) continue;
                                
                                double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                                h_qlcms_mix_2l->Fill(qinv);
                                h_qlcmsCor_mix_2l->Fill(qinv, vectorEventTracksWeight[nEv][p1]*trkWeight[p2]);
                            }
                        }
                    }
                }
                vectorEventTracks4V.push_back(currentEventTracks4V);
                vectorEventTracksCharge.push_back(currentEventTracksCharge);
                vectorEventTracksWeight.push_back(currentEventTracksWeight);
                
                if (vectorEventTracks4V.size() >= poolSize) {
                    vectorEventTracks4V.clear();
                    vectorEventTracksCharge.clear();
                    vectorEventTracksWeight.clear();
                }
            }
        }
        size_t size_after = vectorEventTracks4V.size();
        if (size_after == size_before+1) processedEventsMix++;
        if (size_before == poolSize-1 && size_after == 0) processedEventsMix++;
        
        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        if (currentEventTracks4V.size() <= 1) return; // check if event has 2 or more tracks 

        for (size_t p1 = 0; p1 < currentEventTracks4V.size(); p1++) {
            for (size_t p2 = p1+1; p2 < currentEventTracks4V.size(); p2++) {
                // Checks
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue; // Check if weight is infinity

                // Calculate q_inv
                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                // Use a lock_guard to acquire the mutex. It will be automatically released.
                if (trkCharge[p1] * trkCharge[p2] > 0) { // Fills same charge particle pair histogram
                    h_qlcmsSS_signal_2l->Fill(qinv);
                    h_qlcmsSSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else { // Fills opposite charge particle pair histogram
                    h_qlcmsOS_signal_2l->Fill(qinv);
                    h_qlcmsOSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }

            }
        }

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
        
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;
    
    TH1D *h_qlcmsDivCor_signal_2l = (TH1D *)h_qlcmsSSCor_signal_2l->Clone("h_qlcmsDivCor_signal_2l");
    TH1D *h_qlcmsDiv_signal_2l = (TH1D *)h_qlcmsSS_signal_2l->Clone("h_qlcmsDiv_signal_2l");    

    h_qlcmsDivCor_signal_2l->Divide(h_qlcmsOSCor_signal_2l);
    h_qlcmsDiv_signal_2l->Divide(h_qlcmsOS_signal_2l);
    h_qlcmsDiv_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsDiv_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsDivCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsDivCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsDivCor_signal_2l->Sumw2();

    // Save histograms
    TCanvas *c_sig_SS = new TCanvas("c_sig_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_OS = new TCanvas("c_sig_OS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_Div = new TCanvas("c_sig_cor_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_Div = new TCanvas("c_sig_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_mix = new TCanvas("c_mix", "Mixing Distributions", 1200, 800);
    TCanvas *c_mix_cor = new TCanvas("c_mix_cor", "Mixing Corrected Distributions", 1200, 800);
    
    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_SS->AddEntry(h_qlcmsSS_signal_2l, "Signal SS - One loop", "l");
    legend_sig_OS->AddEntry(h_qlcmsOS_signal_2l, "Signal OS - One loop", "l");
    legend_sig_cor_SS->AddEntry(h_qlcmsSSCor_signal_2l, "Signal SS Cor - One loop", "l");
    legend_sig_cor_OS->AddEntry(h_qlcmsOSCor_signal_2l, "Signal OS Cor - One loop", "l");
    legend_sig_cor_Div->AddEntry(h_qlcmsOSCor_signal_2l, "Signal Div Cor - One loop", "l");
    legend_sig_Div->AddEntry(h_qlcmsOSCor_signal_2l, "Signal Div - One loop", "l");
    legend_mix->AddEntry(h_qlcms_mix_2l, "Mix - Double loop", "l");
    legend_mix_cor->AddEntry(h_qlcmsCor_mix_2l, "Mix Cor - Double loop", "l");
    
    legend_sig_SS->SetFillStyle(0); legend_sig_SS->SetBorderSize(0);
    legend_sig_OS->SetFillStyle(0); legend_sig_OS->SetBorderSize(0);
    legend_sig_cor_SS->SetFillStyle(0); legend_sig_cor_SS->SetBorderSize(0);
    legend_sig_cor_OS->SetFillStyle(0); legend_sig_cor_OS->SetBorderSize(0);
    legend_sig_cor_Div->SetFillStyle(0); legend_sig_cor_Div->SetBorderSize(0);
    legend_sig_Div->SetFillStyle(0); legend_sig_Div->SetBorderSize(0);
    legend_mix->SetFillStyle(0); legend_mix->SetBorderSize(0);
    legend_mix_cor->SetFillStyle(0); legend_mix_cor->SetBorderSize(0);

    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qlcmsSS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS");
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qlcmsOS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS");
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qlcmsSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS Cor");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qlcmsOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS Cor");
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qlcmsDivCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div Cor");
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qlcmsDiv_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div");
    c_mix->cd(); gStyle->SetOptStat(0); h_qlcms_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing");
    c_mix_cor->cd(); gStyle->SetOptStat(0); h_qlcmsCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing Corrected");

    TH1D *histograms[] = {h_qlcmsSS_signal_2l, h_qlcmsSSCor_signal_2l, h_qlcmsOS_signal_2l, h_qlcmsOSCor_signal_2l, h_qlcmsDivCor_signal_2l, h_qlcmsDiv_signal_2l, h_qlcms_mix_2l, h_qlcmsCor_mix_2l};
    Int_t numHistograms = 8;

    TCanvas *canvases[] = {c_sig_SS, c_sig_OS, c_sig_cor_OS, c_sig_cor_SS, c_sig_cor_Div, c_sig_Div, c_mix, c_mix_cor};
    int numCanvases = 8;

    TLegend *legends[] = {legend_sig_SS, legend_sig_OS, legend_sig_cor_OS, legend_sig_cor_SS, legend_sig_cor_Div, legend_sig_Div, legend_mix, legend_mix_cor};
    int numLegends = 8;

    // Prefix
    char prefix[100];
    sprintf(prefix, "sig_qlcms_double_loop_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    // Saving Benchmakrs
    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);

    // Saving histograms
    const char *hpath = "./data/signal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    // Saving image
    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);
    
    AnalysisLog::instance().save("./logs", "makeSignalQLCMS");

    // Closing program
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_qlcms_double_loop_parallel(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, 
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) 
{
    
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
    int thread_count = std::thread::hardware_concurrency();

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMS ===
    double nnscale = 10000, x0 = 0., xt=10.;

    TH1D* h_qlcmsSS_signal_2l = new TH1D("h_qlcmsSS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsSSCor_signal_2l = new TH1D("h_qlcmsSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsOS_signal_2l = new TH1D("h_qlcmsOS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsOSCor_signal_2l = new TH1D("h_qlcmsOSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcms_mix_2l = new TH1D("h_qlcms_mix_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsCor_mix_2l = new TH1D("h_qlcmsCor_mix_2l", "", nnscale, x0, xt);

    h_qlcmsSSCor_signal_2l->Sumw2();
    h_qlcmsOSCor_signal_2l->Sumw2();
    h_qlcmsCor_mix_2l->Sumw2();
    
    h_qlcmsSS_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsSS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsOS_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qlcmsOS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsSSCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsOSCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsOSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcms_mix_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcms_mix_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsCor_mix_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsCor_mix_2l->GetYaxis()->SetTitle("#Pairs");

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;
    
    size_t poolSize = static_cast<size_t>(poolSizeInt);

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        size_t size_before = vectorEventTracks4V.size();
        processMixQLCMS(
            currentEventTracks4V.size(),
            thread_count,
            poolSize,
            vertexDistance,
            pairChargeMult,
            pvZ,
            vectorEventTracks4V,
            vectorEventTracksCharge,
            vectorEventTracksWeight,
            currentEventTracks4V,
            currentEventTracksCharge,
            currentEventTracksWeight,
            h_qlcms_mix_2l,
            h_qlcmsCor_mix_2l
        );
        size_t size_after = vectorEventTracks4V.size();
        if (size_after == size_before+1) processedEventsMix++;
        if (size_before == poolSize-1 && size_after == 0) processedEventsMix++;

        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        // === SIGNAL (Parallelized) ===
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        processSignalQLCMS(
            currentEventTracks4V.size(),
            thread_count,
            trkWeight, 
            trkCharge, 
            currentEventTracks4V, 
            h_qlcmsSS_signal_2l, 
            h_qlcmsSSCor_signal_2l, 
            h_qlcmsOS_signal_2l, 
            h_qlcmsOSCor_signal_2l
        );

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;
    
    TH1D *h_qlcmsDivCor_signal_2l = (TH1D *)h_qlcmsSSCor_signal_2l->Clone("h_qlcmsDivCor_signal_2l");
    TH1D *h_qlcmsDiv_signal_2l = (TH1D *)h_qlcmsSS_signal_2l->Clone("h_qlcmsDiv_signal_2l"); 

    h_qlcmsDivCor_signal_2l->Divide(h_qlcmsOSCor_signal_2l);
    h_qlcmsDiv_signal_2l->Divide(h_qlcmsOS_signal_2l);
    h_qlcmsDiv_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsDiv_signal_2l->GetYaxis()->SetTitle("#Pairs");

    h_qlcmsDivCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsDivCor_signal_2l->GetYaxis()->SetTitle("#Pairs");

    h_qlcmsDivCor_signal_2l->Sumw2();

    // Save histograms
    TCanvas *c_sig_SS = new TCanvas("c_sig_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_OS = new TCanvas("c_sig_OS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_cor_Div = new TCanvas("c_sig_cor_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_sig_Div = new TCanvas("c_sig_Div", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_mix = new TCanvas("c_mix", "Mixing Distributions", 1200, 800);
    TCanvas *c_mix_cor = new TCanvas("c_mix_cor", "Mixing Corrected Distributions", 1200, 800);

    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_SS->AddEntry(h_qlcmsSS_signal_2l, "Signal SS - One loop", "l");
    legend_sig_OS->AddEntry(h_qlcmsOS_signal_2l, "Signal OS - One loop", "l");
    legend_sig_cor_SS->AddEntry(h_qlcmsSSCor_signal_2l, "Signal SS Cor - One loop", "l");
    legend_sig_cor_OS->AddEntry(h_qlcmsOSCor_signal_2l, "Signal OS Cor - One loop", "l");
    legend_sig_cor_Div->AddEntry(h_qlcmsOSCor_signal_2l, "Signal Div Cor - One loop", "l");
    legend_sig_Div->AddEntry(h_qlcmsOSCor_signal_2l, "Signal Div - One loop", "l");
    legend_mix->AddEntry(h_qlcms_mix_2l, "Mix - Double loop", "l");
    legend_mix_cor->AddEntry(h_qlcmsCor_mix_2l, "Mix Cor - Double loop", "l");

    legend_sig_SS->SetFillStyle(0); legend_sig_SS->SetBorderSize(0);
    legend_sig_OS->SetFillStyle(0); legend_sig_OS->SetBorderSize(0);
    legend_sig_cor_SS->SetFillStyle(0); legend_sig_cor_SS->SetBorderSize(0);
    legend_sig_cor_OS->SetFillStyle(0); legend_sig_cor_OS->SetBorderSize(0);
    legend_sig_cor_Div->SetFillStyle(0); legend_sig_cor_Div->SetBorderSize(0);
    legend_sig_Div->SetFillStyle(0); legend_sig_Div->SetBorderSize(0);
    legend_mix->SetFillStyle(0); legend_mix->SetBorderSize(0);
    legend_mix_cor->SetFillStyle(0); legend_mix_cor->SetBorderSize(0);

    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qlcmsSS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS");
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qlcmsOS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS");
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qlcmsSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal SS Cor");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qlcmsOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal OS Cor");
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qlcmsDivCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div Cor");
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qlcmsDiv_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Signal Div");
    c_mix->cd(); gStyle->SetOptStat(0); h_qlcms_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing");
    c_mix_cor->cd(); gStyle->SetOptStat(0); h_qlcmsCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor->Draw();
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", "PbPb 2.76 TeV | Mixing Corrected");

    TH1D *histograms[] = {h_qlcmsSS_signal_2l, h_qlcmsSSCor_signal_2l, h_qlcmsOS_signal_2l, h_qlcmsOSCor_signal_2l, h_qlcmsDivCor_signal_2l, h_qlcmsDiv_signal_2l, h_qlcms_mix_2l, h_qlcmsCor_mix_2l};
    Int_t numHistograms = 8;

    TCanvas *canvases[] = {c_sig_SS, c_sig_OS, c_sig_cor_OS, c_sig_cor_SS, c_sig_cor_Div, c_sig_Div, c_mix, c_mix_cor};
    int numCanvases = 8;

    TLegend *legends[] = {legend_sig_SS, legend_sig_OS, legend_sig_cor_OS, legend_sig_cor_SS, legend_sig_cor_Div, legend_sig_Div, legend_mix, legend_mix_cor};
    int numLegends = 8;

    char prefix[100];
    sprintf(prefix, "sig_qlcms_double_loop_parallel_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);
    
    const char *hpath = "./data/signal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);
    
    AnalysisLog::instance().save("./logs", "makeSignalQLCMSParallel");

    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_qlcms_double_loop_parallel_pt(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT,
    Double_t ptMin = 0.5, 
    Double_t ptMax = -1.0,
    Double_t etaMin = -0.95,
    Double_t etaMax = 0.95,  
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    int poolSizeInt = 10) 
{
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
    int thread_count = std::thread::hardware_concurrency();

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMS ===
    double nnscale = 10000, x0 = 0., xt=10.;

    TH1D* h_qlcmsSSCor_signal_2l = new TH1D("h_qlcmsSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsOSCor_signal_2l = new TH1D("h_qlcmsOSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsSSCor_mix_2l = new TH1D("h_qlcmsSSCor_mix_2l", "", nnscale, x0, xt);
    TH1D* h_qlcmsOSCor_mix_2l = new TH1D("h_qlcmsOSCor_mix_2l", "", nnscale, x0, xt);

    h_qlcmsSSCor_signal_2l->Sumw2();
    h_qlcmsOSCor_signal_2l->Sumw2();
    h_qlcmsSSCor_mix_2l->Sumw2();
    h_qlcmsOSCor_mix_2l->Sumw2();
    
    h_qlcmsSSCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsOSCor_signal_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsOSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsSSCor_mix_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsSSCor_mix_2l->GetYaxis()->SetTitle("#Pairs");
    h_qlcmsOSCor_mix_2l->GetXaxis()->SetTitle("q_{LCMS}[GeV]");
    h_qlcmsOSCor_mix_2l->GetYaxis()->SetTitle("#Pairs");

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;
    
    size_t poolSize = static_cast<size_t>(poolSizeInt);

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        if (std::abs(pvZ) >= vertexDistance) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        // === MIXING === 
        size_t size_before = vectorEventTracks4V.size();
        processMixQLCMS(
            currentEventTracks4V.size(),
            thread_count,
            poolSize,
            vertexDistance,
            pvZ,
            vectorEventTracks4V,
            vectorEventTracksCharge,
            vectorEventTracksWeight,
            currentEventTracks4V,
            currentEventTracksCharge,
            currentEventTracksWeight,
            h_qlcmsSSCor_mix_2l,
            h_qlcmsOSCor_mix_2l
        );
        size_t size_after = vectorEventTracks4V.size();
        if (size_after == size_before+1) processedEventsMix++;
        if (size_before == poolSize-1 && size_after == 0) processedEventsMix++;

        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();
        
        // === SIGNAL (Parallelized) ===
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        processSignalQLCMS(
            currentEventTracks4V.size(),
            thread_count,
            trkWeight, 
            trkCharge,
            currentEventTracks4V, 
            h_qlcmsSSCor_signal_2l, 
            h_qlcmsOSCor_signal_2l
        );

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;

    // Save histograms
    TCanvas *c_sig_cor_SS = new TCanvas("c_sig_cor_SS", "Signal Corrected Distributions", 1200, 800);
    TCanvas *c_sig_cor_OS = new TCanvas("c_sig_cor_OS", "Signal Corrected Distributions", 1200, 800);
    TCanvas *c_mix_cor_SS = new TCanvas("c_mix_cor_SS", "Mixing Corrected Distributions", 1200, 800);
    TCanvas *c_mix_cor_OS = new TCanvas("c_mix_cor_OS", "Mixing Corrected Distributions", 1200, 800);

    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_cor_SS->AddEntry(h_qlcmsSSCor_signal_2l, "Signal SS Cor - Double loop", "l");
    legend_sig_cor_OS->AddEntry(h_qlcmsOSCor_signal_2l, "Signal OS Cor - Double loop", "l");
    legend_mix_cor_SS->AddEntry(h_qlcmsSSCor_mix_2l, "Mix SS Cor- Double loop", "l");
    legend_mix_cor_OS->AddEntry(h_qlcmsOSCor_mix_2l, "Mix OS Cor - Double loop", "l");

    legend_sig_cor_SS->SetFillStyle(0); legend_sig_cor_SS->SetBorderSize(0);
    legend_sig_cor_OS->SetFillStyle(0); legend_sig_cor_OS->SetBorderSize(0);
    legend_mix_cor_SS->SetFillStyle(0); legend_mix_cor_SS->SetBorderSize(0);
    legend_mix_cor_OS->SetFillStyle(0); legend_mix_cor_OS->SetBorderSize(0);

    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qlcmsSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal SS Cor");
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qlcmsOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal OS Cor");
    c_mix_cor_SS->cd(); gStyle->SetOptStat(0); h_qlcmsSSCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor_SS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Mixing SS Cor");
    c_mix_cor_OS->cd(); gStyle->SetOptStat(0); h_qlcmsOSCor_mix_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor_OS->Draw();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Mixing OS Cor");

    TH1D *histograms[] = {h_qlcmsSSCor_signal_2l, h_qlcmsOSCor_signal_2l, h_qlcmsSSCor_mix_2l, h_qlcmsOSCor_mix_2l};
    Int_t numHistograms = 4;

    TCanvas *canvases[] = {c_sig_cor_OS, c_sig_cor_SS, c_mix_cor_SS, c_mix_cor_OS};
    int numCanvases = 4;

    TLegend *legends[] = {legend_sig_cor_OS, legend_sig_cor_SS, legend_mix_cor_SS, legend_mix_cor_OS};
    int numLegends = 4;

    char prefix[100];
    if (ptMax > 0){
        if (ptMin > 0){
            sprintf(prefix,
                    "sig_qlcms_eta-%.2f-to-%.2f_pT-from-%.2f-to-%.2f_%s",
                    etaMin, etaMax, ptMin, ptMax, selectionVarName);
        } else {
            sprintf(prefix,
                    "sig_qlcms_eta-%.2f-to-%.2f_pT-to-%.2f_%s",
                    etaMin, etaMax, ptMax, selectionVarName);
        }
    } else {
        if (ptMin > 0){
            sprintf(prefix,
                    "sig_qlcms_eta-%.2f-to-%.2f_pT-from-%.2f_%s",
                    etaMin, etaMax, ptMin, selectionVarName);
        } else {
            sprintf(prefix,
                    "sig_qlcms_eta-%.2f-to-%.2f_%s",
                    etaMin, etaMax, selectionVarName);
        }
    } 
    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);
    
    const char *hpath = "./data/signal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);

    AnalysisLog::instance().save("./logs", "makeSignalQLCMSParallelPtEta");
    
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_qtqzq0q(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT,
    Double_t ptMin = 0.5, 
    Double_t ptMax = -1.0,
    Double_t etaMin = -0.95,
    Double_t etaMax = 0.95,  
    int test_limit_sig = -1,
    int test_limit_mix = -1,
    Float_t vertexDistance = 2, 
    Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) 
{
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
    int thread_count = std::thread::hardware_concurrency();

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMS ===
    double nnscale = 1000, x0 = 0., xt=0.02;

    TH2D* h_qtqz_SS_cor = new TH2D("h_qtqz_SS_cor", "", nnscale, x0, xt, nnscale, x0, xt);
    TH2D* h_qtqz_OS_cor = new TH2D("h_qtqz_OS_cor", "", nnscale, x0, xt, nnscale, x0, xt);
    TH2D* h_q0q_SS_cor = new TH2D("h_q0q_SS_cor", "", nnscale, x0, xt, nnscale, x0, xt);
    TH2D* h_q0q_OS_cor = new TH2D("h_q0q_OS_cor", "", nnscale, x0, xt, nnscale, x0, xt);

    h_qtqz_SS_cor->Sumw2();
    h_qtqz_OS_cor->Sumw2();
    h_q0q_SS_cor->Sumw2();
    h_q0q_OS_cor->Sumw2();
    
    h_qtqz_SS_cor->GetXaxis()->SetTitle("q_{t} [GeV]");
    h_qtqz_SS_cor->GetYaxis()->SetTitle("q_{z} [GeV]");
    h_qtqz_SS_cor->GetZaxis()->SetTitle("#Pairs");
    h_qtqz_OS_cor->GetXaxis()->SetTitle("q_{t} [GeV]");
    h_qtqz_OS_cor->GetYaxis()->SetTitle("q_{z} [GeV]");
    h_qtqz_OS_cor->GetZaxis()->SetTitle("#Pairs");
    h_q0q_SS_cor->GetXaxis()->SetTitle("q_{0} [GeV]");
    h_q0q_SS_cor->GetYaxis()->SetTitle("|#vec{q}| [GeV]");
    h_q0q_SS_cor->GetZaxis()->SetTitle("#Pairs");
    h_q0q_OS_cor->GetXaxis()->SetTitle("q_{0} [GeV]");
    h_q0q_OS_cor->GetYaxis()->SetTitle("|#vec{q}| [GeV]");
    h_q0q_OS_cor->GetZaxis()->SetTitle("#Pairs");

    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;
    std::vector<std::vector<Float_t>> vectorEventTracksPt;
    
    size_t poolSize = static_cast<size_t>(poolSizeInt);

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        hiBinProxy = hiBin;
        NtrkProxy = Ntrk;
        HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig) break;
        if (processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        processedEventsSig++;
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTracksWeight;
        std::vector<Float_t> currentEventTracksPt;

        for (int j = 0; j < Ntrk; j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTracksWeight.push_back(trkWeight[j]);
            currentEventTracksPt.push_back(trkPt[j]);
        }
        
        
        // === SIGNAL (Parallelized) ===
        auto start_signal_lap = std::chrono::high_resolution_clock::now();
        processSignalqtqzq0q(
            currentEventTracks4V.size(),
            thread_count,
            trkWeight, 
            trkCharge,
            currentEventTracks4V, 
            h_qtqz_SS_cor, 
            h_qtqz_OS_cor, 
            h_q0q_SS_cor, 
            h_q0q_OS_cor
        );

        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;

    // Save histograms
    TCanvas *c_sig_qtqz_cor_SS = new TCanvas("c_sig_qtqz_cor_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_qtqz_cor_OS = new TCanvas("c_sig_qtqz_cor_OS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_q0q_cor_SS = new TCanvas("c_sig_q0q_cor_SS", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_q0q_cor_OS = new TCanvas("c_sig_q0q_cor_OS", "Signal Distributions", 1200, 800);

    c_sig_qtqz_cor_SS->cd(); gStyle->SetOptStat(0); h_qtqz_SS_cor->Draw("COLZ"); gStyle->SetPalette(kLake);
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal SS Cor | (q_{t}, q_{z})");
    c_sig_qtqz_cor_OS->cd(); gStyle->SetOptStat(0); h_qtqz_OS_cor->Draw("COLZ"); gStyle->SetPalette(kLake);
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal OS Cor | (q_{t}, q_{z})");
    c_sig_q0q_cor_SS->cd(); gStyle->SetOptStat(0); h_q0q_SS_cor->Draw("COLZ"); gStyle->SetPalette(kLake);
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal SS Cor | (q_{0}, q_{q})");
    c_sig_q0q_cor_OS->cd(); gStyle->SetOptStat(0); h_q0q_OS_cor->Draw("COLZ"); gStyle->SetPalette(kLake);
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | Signal OS Cor | (q_{0}, q_{q})");

    TH2D *histograms[] = {h_qtqz_SS_cor, h_qtqz_OS_cor, h_q0q_SS_cor, h_q0q_OS_cor};
    Int_t numHistograms = 4;

    TCanvas *canvases[] = {c_sig_qtqz_cor_SS, c_sig_qtqz_cor_OS, c_sig_q0q_cor_OS, c_sig_q0q_cor_SS};
    int numCanvases = 4;

    TLegend *legends[] = {};
    int numLegends = 0;

    char prefix[100];
    if (ptMax > 0){
        if (ptMin > 0){
            sprintf(prefix, "sig_qtqzq0q_eta-%f-to-%f_pT-from-%f-to-%f_%s_%f-%f", etaMin, etaMax, ptMin, ptMax, selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        } else {
            sprintf(prefix, "sig_qtqzq0q_eta-%f-to-%f_pT-to-%f_%s_%f-%f", etaMin, etaMax, ptMax, selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        }
    } else {
        if (ptMin > 0){
            sprintf(prefix, "sig_qtqzq0q_eta-%f-to-%f_pT-from-%f_%s_%f-%f", etaMin, etaMax, ptMin, selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        } else {
            sprintf(prefix, "sig_qtqzq0q_eta-%f-to-%f_%s_%f-%f", etaMin, etaMax, selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        }
    } 
    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix, processedEventsSig, processedEventsMix);
    
    const char *hpath = "./data/signal_mix/";
    save_histograms2d(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/signal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);

    AnalysisLog::instance().save("./logs", "makeSignalqtqzq0qParallelPtEta");
    
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);     
}

void sig_DeltaEtaDeltaPhi(
    const char *fileInput, 
    const char *treeInput, 
    double selVarMoreeq, 
    double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT,
    Float_t ptMin = 0.5, Float_t ptMax = -1.0,
    Float_t etaCut = 0.95,
    Float_t pixHit = 1, Float_t trkHit = 10,
    Float_t chi2 = 0.18,
    Float_t dzSig = 3.0, Float_t dxySig = 3.0,  
    int test_limit_sig = -1, int test_limit_mix = -1,
    Float_t vertexDistance = 15, 
    Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) 
{
    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree(fileInput, treeInput, fr, t);

    // Variables from the tree
    Int_t maxSize = 50000, Ntrk, Npv, hiBin, trkCharge[maxSize]; 
    Float_t HFsumET, pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize], 
            trkNpixLayers[maxSize], trkNvalidHits[maxSize], trkNormChi2[maxSize],
            trkDzSig[maxSize], trkDxySig[maxSize], pionMass = 0.13957039; // Pion mass [GeV] from PDG
    
    double* selectionVar, displaySelVarMoreeq, displaySelVarLess;
    double hiBinProxy, NtrkProxy, HFsumETProxy;
    int thread_count = std::thread::hardware_concurrency();

    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("Npv", &Npv);
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("pvZ", &pvZ);
    t->SetBranchAddress("trkNpixLayers", trkNpixLayers);
    t->SetBranchAddress("trkCharge", trkCharge);
    t->SetBranchAddress("trkWeight", trkWeight);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    t->SetBranchAddress("trkNvalidHits", trkNvalidHits);
    t->SetBranchAddress("trkNormChi2", trkNormChi2);
    t->SetBranchAddress("trkDzSig", trkDzSig);
    t->SetBranchAddress("trkDxySig", trkDxySig);

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    std::cout << "#Events: " << nentries << std::endl;
    
    int processedEventsSig = 0;
    int processedEventsMix = 0;

    // === HISTOGRAMS ===
    double nnscale = 100, x0 = 0., xt=0.15;

    TH2D* hSigSS_DeltaEtaDeltaPhi = new TH2D("hSigSS_DeltaEtaDeltaPhi", "", nnscale, x0, xt, nnscale, x0, xt);
    TH2D* hMixSS_DeltaEtaDeltaPhi = new TH2D("hMixSS_DeltaEtaDeltaPhi", "", nnscale, x0, xt, nnscale, x0, xt);

    hSigSS_DeltaEtaDeltaPhi->Sumw2();
    hMixSS_DeltaEtaDeltaPhi->Sumw2();
    
    hSigSS_DeltaEtaDeltaPhi->GetXaxis()->SetTitle("#Delta#eta");
    hMixSS_DeltaEtaDeltaPhi->GetXaxis()->SetTitle("#Delta#eta");
    hSigSS_DeltaEtaDeltaPhi->GetYaxis()->SetTitle("#Delta#varphi");
    hMixSS_DeltaEtaDeltaPhi->GetYaxis()->SetTitle("#Delta#varphi");
    hSigSS_DeltaEtaDeltaPhi->GetZaxis()->SetTitle("#Pairs");
    hMixSS_DeltaEtaDeltaPhi->GetZaxis()->SetTitle("#Pairs");
    
    double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;

    // 10 buckets for Z-vertex (e.g. -15 to +15 into 3cm bins)
    std::vector<EventData> pool[10]; 
    size_t poolSize = static_cast<size_t>(poolSizeInt);

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    auto start_full = std::chrono::high_resolution_clock::now();
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        if (i % 100 == 0 || i == nentries - 1) {
            printProgressBar(i + 1, nentries);
        }
        
        hiBinProxy = hiBin; NtrkProxy = Ntrk; HFsumETProxy = HFsumET;
        
        if (processedEventsSig == test_limit_sig || processedEventsMix == test_limit_mix) break;
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue;
        //if (Npv != 1) continue;
        if (std::abs(pvZ) >= vertexDistance) continue;
        int zBin = (int)((pvZ + vertexDistance) / ((vertexDistance * 2.0) / 10.0));
        if (zBin < 0) zBin = 0; if (zBin > 9) zBin = 9;
        
        EventData currentEv;
        for (int j = 0; j < Ntrk; j++) {
            // Kinematic cuts
            if (ptMax < 0){
                if (trkPt[j] < ptMin) continue;
            } else {
                if (!(trkPt[j] > ptMin && trkPt[j] < ptMax)) continue;
            }
            
            if (std::abs(trkEta[j]) > etaCut) continue;
            
            // Quality cuts (Hits and Chi2)
            if (trkNpixLayers[j] <= pixHit) continue;
            //if (trkNvalidHits[j] <= trkHit) continue;
            //if (trkNormChi2[j] >= chi2) continue;
            
            // DCA cuts (Distance of Closest Approach)
            if (std::abs(trkDzSig[j]) >= dzSig) continue;
            //if (std::abs(trkDxySig[j]) >= dxySig) continue;
            
            // Check for bugs
            if (std::isinf(trkWeight[j])) continue;
            
            currentEv.tracks.push_back(ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass));
            currentEv.charges.push_back(trkCharge[j]);
        }
        if (currentEv.tracks.size() < 2) continue;
        processedEventsSig++;

        // === MIXING ===
        auto start_mix_lap = std::chrono::high_resolution_clock::now();
        if (!pool[zBin].empty()) {
            processMixDeltaEtaDeltaPhi(thread_count, currentEv, pool[zBin], hMixSS_DeltaEtaDeltaPhi);
            processedEventsMix++;
        }
        auto end_mix_lap = std::chrono::high_resolution_clock::now();
        duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();

        // === SIGNAL ===
        auto start_signal_lap = std::chrono::high_resolution_clock::now();

        processSignalDeltaEtaDeltaPhi(thread_count, currentEv, hSigSS_DeltaEtaDeltaPhi);
        auto end_signal_lap = std::chrono::high_resolution_clock::now();
        duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();

        // === UPDATE POOL ===
        pool[zBin].push_back(currentEv);
        if (pool[zBin].size() > poolSize) {
            pool[zBin].erase(pool[zBin].begin());
        }
    }
    auto end_full = std::chrono::high_resolution_clock::now();
    duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

    std::cout << "#Processed Events(Sig): " << processedEventsSig << std::endl;
    std::cout << "#Processed Events(Mix): " << processedEventsMix << std::endl;

    // Save histograms
    TCanvas *cSigSS_DeltaEtaDeltaPhi = new TCanvas("cSigSS_DeltaEtaDeltaPhi", "SignalSS  #Delta#eta-#Delta#varphi", 1200, 800);
    TCanvas *cMixSS_DeltaEtaDeltaPhi = new TCanvas("cMixSS_DeltaEtaDeltaPhi", "MixedSS  #Delta#eta-#Delta#varphi", 1200, 800);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    // Signal
    cSigSS_DeltaEtaDeltaPhi->cd();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}","PbPb 2.76 TeV | Signal SS | (#Delta#eta, #Delta#varphi)");
    hSigSS_DeltaEtaDeltaPhi->Draw("COLZ");

    // Mixed
    cMixSS_DeltaEtaDeltaPhi->cd();
    drawCMSHeaders("#bf{CMS} #it{Work in Progress}","PbPb 2.76 TeV | Mixed SS | (#Delta#eta, #Delta#varphi)");
    hMixSS_DeltaEtaDeltaPhi->Draw("COLZ");
    

    TH2D *histograms[] = {
        hSigSS_DeltaEtaDeltaPhi, 
        hMixSS_DeltaEtaDeltaPhi
    };
    Int_t numHistograms = 2;

    TCanvas *canvases[] = {
        cSigSS_DeltaEtaDeltaPhi, 
        cMixSS_DeltaEtaDeltaPhi
    };
    Int_t numCanvases = 2;

    TLegend *legends[] = {};
    Int_t numLegends = 0;


    char prefix[256];

    if (ptMax > 0) {
        if (ptMin > 0) {
            sprintf(prefix,
                    "DeltaEtaDeltaPhi_eta-abs%g_pT-%g-to-%g_%s_%g-%g",
                    etaCut, ptMin, ptMax,
                    selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        } else {
            sprintf(prefix,
                    "DeltaEtaDeltaPhi_eta-abs%g_pT-to-%g_%s_%g-%g",
                    etaCut, ptMax,
                    selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        }
    } else {
        if (ptMin > 0) {
            sprintf(prefix,
                    "DeltaEtaDeltaPhi_eta-abs%g_pT-from-%g_%s_%g-%g",
                    etaCut, ptMin,
                    selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        } else {
            sprintf(prefix,
                    "DeltaEtaDeltaPhi_eta-abs%g_%s_%g-%g",
                    etaCut,
                    selectionVarName, displaySelVarMoreeq, displaySelVarLess);
        }
    }

    const char *spath = "benchmarks";
    std::vector<double> durations = {duration_full, duration_signal, duration_mix};
    std::vector<std::string> labels = {"Total Time", "Signal Time", "Mix Time"};
    save_benchmark_chrono(durations, labels, spath, prefix,
                          processedEventsSig, processedEventsMix);

    const char *hpath = "./data/signal_mix/";
    save_histograms2d(histograms, numHistograms, hpath,
                    prefix, displaySelVarMoreeq, displaySelVarLess);

    const char *ipath = "./imgs/test/signal_mix/";
    save_canvas_images(canvases, numCanvases, ipath, prefix, "png");
    save_canvas_images(canvases, numCanvases, ipath, prefix, "pdf");

    AnalysisLog::instance().save("./logs", "sig_DeltaEtaDeltaPhi");

    close_program(canvases, numCanvases,
                  histograms, numHistograms,
                  legends, numLegends, fr);

    }

int main(){
    ROOT::EnableImplicitMT();
    // To compile use:
    // g++ -std=c++17 -pthread makeSignal.cpp -o makeSignal `root-config --cflags --libs` -lGenVector
    // To run use:
    // ./makeSignal
    // --- Analysis Configuration ---
    const char* inputFile = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";  
    const char* treeName  = "demo/TreeMBUCC";

    int test_limit_sig = -1; // -1 for no limit
    int test_limit_mix = -1; // -1 for no limit

    Double_t ptMin = 0.5; // Minimum pT cut for tracks
    Double_t ptMax = -1.0; // Maximum pT cut for tracks (-1 for no max cut)

    Double_t etaCut = 0.95;

    std::vector<double> centralityBins = {3200, 3300}; 
    ControlVar selectedControlVar = ControlVar::CENTHF;
    size_t repeats = 1; 
    
    std::cout << "Starting analysis for " << centralityBins.size() - 1 << " bin(s)." << std::endl;
    for (size_t j = 0; j < repeats; j++){
        for (size_t i = 0; i < centralityBins.size() - 1; ++i){
            double bin_low = centralityBins[i];
            double bin_high = centralityBins[i+1];
            
            std::cout << "\n--- ["<< j << "]Running for centrality bin: " << bin_low << "% - " << bin_high << "% ---" << std::endl;
            sig_DeltaEtaDeltaPhi(inputFile, treeName, bin_low, bin_high, selectedControlVar, ptMin, ptMax, etaCut);
        }
    }
    std::cout << "\nAnalysis finished." << std::endl;

    return 0;
}