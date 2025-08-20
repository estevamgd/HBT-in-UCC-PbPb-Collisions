#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
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
void signal_mix(const char *fileInput, const char *treeInput, double selVarMoreeq, double selVarLess, 
    ControlVar selectionVarType = ControlVar::CENT, Float_t vertexDistance = 2, Int_t pairChargeMult = 1, 
    int poolSizeInt = 10) {
        
    ROOT::EnableImplicitMT();
    auto threadPoolSize = ROOT::GetThreadPoolSize();
    std::cout << "Thread Pool size = " << threadPoolSize << std::endl;

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
    int processedEvents = 0;

    // === HISTOGRAMA ===
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;
    //double nnscale = numBins(ninterval, nlength, nscale), x0 = 0., xt=10.;
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
    
    h_qinvSS_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_signal_1l->GetYaxis()->SetTitle("#Pairs");
 
    h_qinvOS_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOS_signal_1l->GetYaxis()->SetTitle("#Pairs");

    h_qinvSSCor_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_signal_1l->GetYaxis()->SetTitle("#Pairs");

    h_qinvOSCor_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOSCor_signal_1l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinv_mix_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinv_mix_1l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvCor_mix_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvCor_mix_1l->GetYaxis()->SetTitle("#Pairs");

    // Benchmarking
    TStopwatch stopwatchFull, stopwatchMix, stopwatchSignal;
    
    const int numSW = 3;
    TStopwatch* stopWatches[numSW] = {&stopwatchFull, &stopwatchMix, &stopwatchSignal};

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << displaySelVarMoreeq << " to " << displaySelVarLess << std::endl;
    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> vectorEventTracks4V;
    std::vector<std::vector<Int_t>> vectorEventTracksCharge;
    std::vector<std::vector<Float_t>> vectorEventTrackspvZ;
    std::vector<std::vector<Float_t>> vectorEventTracksWeight;

    size_t poolSize = static_cast<size_t>(poolSizeInt);

    stopwatchFull.Start(kFALSE);
    for (Long64_t i = 0; i < nentries; i++){
        t->GetEntry(i);
        
        // update proxies
        hiBinProxy   = hiBin;
        NtrkProxy    = Ntrk;
        HFsumETProxy = HFsumET;
        //std::cout << "pvZ " << fabs(pvZ) << " hibin: " << hiBin << " Ntrk: " << Ntrk << " HFsumET: " << HFsumET << " selectionVar: " <<  *selectionVar << std::endl;
        
        if (!(selVarMoreeq <= *selectionVar && *selectionVar < selVarLess)) continue; 
        processedEvents++;                                                          
        
        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        std::vector<Int_t> currentEventTracksCharge;
        std::vector<Float_t> currentEventTrackspvZ;
        std::vector<Float_t> currentEventTracksWeight;
        
        for (int j = 0; j < Ntrk; j++) {
            // Build Track'ij
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V;
            // The line below makes a 4-vector so we can calculate the q_inv
            currentTrack4V = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            
            // Build Event i
            currentEventTracks4V.push_back(currentTrack4V);
            currentEventTracksCharge.push_back(trkCharge[j]);
            currentEventTrackspvZ.push_back(pvZ);
            currentEventTracksWeight.push_back(trkWeight[j]);
        }
        
        stopwatchMix.Start(kFALSE);
        // === MIXING === 
        if (!currentEventTracks4V.empty()) { //Check if there is any track in the event and if there is any event to mix
            // Acutall mix
            if (vectorEventTracks4V.size() < poolSize  && !vectorEventTracks4V.empty()) { //Check if there is any event to mix and i
                for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); nEv++) {
                    for (size_t k = 0; k < vectorEventTracks4V[nEv].size() * currentEventTracks4V.size(); k++) {
                        size_t i = k / currentEventTracks4V.size();
                        size_t j = k % currentEventTracks4V.size();
                        
                        if (!(fabs(pvZ - vectorEventTrackspvZ[nEv][i]) < vertexDistance)) continue; // Check if their distance respects our setted vertexDistance
                        if (!(vectorEventTracksCharge[nEv][i]*trkCharge[j] == pairChargeMult)) continue; //Check if the pair has same or opposite charge, depending on the configuration
                        if (std::isinf(vectorEventTracksWeight[nEv][i]*trkWeight[j])) continue; // Check if weight is infinity
                        
                        double qinv = GetQ(vectorEventTracks4V[nEv][i], currentEventTracks4V[j]);
                        h_qinv_mix_1l->Fill(qinv);
                        h_qinvCor_mix_1l->Fill(qinv, vectorEventTracksWeight[nEv][i]*trkWeight[j]);
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
        stopwatchMix.Stop();
        
        stopwatchSignal.Start(kFALSE);
        // === SIGNAL ===
        if (currentEventTracks4V.size() > 1) { // check if event has 2 or more tracks 
            // Single Loop
            for (size_t k = 0; k < currentEventTracks4V.size()*(currentEventTracks4V.size() - 1); k++) {
                size_t p1 = k / (currentEventTracks4V.size() - 1);
                size_t p2 = (k % (currentEventTracks4V.size() - 1)) + 1;

                // Checks 
                if (p1 >= p2) continue; // This makes sure we dont repeat a pair
                if (std::isinf(trkWeight[p1]*trkWeight[p2])) continue; // Check if weight is infinity
                
                // Build Histogram
                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);
                
                if (trkCharge[p1]*trkCharge[p2] > 0){ // Fills same charge particle pair histogram
                    h_qinvSS_signal_1l->Fill(qinv);
                    h_qinvSSCor_signal_1l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                } else { // Fills opposite charge particle pair histogram
                    h_qinvOS_signal_1l->Fill(qinv);
                    h_qinvOSCor_signal_1l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                }
            }
        }
        stopwatchSignal.Stop();
    }
    stopwatchFull.Stop();

    std::cout << "#Processed Events: " << processedEvents << std::endl;
    
    TH1D *h_qinvDivCor_signal_1l = (TH1D *)h_qinvSSCor_signal_1l->Clone("h_qinvDivCor_signal_1l");
    TH1D *h_qinvDiv_signal_1l = (TH1D *)h_qinvSS_signal_1l->Clone("h_qinvDiv_signal_1l");    

    h_qinvDivCor_signal_1l->Divide(h_qinvOSCor_signal_1l);
    h_qinvDiv_signal_1l->Divide(h_qinvOS_signal_1l);


    h_qinvDiv_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDiv_signal_1l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvDivCor_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvDivCor_signal_1l->GetYaxis()->SetTitle("#Pairs");

    h_qinvDivCor_signal_1l->Sumw2();

    // Save histograms
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

    h_qinvSS_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");
    h_qinvOS_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");
    h_qinvSSCor_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Corrected");
    h_qinvOSCor_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Corrected");
    
    h_qinvDivCor_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division Corrected");
    h_qinvDiv_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division");
    
    h_qinv_mix_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division");
    h_qinvCor_mix_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division Corrected");
    
    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);

    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    
    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_SS->AddEntry(h_qinvSS_signal_1l, "Signal SS - One loop", "l");
    legend_sig_OS->AddEntry(h_qinvOS_signal_1l, "Signal OS - One loop", "l");
    legend_sig_cor_SS->AddEntry(h_qinvSSCor_signal_1l, "Signal SS Cor - One loop", "l");
    legend_sig_cor_OS->AddEntry(h_qinvOSCor_signal_1l, "Signal OS Cor - One loop", "l");
    
    legend_sig_cor_Div->AddEntry(h_qinvOSCor_signal_1l, "Signal Div Cor - One loop", "l");
    legend_sig_Div->AddEntry(h_qinvOSCor_signal_1l, "Signal Div - One loop", "l");
   
    legend_mix->AddEntry(h_qinv_mix_1l, "Mix - One loop", "l");
    legend_mix_cor->AddEntry(h_qinvCor_mix_1l, "Mix Cor - One loop", "l");

    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qinvSS_signal_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qinvOS_signal_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_signal_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_signal_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qinvDivCor_signal_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qinvDiv_signal_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();
    
    c_mix->cd(); gStyle->SetOptStat(0); h_qinv_mix_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix->Draw();
    c_mix_cor->cd(); gStyle->SetOptStat(0); h_qinvCor_mix_1l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_mix_cor->Draw();
    
    Int_t numHistograms = 8;
    TH1D *histograms[] = {h_qinvSS_signal_1l, h_qinvSSCor_signal_1l, h_qinvOS_signal_1l, h_qinvOSCor_signal_1l, h_qinvDivCor_signal_1l, h_qinvDiv_signal_1l, h_qinv_mix_1l, h_qinvCor_mix_1l};

    // Prefix


    char prefix[50];
    sprintf(prefix, "signal_mix_%s_%f-%f", selectionVarName, displaySelVarMoreeq, displaySelVarLess);

    // Saving Benchmakrs
    const char *spath = "benchmarks";
    save_benchmark(stopWatches, numSW, spath, prefix, processedEvents);

    // Saving histograms
    const char *hpath = "./data/sigal_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, displaySelVarMoreeq, displaySelVarLess);

    // Saving image
    const char *ipath = "./imgs/test/sigal_mix/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);
    
    // Closing program
    close_program(canvases, numCanvases, histograms, numHistograms, fr);     
}

void executeSigMix(){
    // --- Analysis Configuration ---
    const char* inputFile = "data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root";  
    const char* treeName  = "demo/TreeMBUCC";

    // Define the centrality bins to analyze. The code will loop over adjacent pairs.
    // This is easy to modify - just add or remove numbers from the list.
    std::vector<double> centralityBins = {50, 70}; 
    // Example for more bins: std::vector<double> centralityBins = {0, 10, 30, 50, 70, 90};

    ControlVar selectedControlVar = ControlVar::CENT;
    Float_t selectedVertexDistance = 2.0f;
    int poolSize = 10;   
    Int_t pairChargeMult = 1; // 1 for Same-Sign, -1 for Opposite-Sign

    // --- Loop over Bins and Execute Analysis ---
    std::cout << "Starting analysis for " << centralityBins.size() - 1 << " bin(s)." << std::endl;
    
    // Loop automatically handles the number of bins. Use size_t to avoid compiler warnings.
    for (size_t i = 0; i < centralityBins.size() - 1; ++i){
        double bin_low = centralityBins[i];
        double bin_high = centralityBins[i+1];
        
        std::cout << "\n--- Running for centrality bin: " << bin_low << "% - " << bin_high << "% ---" << std::endl;
        
        // Pass all the configured parameters to the analysis function
        signal_mix(inputFile, treeName, bin_low, bin_high, 
                   selectedControlVar, selectedVertexDistance, pairChargeMult, poolSize);
    }
    
    std::cout << "\nAnalysis finished." << std::endl;
}