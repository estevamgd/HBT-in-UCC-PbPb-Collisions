#ifndef SIG_MIX_H
#define SIG_MIX_H

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

#include "data_func.h"
#include "my_func.h"
#include "sig_mix.h"

void sig_mix(const char *fileInput, const char *treeInput, 
    int selectionVarI = 100, int selectionVarF = 375, 
    ControlVar choice = ControlVar::CENT, int poolSize = 10) {
    
    /*ROOT::EnableImplicitMT();
    auto threadPoolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << threadPoolSize << std::endl;*/

    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree(fileInput, treeInput, fr, t);

    // Variables from the tree
    Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
    Float_t HFsumET, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize], 
            pionMass = 0.13957039; // Pion mass [GeV] from PDG

    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("trkCharge", trkCharge);
    t->SetBranchAddress("trkWeight", trkWeight);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    
    // Set the type of control var, centrality, multiplicity or vertex
    int choices[nControlVar] = {hiBin, Ntrk, maxSize}; // maxsize holding position for vertex
    const char * choicesNames[nControlVar] = {"cent", "mult", "vertex"}; 
    
    int selectionVar = choices[0];
    const char * selectionVarName = choicesNames[0];

    int intchoice = int(choice);

    if (intchoice <= nControlVar){
        selectionVar = choices[intchoice];
        selectionVarName = choicesNames[intchoice];
    } 

    // Getting how many entries
    Long64_t nentries = t->GetEntries(), trackvec_max_size = 0;

    // Histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;
    
    double nnscale = numBins(ninterval, nlength, nscale), x0 = 0., xt=4;
    
    TH1D* h_qinvSS_signal_1l = new TH1D("h_qinvSS_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_signal_1l = new TH1D("h_qinvSSCor_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_signal_1l = new TH1D("h_qinvOS_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_signal_1l = new TH1D("h_qinvOSCor_signal_1l", "", nnscale, x0, xt);
    TH1D* h_qinvSS_mixing_1l = new TH1D("h_qinvSS_mixing_1l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_mixing_1l = new TH1D("h_qinvSSCor_mixing_1l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_mixing_1l = new TH1D("h_qinvOS_mixing_1l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_mixing_1l = new TH1D("h_qinvOSCor_mixing_1l", "", nnscale, x0, xt);
    TH1D* h_qinvSS_signal_2l = new TH1D("h_qinvSS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_signal_2l = new TH1D("h_qinvSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_signal_2l = new TH1D("h_qinvOS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_signal_2l = new TH1D("h_qinvOSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvSS_mixing_2l = new TH1D("h_qinvSS_mixing_2l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_mixing_2l = new TH1D("h_qinvSSCor_mixing_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_mixing_2l = new TH1D("h_qinvOS_mixing_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_mixing_2l = new TH1D("h_qinvOSCor_mixing_2l", "", nnscale, x0, xt);
        
    h_qinvSSCor_signal_1l->Sumw2();
    h_qinvOSCor_signal_1l->Sumw2();
    h_qinvSSCor_mixing_1l->Sumw2();
    h_qinvOSCor_mixing_1l->Sumw2();
    h_qinvSSCor_signal_2l->Sumw2();
    h_qinvOSCor_signal_2l->Sumw2();
    h_qinvSSCor_mixing_2l->Sumw2();
    h_qinvOSCor_mixing_2l->Sumw2();
    
    h_qinvSS_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_signal_1l->GetYaxis()->SetTitle("#Pairs");

    h_qinvSS_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvSS_mixing_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_mixing_1l->GetYaxis()->SetTitle("#Pairs");

    h_qinvSS_mixing_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_mixing_2l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvSSCor_signal_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_signal_1l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvSSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvSSCor_mixing_1l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_mixing_1l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvSSCor_mixing_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_mixing_2l->GetYaxis()->SetTitle("#Pairs");
    
    Int_t numHistograms = 16;
    TH1D *histograms[] = {h_qinvSS_signal_1l, h_qinvSSCor_signal_1l, h_qinvOS_signal_1l, h_qinvOSCor_signal_1l, 
        h_qinvSS_mixing_1l, h_qinvSSCor_mixing_1l, h_qinvOS_mixing_1l, h_qinvOSCor_mixing_1l, 
        h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l, 
        h_qinvSSCor_mixing_2l, h_qinvSS_mixing_2l, h_qinvOSCor_mixing_2l, h_qinvOS_mixing_2l};
        
    // Benchmarking
    TStopwatch stopwatchFull, stopwatchMix, stopwatchMixPhi, stopwatchSignal, stopwatchSingle, stopwatchSingleMix, 
        stopwatchSingleMixPhi, stopwatchSingleSignal,  stopwatchDouble, stopwatchDoubleMix, stopwatchDoubleMixPhi, 
        stopwatchDoubleSignal;

    const int numSW = 12;
    TStopwatch* stopWatches[numSW] = {&stopwatchFull, &stopwatchMix, &stopwatchMixPhi, &stopwatchSignal, &stopwatchSingle, &stopwatchSingleMix, 
        &stopwatchSingleMixPhi, &stopwatchSingleSignal,  &stopwatchDouble, &stopwatchDoubleMix, &stopwatchDoubleMixPhi, 
        &stopwatchDoubleSignal};
    
    // --- POOL SETTINGS ---
    std::vector<std::vector<TrackInfo>> poolEvents;
    std::vector<std::vector<TrackInfo>> phiEvents;

    int poolMaxTracks = 0;
    int phiMaxTracks = 0;
    int phiSize = poolSize - (nentries % poolSize);
    bool phiSwitch = true;
    if (phiSize < 1) phiSwitch = false;
        
    int phiSteps;
    
    if (phiSwitch) {
        phiSteps = (nentries - (poolSize - phiSize)) / phiSize; 
    } else {
        phiSteps = 1;
    }
    
    int x = 0;
    // ------

    stopwatchFull.Start(kFALSE);
    // Create 4-vector
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
 
        std::vector<TrackInfo> currentEventTracks;
        for (int j = 0; j < Ntrk; j++) {
            // Build Track'ij
            TrackInfo currentTrack;
            currentTrack.p4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            currentTrack.charge = trkCharge[j];
            currentTrack.weight = trkWeight[j];
            
            // === MIXING ===
            stopwatchMix.Start(kFALSE);

            if (poolEvents.size() == poolSize - 1) {
                // Double Loop
                stopwatchDouble.Start(kFALSE);
                stopwatchDoubleMix.Start(kFALSE);

                for (int k = 0; k < poolEvents.size(); k++) {
                    for (int z = 0; z < poolEvents[k].size(); z++) {
                        // Checks
                        if (selectionVar < selectionVarI || selectionVar > selectionVarF) continue;
                        
                        // Build Histogram
                        double qinv = GetQ(currentTrack.p4, poolEvents[k][z].p4);
                        if (trkCharge[j]*poolEvents[k][z].charge>0){
                            h_qinvSS_mixing_2l->Fill(qinv);
                            h_qinvSSCor_mixing_2l->Fill(qinv, trkWeight[j]*poolEvents[k][z].weight);
                        } else {
                            h_qinvOS_mixing_2l->Fill(qinv);
                            h_qinvOSCor_mixing_2l->Fill(qinv, trkWeight[j]*poolEvents[k][z].weight);
                        }
                    }
                }
                stopwatchDouble.Stop();
                stopwatchDoubleMix.Stop();

                // Single Loop
                stopwatchSingle.Start(kFALSE);
                stopwatchSingleMix.Start(kFALSE);

                for (int w = 0; w < (poolMaxTracks)*(poolEvents.size()); w++) { 
                    int k = w / poolMaxTracks; 
                    int z = w % poolMaxTracks; 

                    // Checks
                    if (z > poolEvents[k].size() - 1 || poolEvents[k].size() == 0) continue;
                    if (selectionVar < selectionVarI || selectionVar > selectionVarF) continue;
                    
                    // Build Histogram
                    double qinv = GetQ(currentTrack.p4, poolEvents[k][z].p4);
                    if (trkCharge[j]*poolEvents[k][z].charge>0){
                        h_qinvSS_mixing_1l->Fill(qinv);
                        h_qinvSSCor_mixing_1l->Fill(qinv, trkWeight[j]*poolEvents[k][z].weight);
                    } else {
                        h_qinvOS_mixing_1l->Fill(qinv);
                        h_qinvOSCor_mixing_1l->Fill(qinv, trkWeight[j]*poolEvents[k][z].weight);
                    }
                }
                stopwatchSingle.Stop();
                stopwatchSingleMix.Stop();
                
            }
            stopwatchMix.Stop();

            // Build Event i
            currentEventTracks.push_back(currentTrack);
            
        }

        // --- Phi Handler ---
        if (x == phiSteps && phiSwitch && phiEvents.size() < phiSize) {
            phiEvents.push_back(currentEventTracks);
            if (currentEventTracks.size() > phiMaxTracks) phiMaxTracks = currentEventTracks.size();
            x = 0;
        }

        // --- Pool Handler ---
        poolEvents.push_back(currentEventTracks);
        if (currentEventTracks.size() > poolMaxTracks) poolMaxTracks = currentEventTracks.size();

        if (poolEvents.size() == poolSize) {
            poolEvents.clear();
            poolMaxTracks = 0;
        }

        // === SIGNAL ===
        stopwatchSignal.Start(kFALSE);

        if (currentEventTracks.size() > 1) {
            // Double Loop
            stopwatchDouble.Start(kFALSE);
            stopwatchDoubleSignal.Start(kFALSE);

            for (size_t p1 = 0; p1 < currentEventTracks.size(); p1++) {
                for (size_t p2 = p1 + 1; p2 < currentEventTracks.size(); p2++) {
                    // Build Histogram
                    double qinv = GetQ(currentEventTracks[p1].p4, currentEventTracks[p2].p4);
                    if (trkCharge[p1]*trkCharge[p2]>0){
                        h_qinvSS_signal_2l->Fill(qinv);
                        h_qinvSSCor_signal_2l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                    } else {
                        h_qinvOS_signal_2l->Fill(qinv);
                        h_qinvOSCor_signal_2l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                    }
                }
            }
            stopwatchDouble.Stop();
            stopwatchDoubleSignal.Stop();

            // Single Loop
            stopwatchSingle.Start(kFALSE);
            stopwatchSingleSignal.Start(kFALSE);
            /*
            int s = currentEventTracks.size();
            int wt = (s * (s - 1)) / 2;

            for (int w = 0; w < wt; w++) {
                int p1 = s - 2 - (int)(std::sqrt(-8 * w + 4 * s * (s - 1) - 7) / 2.0 - 0.5);
                int t = (p1 * (2 * s - p1 - 1)) / 2;
                int p2 = w - t + p1 + 1;
                
                // Build Histogram
                double qinv = GetQ(currentEventTracks[p1].p4, currentEventTracks[p2].p4);
                if (trkCharge[p1]*trkCharge[p2]>0){
                    h_qinvSS_signal_1l->Fill(qinv);
                    h_qinvSSCor_signal_1l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                } else {
                    h_qinvOS_signal_1l->Fill(qinv);
                    h_qinvOSCor_signal_1l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                }
            }
            */
            for (int k = 0; k < currentEventTracks.size()*(currentEventTracks.size() - 1); k++) {
                int p1 = k / (currentEventTracks.size() - 1);
                int p2 = (k % (currentEventTracks.size() - 1)) + 1;

                // Checks 
                if (p1 >= p2) continue;
                
                // Build Histogram
                double qinv = GetQ(currentEventTracks[p1].p4, currentEventTracks[p2].p4);
                if (trkCharge[p1]*trkCharge[p2]>0){
                    h_qinvSS_signal_1l->Fill(qinv);
                    h_qinvSSCor_signal_1l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                } else {
                    h_qinvOS_signal_1l->Fill(qinv);
                    h_qinvOSCor_signal_1l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                }
            }
            stopwatchSingle.Stop();
            stopwatchSingleSignal.Stop();
        }
        stopwatchSignal.Stop();

        x++;
    }
    
    // === PHI MIXING ===
    stopwatchMixPhi.Start(kFALSE);

    if (phiSwitch) {
        for (int event = 0; event < poolEvents.size(); event++) {
            phiEvents.push_back(poolEvents[event]);
        }

        if (poolMaxTracks > phiMaxTracks) phiMaxTracks = poolMaxTracks;

        int lastIndex = phiEvents.size() - 1; if (phiEvents.size() == 0) lastIndex = 0;
        int lastEventSize = phiEvents[lastIndex].size();

        for (int lastEventTrackID = 0; lastEventTrackID < lastEventSize; lastEventTrackID++) {
            // Double loop
            stopwatchDouble.Start(kFALSE);
            stopwatchDoubleMixPhi.Start(kFALSE);
            for (int k = 0; k < phiEvents.size() - 1; k++) {
                for (int z = 0; z < phiEvents[k].size(); z++) {
                    // Checks 
                    if (selectionVar < selectionVarI || selectionVar > selectionVarF) continue;
                
                    // Building Histograms
                    double qinv = GetQ(phiEvents[lastIndex][lastEventTrackID].p4, phiEvents[k][z].p4);
                    if (phiEvents[lastIndex][lastEventTrackID].charge*phiEvents[k][z].charge > 0) {
                        h_qinvSS_mixing_2l->Fill(qinv);
                        h_qinvSSCor_mixing_2l->Fill(qinv, phiEvents[lastIndex][lastEventTrackID].weight*phiEvents[k][z].weight);
                    } else {
                        h_qinvOS_mixing_2l->Fill(qinv);
                        h_qinvOSCor_mixing_2l->Fill(qinv, phiEvents[lastIndex][lastEventTrackID].weight*phiEvents[k][z].weight);
                    }
                    
                }
            }
            stopwatchDouble.Stop();
            stopwatchDoubleMixPhi.Stop();

            // Single loop
            stopwatchSingle.Start(kFALSE);
            stopwatchSingleMixPhi.Start(kFALSE);
            for (int w = 0; w < phiMaxTracks*(phiEvents.size()-1); w++) {
                int k = w / phiMaxTracks;
                int z = w % phiMaxTracks;
                
                // Checks
                if (z > poolEvents[k].size() - 1 || poolEvents[k].size() == 0) continue;
                if (selectionVar < selectionVarI || selectionVar > selectionVarF) continue;

                // Building Histograms
                double qinv = GetQ(phiEvents[lastIndex][lastEventTrackID].p4, phiEvents[k][z].p4);
                if (phiEvents[lastIndex][lastEventTrackID].charge*phiEvents[k][z].charge > 0) {
                    h_qinvSS_mixing_1l->Fill(qinv);
                    h_qinvSSCor_mixing_1l->Fill(qinv, phiEvents[lastIndex][lastEventTrackID].weight*phiEvents[k][z].weight);
                } else {
                    h_qinvOS_mixing_1l->Fill(qinv);
                    h_qinvOSCor_mixing_1l->Fill(qinv, phiEvents[lastIndex][lastEventTrackID].weight*phiEvents[k][z].weight);
                }
            }
            stopwatchSingle.Stop();
            stopwatchSingleMixPhi.Stop();
        }
    }
    stopwatchMixPhi.Stop();
    stopwatchFull.Stop();
    
    // Save histograms
    TCanvas *c_sig = new TCanvas("c_sig", "Signal Distributions", 1200, 800);
    TCanvas *c_sig_cor = new TCanvas("c_sig_cor", "Signal Correlation Distributions", 1200, 800);
    TCanvas *c_mix = new TCanvas("c_mix", "Mixing Distributions", 1200, 800);
    TCanvas *c_mix_cor = new TCanvas("c_mix_cor", "Mixing Correlation Distributions", 1200, 800);
    
    TCanvas *canvases[] = {c_sig, c_sig_cor, c_mix, c_mix_cor};
    int numCanvases = 4;

    c_sig->cd();
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    h_qinvSS_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");

    h_qinvSS_signal_1l->Draw("HIST PLC");
    h_qinvOS_signal_1l->Draw("HIST SAME PLC");
    h_qinvSS_signal_2l->Draw("HIST SAME PLC");
    h_qinvOS_signal_2l->Draw("HIST SAME PLC");

    TLegend *legend_sig = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig->AddEntry(h_qinvSS_signal_1l, "Signal SS - One loop", "l");
    legend_sig->AddEntry(h_qinvOS_signal_1l, "Signal OS - One loop", "l");
    legend_sig->AddEntry(h_qinvSS_signal_2l, "Signal SS - Two loops", "l");
    legend_sig->AddEntry(h_qinvOS_signal_2l, "Signal OS - Two loops", "l");
    legend_sig->Draw();
    
    c_sig->Update();
    c_sig_cor->cd();
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    h_qinvSS_signal_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Correlation");

    h_qinvSSCor_signal_1l->Draw("HIST PLC");
    h_qinvOSCor_signal_1l->Draw("HIST SAME PLC");    
    h_qinvSSCor_signal_2l->Draw("HIST SAME PLC");
    h_qinvOSCor_signal_2l->Draw("HIST SAME PLC");

    TLegend *legend_sig_cor = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_sig_cor->AddEntry(h_qinvSSCor_signal_1l, "Signal SS Cor - One loop", "l");
    legend_sig_cor->AddEntry(h_qinvOSCor_signal_1l, "Signal OS Cor- One loop", "l");
    legend_sig_cor->AddEntry(h_qinvSSCor_signal_2l, "Signal SS Cor - Two loops", "l");
    legend_sig_cor->AddEntry(h_qinvOSCor_signal_2l, "Signal OS Cor- Two loops", "l");
    legend_sig_cor->Draw();

    c_sig_cor->Update();
    c_mix->cd();
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    h_qinvSS_mixing_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Mixing");

    h_qinvSS_mixing_1l->Draw("HIST PLC");
    h_qinvOS_mixing_1l->Draw("HIST SAME PLC");
    h_qinvSS_mixing_2l->Draw("HIST SAME PLC");
    h_qinvOS_mixing_2l->Draw("HIST SAME PLC");

    TLegend *legend_mix = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_mix->AddEntry(h_qinvSS_mixing_1l, "Mixing SS - One loop", "l");
    legend_mix->AddEntry(h_qinvOS_mixing_1l, "Mixing OS - One loop", "l");
    legend_mix->AddEntry(h_qinvSS_mixing_2l, "Mixing SS - Two loops", "l");
    legend_mix->AddEntry(h_qinvOS_mixing_2l, "Mixing OS - Two loops", "l");
    legend_mix->Draw();
    
    c_mix->Update();
    c_mix_cor->cd();
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLake);

    h_qinvSS_mixing_1l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Mixing Correlation");

    h_qinvSSCor_mixing_1l->Draw("HIST PLC");
    h_qinvOSCor_mixing_1l->Draw("HIST SAME PLC");    
    h_qinvSSCor_mixing_2l->Draw("HIST SAME PLC");
    h_qinvOSCor_mixing_2l->Draw("HIST SAME PLC");

    TLegend *legend_mix_cor = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_mix_cor->AddEntry(h_qinvSSCor_mixing_1l, "Mixing SS Cor - One loop", "l");
    legend_mix_cor->AddEntry(h_qinvOSCor_mixing_1l, "Mixing OS Cor- One loop", "l");
    legend_mix_cor->AddEntry(h_qinvSSCor_mixing_2l, "Mixing SS Cor - Two loops", "l");
    legend_mix_cor->AddEntry(h_qinvOSCor_mixing_2l, "Mixing OS Cor- Two loops", "l");
    legend_mix_cor->Draw();

    c_mix_cor->Update();
    
    // Prefix
    char prefix[30];
    sprintf(prefix, "sigmix_%s", selectionVarName);

    // Saving Benchmakrs
    const char *spath = "benchmarks";
    save_benchmark(stopWatches, numSW, spath, prefix);

    // Saving histograms
    const char *hpath = "./data/Sig_mix/";
    save_histograms(histograms, numHistograms, hpath, prefix, selectionVarI, selectionVarF);

    // Saving image
    const char *ipath = "./imgs/test/";
    const char *ifile_type = "png";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type);
    const char *ifile_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, ipath, prefix, ifile_type2);
    
    // Closing program
    close_program(canvases, numCanvases, histograms, numHistograms, fr);     
}

#endif 