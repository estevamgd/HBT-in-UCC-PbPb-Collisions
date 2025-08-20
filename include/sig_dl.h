#ifndef SIG_DL_H
#define SIG_DL_H

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
#include "sig_dl.h"

void sig_dl(const char *fileInput, const char *treeInput, int selectionVarI, int selectionVarF) {
    ROOT::EnableImplicitMT();
    auto threadPoolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << threadPoolSize << std::endl;

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
    
    const char * selectionVarName = "cent";

    // Getting how many entries
    Long64_t nentries = t->GetEntries(), trackvec_max_size = 0;
    std::cout << "#Events: " << nentries << std::endl;
    
    // See how many processed events
    int processedEvents = 0;

    // Histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;
    
    //double nnscale = numBins(ninterval, nlength, nscale), x0 = 0., xt=10.;
    double nnscale = 10000, x0 = 0., xt=10.;
    
    TH1D* h_qinvSS_signal_2l = new TH1D("h_qinvSS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvSSCor_signal_2l = new TH1D("h_qinvSSCor_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOS_signal_2l = new TH1D("h_qinvOS_signal_2l", "", nnscale, x0, xt);
    TH1D* h_qinvOSCor_signal_2l = new TH1D("h_qinvOSCor_signal_2l", "", nnscale, x0, xt);
        
    h_qinvSSCor_signal_2l->Sumw2();
    h_qinvOSCor_signal_2l->Sumw2();
    
    h_qinvSS_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSS_signal_2l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvOS_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOS_signal_2l->GetYaxis()->SetTitle("#Pairs");

    h_qinvSSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvSSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");
    
    h_qinvOSCor_signal_2l->GetXaxis()->SetTitle("q_{inv}[GeV]");
    h_qinvOSCor_signal_2l->GetYaxis()->SetTitle("#Pairs");

        
    // Benchmarking
    TStopwatch stopwatchFull;

    const int numSW = 1;
    TStopwatch* stopWatches[numSW] = {&stopwatchFull};

    std::cout << "Processing " << fileInput << "/" << treeInput << " events with centrality from " << selectionVarI << " to " << selectionVarF << std::endl;

    stopwatchFull.Start(kFALSE);
    // Create 4-vector
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        
        if (selectionVarI > hiBin || hiBin > selectionVarF) continue;
        processedEvents++;

        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> currentEventTracks4V;
        for (int j = 0; j < Ntrk; j++) {
            // Build Track'ij
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> currentTrack4V;
            // The line below makes a 4-vector so we can calculate the q_inv
            currentTrack4V = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass);
            
            // Build Event i
            currentEventTracks4V.push_back(currentTrack4V);
        }

        // === SIGNAL ===
        if (currentEventTracks4V.size() > 1) {
            // Double Loop
            for (size_t p1 = 0; p1 < currentEventTracks4V.size(); p1++) {
                for (size_t p2 = p1 + 1; p2 < currentEventTracks4V.size(); p2++) {
                    //Checks
                    if (std::isinf(trkWeight[p1]*trkWeight[p2])) continue; // Check if weight is infinity
                    
                    // Build Histogram
                    double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);
                    
                    if (trkCharge[p1]*trkCharge[p2] > 0){ // Fills same charge particle pair histogram
                        h_qinvSS_signal_2l->Fill(qinv);
                        h_qinvSSCor_signal_2l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                    } else { // Fills opposite  charge particle pair histogram
                        h_qinvOS_signal_2l->Fill(qinv);
                        h_qinvOSCor_signal_2l->Fill(qinv, trkWeight[p1]*trkWeight[p2]);
                    }
                }
            }
        }

        currentEventTracks4V.clear();
    }
    stopwatchFull.Stop();

    std::cout << "#Processed Events: " << processedEvents << std::endl;
    
    TH1D *h_qinvDivCor_signal_2l = (TH1D *)h_qinvSSCor_signal_2l->Clone("h_qinvDivCor_signal_1l");
    TH1D *h_qinvDiv_signal_2l = (TH1D *)h_qinvSS_signal_2l->Clone("h_qinvDiv_signal_1l");    

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
    
    TCanvas *canvases[] = {c_sig_SS, c_sig_OS, c_sig_cor_OS, c_sig_cor_SS, c_sig_cor_Div, c_sig_Div};
    int numCanvases = 6;

    h_qinvSS_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");
    h_qinvOS_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal");
    h_qinvSSCor_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Correlation");
    h_qinvOSCor_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Correlation");
    
    h_qinvDivCor_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division Correlation");
    h_qinvDiv_signal_2l->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV | Signal Division");
    
    TLegend *legend_sig_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_OS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_SS = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_cor_OS = new TLegend(0.7, 0.7, 0.9, 0.9);

    TLegend *legend_sig_cor_Div = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend_sig_Div = new TLegend(0.7, 0.7, 0.9, 0.9);

    legend_sig_SS->AddEntry(h_qinvSS_signal_2l, "Signal SS - One loop", "l");
    legend_sig_OS->AddEntry(h_qinvOS_signal_2l, "Signal OS - One loop", "l");
    legend_sig_cor_SS->AddEntry(h_qinvSSCor_signal_2l, "Signal SS Cor - One loop", "l");
    legend_sig_cor_OS->AddEntry(h_qinvOSCor_signal_2l, "Signal OS Cor - One loop", "l");
    
    legend_sig_cor_Div->AddEntry(h_qinvOSCor_signal_2l, "Signal Div Cor - One loop", "l");
    legend_sig_Div->AddEntry(h_qinvOSCor_signal_2l, "Signal Div - One loop", "l");

    c_sig_SS->cd(); gStyle->SetOptStat(0); h_qinvSS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_SS->Draw();
    c_sig_OS->cd(); gStyle->SetOptStat(0); h_qinvOS_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_OS->Draw();
    c_sig_cor_SS->cd(); gStyle->SetOptStat(0); h_qinvSSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_SS->Draw();
    c_sig_cor_OS->cd(); gStyle->SetOptStat(0); h_qinvOSCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_OS->Draw();
    
    c_sig_cor_Div->cd(); gStyle->SetOptStat(0); h_qinvDivCor_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_cor_Div->Draw();
    c_sig_Div->cd(); gStyle->SetOptStat(0); h_qinvDiv_signal_2l->Draw("HIST PLC"); gStyle->SetPalette(kLake); legend_sig_Div->Draw();
    
    Int_t numHistograms = 6;
    TH1D *histograms[] = {h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l, h_qinvDivCor_signal_2l, h_qinvDiv_signal_2l};
    
    // Prefix
    char prefix[50]; 
    sprintf(prefix, "sig_dl_%s_%d-%d", selectionVarName, selectionVarI, selectionVarF);

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