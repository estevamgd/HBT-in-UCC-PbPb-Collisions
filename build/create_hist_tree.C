#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TH1D.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <TStyle.h>
#include "TLegend.h"
#include <TText.h>
#include <TBenchmark.h>
#include <TProfile.h>
#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"


void create_hist_tree() {
    ROOT::EnableImplicitMT();
    gBenchmark->Start("create_hist_tree");
    auto poolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << poolSize << std::endl;

    // Setting up, copy from HERE to...
    // Getting TTree
    TFile *fr = nullptr;
    TTree *t = nullptr;

    getFileTree("data/HiForestAOD_UPC.root", "demo/HBT", fr, t);
    
    // Variables
    Int_t maxSize = 17100;
    Int_t Ntrk = maxSize, NSSpair = maxSize, NOSpair = maxSize;
    Float_t HFsumET;
    Float_t coulombWOS[NOSpair], coulombWSS[NSSpair], qinvSigSS[NSSpair], qinvSigOS[NOSpair], 
            trkPt[Ntrk], trkEta[Ntrk], trkPhi[Ntrk], trkPtRes[Ntrk], 
            trkDxySig[Ntrk], trkNpixLayers[Ntrk], trkDzSig[Ntrk];
    
    // Arrays of variables
    void* variables[] = {
        &HFsumET, &Ntrk, &NSSpair, &NOSpair, coulombWOS, coulombWSS, qinvSigOS, qinvSigSS, 
        trkPt, trkEta, trkPhi, trkPtRes, trkDxySig, trkDzSig, trkNpixLayers
    };

    // Branch names
    const char* branchNames[] = {
        "HFsumET", "Ntrk", "NSSpair", "NOSpair", "coulombWOS", "coulombWSS", 
        "qinvSigOS", "qinvSigSS", "trkPt", "trkEta", "trkPhi", "trkPtRes", 
        "trkDxySig", "trkDzSig", "trkNpixLayers"
    };

    // Setting addresses
    int numBranches = sizeof(variables) / sizeof(variables[0]);
    for (int i = 0; i < numBranches; i++) {
        t->SetBranchAddress(branchNames[i], variables[i]);
    }

    // ...HERE
    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    double d_nentries = t->GetEntries();

    // Setting canvases
    TCanvas *c1 = new TCanvas("c1", "", 1920, 1080);

    TCanvas *canvases[] = {c1};
    int numCanvases = 1;

    // Setting histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;

    TH1D *h1 = cHist("h1", "", "",  0., 1., ninterval, nlength, nscale, 0, 0, 0, 932, 1, 1);
    TH1D *h2 = cHist("h2", "", "",  0., 1., ninterval, nlength, nscale, 0, 0, 0, 932, 1, 1);
    TH1D *h3 = cHist("h3", "", "",  0., 1., ninterval, nlength, nscale, 0, 0, 0, 932, 1, 1);
    TH1D *h4 = cHist("h4", "", "",  0., 1., ninterval, nlength, nscale, 0, 0, 0, 932, 1, 1);

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        if (HFsumET > 100 && HFsumET < 375) { // This selects a centrality
            for (int k = 0; k < NOSpair; k++) {
                h1->Fill(qinvSigOS[k]);
                h3->Fill(qinvSigOS[k], coulombWOS[k]);

            }
            for (int l = 0; l < NSSpair; l++) {
                h2->Fill(qinvSigSS[l]);
                h4->Fill(qinvSigSS[l], coulombWSS[l]);
            }    
        }

        // Display progress
        float progress = (float)(i + 1) / nentries * 100;
        std::cout << "\rFilling Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }

    std::cout << std::endl;

    TH1D *h1_normalized = (TH1D *)h1->Clone("h1_normalized");
    TH1D *h2_normalized = (TH1D *)h2->Clone("h2_normalized");
    TH1D *h3_normalized = (TH1D *)h3->Clone("h3_normalized");
    TH1D *h4_normalized = (TH1D *)h4->Clone("h4_normalized");
    
    TH1D *tonormhist[] = {h1_normalized, h2_normalized, h3_normalized, h4_normalized};
    int numTonorm = 4;

    // Normalize histograms
    Double_t scale = 1;
    normalizer(tonormhist, numTonorm, scale);
   
    // Dividing SS/OS
    TH1D *sr = (TH1D *)h2_normalized->Clone("sr");
    TH1D *sr_cor = (TH1D *)h4_normalized->Clone("sr_cor");
    
    sr->Divide(h1_normalized);
    sr->SetLineColor(kRed);

    sr_cor->Divide(h3_normalized);
    sr_cor->SetLineColor(kBlue);

    // Saving into a Tree
    TFile save_f2("data/50_70_qinv_normqinv_sr.root","RECREATE");
    /*    auto save_t = new TTree("save_t", "some_hist");
    
    save_t->Branch("h1", "TH1D", &h1);
    save_t->Branch("h2", "TH1D", &h2);
    save_t->Branch("h1_normalized", "TH1D", &h1_normalized);
    save_t->Branch("h2_normalized", "TH1D", &h2_normalized);
    save_t->Branch("sr", "TH1D", &sr);

    save_t->Fill();
    save_t->Write();*/

    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h1_normalized->Write();
    h2_normalized->Write();
    h3_normalized->Write();
    h4_normalized->Write();
    sr->Write();
    sr_cor->Write();

    save_f2.Close();
    
    // Closing program
    TH1D *histograms[] = {h1, h2, h3, h4, h1_normalized, h2_normalized, h3_normalized, h4_normalized, sr, sr_cor};
    int numHistograms = 10;

    close_program(canvases, numCanvases, histograms, numHistograms);
    gBenchmark->Show("create_hist_tree");
}