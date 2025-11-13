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

void single_ratio() {
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
    TCanvas *c1 = new TCanvas("c1", "", 7680, 4320);

    TCanvas *canvases[] = {c1};
    int numCanvases = 1;

    // Setting histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;

    TH1D *h1 = cHist("h1", "qinv[GeV]", "#Pairs/bin", 0., 1., ninterval, nlength, nscale, 920, 1., 1001, 632, 1, 1);
    TH1D *h2 = cHist("h2", "", "", 0., 1., ninterval, nlength, nscale, 920, 1., 1001, 632, 1, 1);

    TH1D *tonormhist[] = {h1, h2};
    int numTonorm = 2;

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        if (HFsumET > 100 && HFsumET < 375) { // This selects a centrality
            for (int k = 0; k < NOSpair; k++) {
                h1->Fill(qinvSigOS[k]);

            }
            for (int l = 0; l < NSSpair; l++) {
                h2->Fill(qinvSigSS[l]);
            }    
        }

        // Display progress
        float progress = (float)(i + 1) / nentries * 100;
        std::cout << "\rFilling Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }

    // Normalize histograms
    Double_t scale = 1;
    normalizer(tonormhist, numTonorm, scale);
   
    // Dividing SS/OS
    TH1D *sr = (TH1D *)h2->Clone("sr");
    sr->Divide(h1);

    // Adding labels
    sr->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV");
    sr->GetXaxis()->SetTitle("qinv[GeV]");
    sr->GetYaxis()->SetTitle("Single Ratio SS/OS");
    
    TText *text = new TText(0.8, 1.5, "50-70%");
    text->SetTextSize(0.03);
    
    // Overwriting bin 1 with bin 2 to remove big eror
    int nBin = difZeroBin(sr);
    
    // Setting y range to 0.95<y<1.6
    sr->GetYaxis()->SetRangeUser(0.95, 1.6);
    
    // Removing statistics box
    sr->SetStats(0);

    // Drawing the single ratio
    c1->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); sr->Draw("HIST"); text->Draw();
    // Saving image
    const char *path = "./imgs/teste/";
    const char *prefix = "teste-normalizar-sr";
    const char *file_type = "png";
    save_canvas_images(canvases, numCanvases, path, prefix, file_type);

    // Closing program
    TH1D *histograms[] = {h1, h2, sr};
    int numHistograms = 3;

    close_program(canvases, numCanvases, histograms, numHistograms, fr);
}
