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

void coulomb_qinv() {
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
    TCanvas *c2 = new TCanvas("c2", "", 7680, 4320);

    TCanvas *canvases[] = {c1, c2};
    int numCanvases = 2;

    double bin_scale = 50.;

    // Setting histograms
    TH2D *hOS = new TH2D("hOS", "CMS Open Data 2011 - PbPb 2.76 TeV", bin_scale, 0., 1., 100, 0, 2);
    TH2D *hSS = new TH2D("hSS", "CMS Open Data 2011 - PbPb 2.76 TeV", bin_scale, 0., 1., 100, 0, 2);

    TH2D *histograms[] = {hOS, hSS};
    int numHistograms = 2;

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        if (HFsumET > 100 && HFsumET < 375) { // This selects a centrality
            for (int k = 0; k < NOSpair; k++) {
                hOS->Fill(qinvSigOS[k], coulombWOS[k]);

            }
            for (int l = 0; l < NSSpair; l++) {
                hSS->Fill(qinvSigSS[l], coulombWSS[l]);
            }    
        }

        // Display progress
        float progress = (float)(i + 1) / nentries * 100;
        std::cout << "\rFilling Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }

    // Normalize histograms
    Double_t scale = 1;
    normalizer2d(histograms, numHistograms, scale);

    TProfile *profhOSX, *profhSSX;
    profhOSX = hOS->ProfileX();
    profhSSX = hSS->ProfileX();

    profhOSX->SetFillColor(kGray);
    profhSSX->SetFillColor(kGray);

    profhOSX->SetLineColor(kBlack);
    profhSSX->SetLineColor(kBlack);

    profhOSX->GetXaxis()->SetTitle("OS - qinv[GeV]"); profhOSX->GetYaxis()->SetTitle("OS - coulomb correction");
    profhSSX->GetXaxis()->SetTitle("SS - qinv[GeV]"); profhSSX->GetYaxis()->SetTitle("SS - coulomb correction");

    // Drawing
    c1->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); profhOSX->Draw("BAR");
    c2->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); profhSSX->Draw("BAR");
    
    // Saving image
    const char *path = "./imgs/final/";
    const char *prefix = "final-coulomb-qinv";
    const char *file_type = "png";
    save_canvas_images(canvases, numCanvases, path, prefix, file_type);

    // Closing program
    close_program(canvases, numCanvases, nullptr, numHistograms, fr, histograms);
}
