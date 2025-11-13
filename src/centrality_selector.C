#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TH1D.h"
#include <TStyle.h>
#include "TLegend.h"
#include <TText.h>
#include <TBenchmark.h>

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"


void centrality_selector(){
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

    // Getting how many entries
    Long64_t nentries = t->GetEntries();

    // Setting canvases
    TCanvas *c1 = new TCanvas("c1", "Invariant Difference in 4-Momentum", 7680, 4320);
    TCanvas *c2 = new TCanvas("c2", "Invariant Difference in 4-Momentum", 7680, 4320);
    TCanvas *c3 = new TCanvas("c3", "Invariant Difference in 4-Momentum", 7680, 4320);

    TCanvas *canvases[] = {c1, c2, c3};
    int numCanvases = 3;
    
    // Setting histograms
    TH1D *h1 = cHist("h1", "qinv[GeV]", "#Pairs/bin", nentries, -0.1, 1.3, 0, 0, 0, 632, 1, 1);
    TH1D *h2 = cHist("h2", "", "", nentries, -0.1, 1.3, 0, 0, 0, 600, 1, 1); // <- tentar 0.3, 0.5 ou tentar outro tipo de fill style/ sem fill style
    TH1D *h3 = cHist("h3", "", "", nentries, -0.1, 1.3, 0, 0, 0, 632, 1, 2);
    TH1D *h4 = cHist("h4", "", "", nentries, -0.1, 1.3, 0, 0, 0, 600, 1, 2);
    TH1D *h5 = cHist("h5", "HFsumET[GeV]", "#Events/bin", nentries, 80, 390, 632);

    TH1D *histograms[] = {h1, h2, h3, h4, h5};
    int numHistograms = 5;

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        if (HFsumET > 100 && HFsumET < 375) {
            h5->Fill(HFsumET);

            for (int k = 0; k < NOSpair; k++) {
                h1->Fill(qinvSigOS[k]);

            }
            for (int l = 0; l < NSSpair; l++) {
                h2->Fill(qinvSigSS[l]);
            }    
        }

        for (int k = 0; k < NOSpair; k++) {
            h3->Fill(qinvSigOS[k]);

        }
        for (int l = 0; l < NSSpair; l++) {
            h4->Fill(qinvSigSS[l]);
        } 

        // Display progress
        float progress = (float)(i + 1) / nentries * 100;
        std::cout << "\rFilling Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }
    
    // Create the legend
    TLegend *legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend->AddEntry(h1, "OS qinvSig (100 < HFsumET < 375 GeV)", "l");
    legend->AddEntry(h2, "SS qinvSig (100 < HFsumET < 375 GeV)", "l");
    legend->AddEntry(h3, "OS qinvSig (All HFsumET)", "l");
    legend->AddEntry(h4, "SS qinvSig (All HFsumET)", "l");

    TLegend *legend2 = new TLegend(0.2, 0.8, 0.4, 0.9);
    legend2->AddEntry(h1, "OS qinvSig (100 < HFsumET < 375 GeV)", "l");
    legend2->AddEntry(h2, "SS qinvSig (100 < HFsumET < 375 GeV)", "l");

    c1->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); h1->Draw(); h2->Draw("same"); h3->Draw("same"); h4->Draw("same");
    legend->Draw();

    c2->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); h1->Draw(); h2->Draw("same"); 
    legend2->Draw();

    c3->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); h5->Draw();

    const char *path = "./imgs/final/";
    const char *prefix = "final-central-selec";
    const char *file_type = "png";
    save_canvas_images(canvases, numCanvases, path, prefix, file_type);
    close_program(canvases, numCanvases, histograms, numHistograms, fr);
}