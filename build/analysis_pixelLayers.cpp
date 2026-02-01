#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include "TCanvas.h"
#include "TH1D.h"
#include <iostream>
#include <iomanip>
#include "../include/my_func.h"
#include "../include/analyze_tools.h"

void analysis_pixelLayers() {
    // 1. Setup File and Tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    
    const char* fileName = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";
    
    getFileTree(fileName, "demo/TreeMBUCC", fr, t);

    // --- SAFETY CHECK ---
    if (!t || !fr) {
        std::cerr << "ERROR: Could not open file or tree at: " << fileName << std::endl;
        std::cerr << "Please check that the file exists and the path is correct." << std::endl;
        return; 
    }

    // 2. Define Variables
    const Int_t maxSize = 50000;
    Int_t Ntrk;
    Float_t trkNpixLayers[maxSize]; 

    // 3. Set Branch Addresses
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("trkNpixLayers", trkNpixLayers);

    // 4. Setup Histogram and Canvas
    TH1D* htrkNpixLayers = cHist("htrkNpixLayers", "Number of Pixel Layers", "#Tracks/bin", 0, 6, 1., 0.02, 1.);
    
    TCanvas *c1 = new TCanvas("c1", "Pixel Layers", 800, 600);

    // 5. Loop and Fill
    Long64_t nentries = t->GetEntries();
    
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        for (int j = 0; j < Ntrk; j++) {
            htrkNpixLayers->Fill(trkNpixLayers[j]);
        }
        
        if (i % 1000 == 0) {
            float progress = (float)(i + 1) / nentries * 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
        }
    }
    std::cout << std::endl;

    // --- PRINT TRACKS WITH 0 LAYERS ---
    // We look for the bin containing the value 0.
    int bin0 = htrkNpixLayers->FindBin(0.0);
    double tracksAtZero = htrkNpixLayers->GetBinContent(bin0);
    
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Tracks with Number of Pixel Layers = 0: " << (long long)tracksAtZero << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    // 6. Draw and Clean up
    c1->cd();
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy();
    
    htrkNpixLayers->SetFillColor(kBlue-9);
    htrkNpixLayers->SetLineColor(kBlue+2);
    htrkNpixLayers->Draw();

    // Arrays for your custom close/save functions
    TCanvas *canvases[] = {c1};
    TH1D* histograms[] = {htrkNpixLayers};
    TLegend *legends[] = {}; 

    // Save image
    save_canvas_images(canvases, 1, "./imgs/test/", "pixel-layers", "png");

    // Close program
    close_program(canvases, 1, histograms, 1, legends, 0, fr);
}