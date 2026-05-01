#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>

void tracks_vs_hf() {

    // --- Input ---
    const char* fileName = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";
    const char* treeName = "demo/TreeMBUCC";

    std::vector<double> bins = {3300, 3400, 3500, 3600, 3700, 3800};
    int nbins = bins.size() - 1;

    // --- Open file ---
    TFile *f = TFile::Open(fileName);
    TTree *t = (TTree*)f->Get(treeName);

    // --- Branch variables ---
    Float_t HFsumET;
    Int_t Ntrk;

    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);

    // --- Histograms ---
    TH1F *hSum = new TH1F("hSum", "", nbins, &bins[0]);
    TH1F *hCount = new TH1F("hCount", "", nbins, &bins[0]);

    // --- Event loop ---
    Long64_t nentries = t->GetEntries();
    std::cout << "Entries: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        hSum->Fill(HFsumET, Ntrk);  // sum of tracks
        hCount->Fill(HFsumET);      // number of events
    }

    // --- Mean tracks ---
    TH1F *hMean = (TH1F*)hSum->Clone("hMean");
    hMean->SetTitle("<N_{trk}> vs HFsumET;HFsumET;<N_{trk}>");
    hMean->Divide(hCount);

    // --- Plot ---
    TCanvas *c = new TCanvas("c", "", 800, 600);
    hMean->SetMarkerStyle(20);
    hMean->Draw("E");

    c->SaveAs("tracks_vs_hf.png");

    // --- Print values ---
    std::cout << "\nResults:\n";
    for (int i = 1; i <= nbins; i++) {
        std::cout << "[" << bins[i-1] << ", " << bins[i] << "] : "
                  << hMean->GetBinContent(i) << std::endl;
    }
}