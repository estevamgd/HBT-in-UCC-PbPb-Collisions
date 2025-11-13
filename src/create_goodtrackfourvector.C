#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include "TLegend.h"

#include "../include/data_func.h"
#include "../include/my_func.h"

void create_goodtrackfourvector() {
    ROOT::EnableImplicitMT();
    auto poolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << poolSize << std::endl;

    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree("data/HiForestAOD_UPC.root", "demo/HBT", fr, t);

    // Variables from the tree
    Int_t maxSize = 1700, Ntrk;
    Float_t HFsumET, trkPt[maxSize], trkPtRes[maxSize], trkEta[maxSize], trkPhi[maxSize], pionMass = 0.13957;

    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    t->SetBranchAddress("trkPhi", trkPtRes);

    // Getting how many entries
    Long64_t nentries = t->GetEntries();

    // Creating a new ROOT file to store the output
    TFile *outputFile = new TFile("data/output_tracks.root", "RECREATE");

    // Creating a new TTree to store particles
    TTree *tree = new TTree("tracks", "Tree with tracks");

    // Defining the structure to store the TLorentzVector objects
    std::vector<TLorentzVector> TrackFourVector;
    //std::vector<int> TrackCharge;

    tree->Branch("TrackFourVector", &TrackFourVector);
    //tree->Branch("TrackCharge", &TrackCharge);

    // Setting Lorentz Vector
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        TrackFourVector.clear();
        //TrackCharge.clear();
 
        if (HFsumET > 100 && HFsumET < 375) continue;
        
        for (int j = 0; j < Ntrk; j++) {
            TLorentzVector vec;
            /*int charge = 1;
            
            if (trkPt[j] < 0) { 
                charge = -1;
            }*/
            vec.SetPtEtaPhiM(trkPt[j], trkEta[j], trkPhi[j], pionMass);

            TrackFourVector.push_back(vec);
            //TrackCharge.push_back(charge);
            std::cout << "nentries: " << nentries << " i: " << i << " Ntrk: " << Ntrk << " j: " << j << std::endl;
        }

        if (TrackFourVector.empty()) {
            std::cout << "Warning: TrackFourVector is empty for event " << i << std::endl;
        } else {
            tree->Fill();
        }

    }
    // Writing the tracks in the ROOT file
    tree->Write();
    outputFile->Close();
    
    std::cout << "tracks saved on the ROOT file!" << std::endl;
}