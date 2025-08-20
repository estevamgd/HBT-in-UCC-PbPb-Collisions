#include <iostream>
#include <vector>
#include <cmath>

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
#include <TSystem.h>

void testHFsumET(const char *fileInput, const char *treeInput) { 
    ROOT::EnableImplicitMT();
    auto threadPoolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << threadPoolSize << std::endl;

    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree(fileInput, treeInput, fr, t);

    // Variables from the tree
    Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
    Float_t HFsumET, pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize], 
            pionMass = 0.13957039; // Pion mass [GeV] from PDG

    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("pvZ", &pvZ); // Z in cm; Use fabs(pvZ) for absolute value of the vertex position in z axis
    t->SetBranchAddress("trkCharge", trkCharge);
    t->SetBranchAddress("trkWeight", trkWeight);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    
    for (int i = 0; i < 5; i++){
        t->GetEntry(i);
        std::cout << "testHFsumEt> hiBin: " << hiBin << " HFsumET: " << HFsumET << std::endl;
    }
}

void testHiBin(const char *fileInput, const char *treeInput) { 
    ROOT::EnableImplicitMT();
    auto threadPoolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << threadPoolSize << std::endl;

    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree(fileInput, treeInput, fr, t);

    // Variables from the tree
    Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
    Float_t HFsumET, pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize], 
            pionMass = 0.13957039; // Pion mass [GeV] from PDG

    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("pvZ", &pvZ); // Z in cm; Use fabs(pvZ) for absolute value of the vertex position in z axis
    t->SetBranchAddress("trkCharge", trkCharge);
    t->SetBranchAddress("trkWeight", trkWeight);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    
    for (int i = 0; i < 5; i++){
        t->GetEntry(i);
        std::cout << "testHibin> hiBin: " << hiBin << " HFsumET: " << HFsumET << std::endl;
    }
}

void exectest(){
    const char* inputFile = "data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root";  
    const char* treeName  = "demo/TreeMBUCC";

    testHFsumET(inputFile, treeName);
    testHiBin(inputFile, treeName);
}