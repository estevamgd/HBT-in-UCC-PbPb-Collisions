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

#include "../include/data_func.h"
#include "../include/my_func.h"

struct TrackInfo {
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> p4;
    Int_t charge;    
    Float_t weight;  
};

void sigmix_test(){
    ROOT::EnableImplicitMT();
    auto threadPoolSize = ROOT::GetThreadPoolSize();
    std::cout << "Pool size = " << threadPoolSize << std::endl;

    // Load the ROOT file and tree
    TFile *fr = nullptr;
    TTree *t = nullptr;
    getFileTree("data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root", "demo/TreeMBUCC", fr, t);

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

    // Getting how many entries
    Long64_t nentries = t->GetEntries(), trackvec_max_size = 0;

    // --- POOL SETTINGS ---
    std::vector<std::vector<TrackInfo>> poolEvents;
    std::vector<std::vector<TrackInfo>> phiEvents;

    int poolSize = 10;
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

    int hiBinLimInf = 103, hiBinLimSup = 104;

    for (Long64_t globalEventIndex = 0; globalEventIndex < nentries; globalEventIndex++){
        t->GetEntry(globalEventIndex);
        
        
        std::vector<TrackInfo> currentEventTracks;

        for (int globalTrackIndex = 0; globalTrackIndex < Ntrk; globalTrackIndex++) {
            // Build Track'ij
            TrackInfo currentTrack;
            currentTrack.p4 = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[globalTrackIndex], trkEta[globalTrackIndex], trkPhi[globalTrackIndex], pionMass);
            currentTrack.charge = trkCharge[globalTrackIndex];
            currentTrack.weight = trkWeight[globalTrackIndex];
            
            // === MIXING ===
            if (poolEvents.size() == poolSize - 1) {
                // Double Loop
                for (int poolEventIndex = 0; poolEventIndex < poolEvents.size(); poolEventIndex++) {
                    for (int poolTrackIndex = 0; poolTrackIndex < poolEvents[poolEventIndex].size(); poolTrackIndex++) {
                        
                    }
                }

                // Single Loop
                for (int poolEventIndexIndex = 0; poolEventIndexIndex < poolMaxTracks*poolEvents.size(); poolEventIndexIndex++) {
                    int poolEventIndex = poolEventIndexIndex / poolMaxTracks;
                    int poolTrackIndex = poolEventIndexIndex % poolMaxTracks;

                    if (poolTrackIndex >= poolEvents[poolEventIndex].size()) continue;
                    if (hiBin < hiBinLimInf || hiBin > hiBinLimSup) continue;
                    std::cout << "mix currentTrack " << globalEventIndex << " " << globalTrackIndex << ": " << currentTrack.p4 << std::endl;
                    std::cout << "mix poolEvents[" << poolEventIndex << "][" << poolTrackIndex << "]: " << poolEvents[poolEventIndex][poolTrackIndex].p4 << std::endl;
                    std::cout << std::endl;
                }
            }
            
            // Build Event i
            currentEventTracks.push_back(currentTrack);
        }
        // === SIGNAL ===
        if (currentEventTracks.size() > 1) {
            // Double Loop
            for (int signalTrackIndex = 0; signalTrackIndex < currentEventTracks.size(); signalTrackIndex++){
                for (int singalTrack2Index = signalTrackIndex + 1; singalTrack2Index < currentEventTracks.size(); singalTrack2Index++){
                    
                }
            }
            // Single Loop
            for (int signalTrackIndexIndex = 0; signalTrackIndexIndex < currentEventTracks.size()*(currentEventTracks.size()-1); signalTrackIndexIndex++) {
                int singalTrackIndex = signalTrackIndexIndex / (currentEventTracks.size() - 1);
                int singalTrack2Index = (signalTrackIndexIndex % (currentEventTracks.size() - 1)) + 1;

                if (singalTrack2Index <= singalTrackIndex) continue;
                std::cout << "sig currentEventTracks[" << singalTrackIndex << "]: " << currentEventTracks[singalTrackIndex].p4 << std::endl;
                std::cout << "sig currentEventTracks[" << singalTrack2Index << "]: " << currentEventTracks[singalTrack2Index].p4 << std::endl;
                std::cout << std::endl;
            }
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

        x++;
    }
    // === PHI MIXING ===
    if (phiSwitch) {
        for (int event = 0; event < poolEvents.size(); event++) {
            phiEvents.push_back(poolEvents[event]);
        }

        if (poolMaxTracks > phiMaxTracks) phiMaxTracks = poolMaxTracks;

        int lastIndex = phiEvents.size() - 1;
        int lastEventSize = phiEvents[lastIndex].size();

        for (int lastEventTrackID = 0; lastEventTrackID < lastEventSize; lastEventTrackID++) {
            // Double Loop
            for (int k = 0; k < phiEvents.size() - 1; k++) {
                for (int z = 0; z < phiEvents[k].size(); z++) {
                    
                }
            }
            // Single Loop
            for (int w = 0; w < phiMaxTracks*(phiEvents.size()-1); w++) {
                int k = w / phiMaxTracks;
                int z = w % phiMaxTracks;

                if (z > phiEvents[k].size() - 1) continue;
                if (hiBin < hiBinLimInf || hiBin > hiBinLimSup) continue;
                std::cout << "mix2 phiEvents[" << lastIndex << "][" << lastEventTrackID << "]: " << phiEvents[lastIndex][lastEventTrackID].p4 << std::endl;
                std::cout << "mix2 phiEvents[" << k << "][" << z << "]: "<< phiEvents[k][z].p4 << std::endl;
                std::cout << std::endl;
            }

        }
    }

}