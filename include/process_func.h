#ifndef PROCESS_FUNC_H
#define PROCESS_FUNC_H

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <TStyle.h>
#include <fstream>
#include <ctime>
#include <sstream>
#include <filesystem>
#include "my_func.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/Boost.h"
#include "TMath.h"
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>
#include <chrono>
#include <string>
#include "data_func.h"
#include <TSystem.h>
#include "TVector2.h"


#define MIN_PT_CUT 0.5 // GeV/c
#define MAX_PT_CUT 2.0 // GeV/c
#define ABS_ETA_CUT 0.95
#define N_HITS 1.0

// ===== Process Signal Qinv =====
// ----- Signal -----
// Serialized version
void processSignalQinv(
    size_t i_n_tracks,
    size_t j_n_tracks,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSS_signal_2l,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOS_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)

{
    if (i_n_tracks <= 1) return; // check if event has 2 or more tracks

    for (size_t p1 = 0; p1 < i_n_tracks; p1++) {
        for (size_t p2 = p1+1; p2 < j_n_tracks; p2++) {
            // Checks
            if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue; // Check if weight is infinity
            // Calculate q_inv
            double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);
            if (trkCharge[p1] * trkCharge[p2] > 0) { // Fills same charge particle pair histogram
                h_qinvSS_signal_2l->Fill(qinv);
                h_qinvSSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
            } else { // Fills opposite charge particle pair histogram
                h_qinvOS_signal_2l->Fill(qinv);
                h_qinvOSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
            }
        }
    }
}

// Miltithreaded version
void processSignalQinv(
    size_t n_tracks,
    int thread_count,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSS_signal_2l,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOS_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)

{
    if (n_tracks <= 1) return; // check if event has 2 or more tracks
    if (thread_count == 0) thread_count = 1;

    // A simple struct to hold the set of histograms for one thread
    struct LocalHistograms {
        TH1D* ss;
        TH1D* ss_cor;
        TH1D* os;
        TH1D* os_cor;

        // Constructor to clone the main histograms
        LocalHistograms(const TH1D* base_ss, const TH1D* base_ss_cor, const TH1D* base_os, const TH1D* base_os_cor) {
            // Give each clone a unique name to avoid ROOT conflicts
            ss      = (TH1D*)base_ss->Clone(TString::Format("%s_clone_%p", base_ss->GetName(), this));
            ss_cor  = (TH1D*)base_ss_cor->Clone(TString::Format("%s_clone_%p", base_ss_cor->GetName(), this));
            os      = (TH1D*)base_os->Clone(TString::Format("%s_clone_%p", base_os->GetName(), this));
            os_cor  = (TH1D*)base_os_cor->Clone(TString::Format("%s_clone_%p", base_os_cor->GetName(), this));

            // Reset clones to be empty
            ss->Reset();
            ss_cor->Reset();
            os->Reset();
            os_cor->Reset();
        }
        // Destructor to clean up cloned histograms
        ~LocalHistograms() {
            delete ss;
            delete ss_cor;
            delete os;
            delete os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    // Create a set of local histograms for each thread
    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l));
    }
    // The task for each thread
    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;

                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    hists->ss->Fill(qinv);
                    hists->ss_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else {
                    hists->os->Fill(qinv);
                    hists->os_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
    };
    // Distribute the work among threads
    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;

        // Each thread gets a pointer to its own set of histograms
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());

    }
    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Merge the results from local histograms into the main ones
    for (int t = 0; t < thread_count; ++t) {
        h_qinvSS_signal_2l->Add(local_hists[t]->ss);
        h_qinvSSCor_signal_2l->Add(local_hists[t]->ss_cor);
        h_qinvOS_signal_2l->Add(local_hists[t]->os);
        h_qinvOSCor_signal_2l->Add(local_hists[t]->os_cor);

    }
}

// Miltithreaded version with pT cut
void processSignalQinv(
    size_t n_tracks,
    int thread_count,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const Float_t* trkPt,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSS_signal_2l,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOS_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)

{
    if (n_tracks <= 1) return; // check if event has 2 or more tracks
    if (thread_count == 0) thread_count = 1;

    // A simple struct to hold the set of histograms for one thread
    struct LocalHistograms {
        TH1D* ss;
        TH1D* ss_cor;
        TH1D* os;
        TH1D* os_cor;

        // Constructor to clone the main histograms
        LocalHistograms(const TH1D* base_ss, const TH1D* base_ss_cor, const TH1D* base_os, const TH1D* base_os_cor) {
            // Give each clone a unique name to avoid ROOT conflicts
            ss      = (TH1D*)base_ss->Clone(TString::Format("%s_clone_%p", base_ss->GetName(), this));
            ss_cor  = (TH1D*)base_ss_cor->Clone(TString::Format("%s_clone_%p", base_ss_cor->GetName(), this));
            os      = (TH1D*)base_os->Clone(TString::Format("%s_clone_%p", base_os->GetName(), this));
            os_cor  = (TH1D*)base_os_cor->Clone(TString::Format("%s_clone_%p", base_os_cor->GetName(), this));

            // Reset clones to be empty
            ss->Reset();
            ss_cor->Reset();
            os->Reset();
            os_cor->Reset();
        }
        // Destructor to clean up cloned histograms
        ~LocalHistograms() {
            delete ss;
            delete ss_cor;
            delete os;
            delete os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    // Create a set of local histograms for each thread
    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l));
    }
    // The task for each thread
    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;
                if (trkPt[p1] <= MIN_PT_CUT || trkPt[p2] <= MIN_PT_CUT) continue; // pT cut

                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    hists->ss->Fill(qinv);
                    hists->ss_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else {
                    hists->os->Fill(qinv);
                    hists->os_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
    };
    // Distribute the work among threads
    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;

        // Each thread gets a pointer to its own set of histograms
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());

    }
    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Merge the results from local histograms into the main ones
    for (int t = 0; t < thread_count; ++t) {
        h_qinvSS_signal_2l->Add(local_hists[t]->ss);
        h_qinvSSCor_signal_2l->Add(local_hists[t]->ss_cor);
        h_qinvOS_signal_2l->Add(local_hists[t]->os);
        h_qinvOSCor_signal_2l->Add(local_hists[t]->os_cor);

    }
}

// Multithreaded version with pT cut and eta cut, only Corrected hists
void processSignalQinv(
    size_t n_tracks,
    int thread_count,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)
{
    if (n_tracks <= 1) return; // check if event has 2 or more tracks
    if (thread_count == 0) thread_count = 1;

    // A simple struct to hold the set of histograms for one thread
    struct LocalHistograms {
        TH1D* ss_cor;
        TH1D* os_cor;

        // Constructor to clone the main histograms
        LocalHistograms(const TH1D* base_ss_cor, const TH1D* base_os_cor) {
            // Give each clone a unique name to avoid ROOT conflicts
            ss_cor  = (TH1D*)base_ss_cor->Clone(TString::Format("%s_clone_%p", base_ss_cor->GetName(), this));
            os_cor  = (TH1D*)base_os_cor->Clone(TString::Format("%s_clone_%p", base_os_cor->GetName(), this));

            // Reset clones to be empty
            ss_cor->Reset();
            os_cor->Reset();
        }
        // Destructor to clean up cloned histograms
        ~LocalHistograms() {
            delete ss_cor;
            delete os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    // Create a set of local histograms for each thread
    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSSCor_signal_2l, h_qinvOSCor_signal_2l));
    }
    // The task for each thread
    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;
                if (currentEventTracks4V[p1].Pt() <= MIN_PT_CUT 
                    || currentEventTracks4V[p2].Pt() <= MIN_PT_CUT) continue; // pT cut
                if (std::abs(currentEventTracks4V[p1].Eta()) >= ABS_ETA_CUT 
                    || std::abs(currentEventTracks4V[p2].Eta()) >= ABS_ETA_CUT) continue; // eta cut

                double qinv = GetQ(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    hists->ss_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else {
                    hists->os_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
    };
    // Distribute the work among threads
    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;

        // Each thread gets a pointer to its own set of histograms
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());

    }
    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Merge the results from local histograms into the main ones
    for (int t = 0; t < thread_count; ++t) {
        h_qinvSSCor_signal_2l->Add(local_hists[t]->ss_cor);
        h_qinvOSCor_signal_2l->Add(local_hists[t]->os_cor);

    }
}

// ----- Mixing -----
// Serialized version
void processMixQinv(
    size_t poolSize,
    Float_t vertexDistance,
    Int_t pairChargeMult,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    TH1D* h_qinv_mix_2l,
    TH1D* h_qinvCor_mix_2l)
{
    if (currentEventTracks4V.empty()) return;
    if (vectorEventTracks4V.empty()) return;
    if (vectorEventTracks4V.size() >= poolSize) return;
    if (pvZ >= vertexDistance) return;

    for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
        for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
            for (size_t p2 = 0; p2 < currentEventTracks4V.size(); ++p2) {
                if (vectorEventTracksCharge[nEv][p1] * currentEventTracksCharge[p2] != pairChargeMult) continue;
                if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;

                double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                h_qinv_mix_2l->Fill(qinv);
                h_qinvCor_mix_2l->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
            }
        }
    }
    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);

    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
    }
}

// Multithreaded version
void processMixQinv(
    Int_t n_tracks,
    int thread_count,
    size_t poolSize,
    Float_t vertexDistance,
    Int_t pairChargeMult,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    TH1D* h_qinv_mix_2l,
    TH1D* h_qinvCor_mix_2l)
{   
    if (pvZ >= vertexDistance) return;
    if (!currentEventTracks4V.empty()){
        if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()){
            if (thread_count == 0) thread_count = 1;
            struct LocalHistograms {
                TH1D* mix;
                TH1D* mix_cor;
        
                LocalHistograms(const TH1D* base_mix, const TH1D* base_mix_cor) {
                    mix      = (TH1D*)base_mix->Clone(TString::Format("%s_clone_%p", base_mix->GetName(), this));
                    mix_cor  = (TH1D*)base_mix_cor->Clone(TString::Format("%s_clone_%p", base_mix_cor->GetName(), this));
                    mix->Reset();
                    mix_cor->Reset();
                }

                ~LocalHistograms() {
                    delete mix;
                    delete mix_cor;
                }
            };
        
            std::vector<std::thread> threads;
            std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        
            for (int i = 0; i < thread_count; ++i) {
                local_hists.push_back(std::make_unique<LocalHistograms>(h_qinv_mix_2l, h_qinvCor_mix_2l));
            }
        
            auto thread_task = [&](size_t start_p1, size_t end_p1, size_t nEv, LocalHistograms* hists) {
                for (size_t p2 = start_p1; p2 < end_p1; ++p2) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
                        if (vectorEventTracksCharge[nEv][p1] * currentEventTracksCharge[p2] != pairChargeMult) continue;
                        if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;
        
                        double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                        hists->mix->Fill(qinv);
                        hists->mix_cor->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
                    }
                }
            };
        
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
                size_t chunk_size = n_tracks / thread_count;
                for (int t = 0; t < thread_count; t++) {
                    size_t start_i = t * chunk_size;
                    size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
            
                    threads.emplace_back(thread_task, start_i, end_i, nEv, local_hists[t].get());
                }
                for (auto& thread : threads) {
                    thread.join();
                }

                threads.clear();

                for (int t = 0; t < thread_count; ++t) {
                    h_qinv_mix_2l->Add(local_hists[t]->mix);
                    h_qinvCor_mix_2l->Add(local_hists[t]->mix_cor);
                    local_hists[t]->mix->Reset();
                    local_hists[t]->mix_cor->Reset();
                }
            }
        }
    }

    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);
    
    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
    }
}

// Multithreaded version with pT cut
void processMixQinv(
    Int_t n_tracks,
    int thread_count,
    size_t poolSize,
    Float_t vertexDistance,
    Int_t pairChargeMult,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    std::vector<std::vector<Float_t>>& vectorEventTracksPt,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    const std::vector<Float_t>& currentEventTracksPt,
    TH1D* h_qinv_mix_2l,
    TH1D* h_qinvCor_mix_2l)
{   
    if (pvZ >= vertexDistance) return;
    if (!currentEventTracks4V.empty()){
        if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()){
            if (thread_count == 0) thread_count = 1;
            struct LocalHistograms {
                TH1D* mix;
                TH1D* mix_cor;
        
                LocalHistograms(const TH1D* base_mix, const TH1D* base_mix_cor) {
                    mix      = (TH1D*)base_mix->Clone(TString::Format("%s_clone_%p", base_mix->GetName(), this));
                    mix_cor  = (TH1D*)base_mix_cor->Clone(TString::Format("%s_clone_%p", base_mix_cor->GetName(), this));
                    mix->Reset();
                    mix_cor->Reset();
                }

                ~LocalHistograms() {
                    delete mix;
                    delete mix_cor;
                }
            };
        
            std::vector<std::thread> threads;
            std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        
            for (int i = 0; i < thread_count; ++i) {
                local_hists.push_back(std::make_unique<LocalHistograms>(h_qinv_mix_2l, h_qinvCor_mix_2l));
            }
        
            auto thread_task = [&](size_t start_p1, size_t end_p1, size_t nEv, LocalHistograms* hists) {
                for (size_t p2 = start_p1; p2 < end_p1; ++p2) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
                        if (vectorEventTracksCharge[nEv][p1] * currentEventTracksCharge[p2] != pairChargeMult) continue;
                        if (vectorEventTracksPt[nEv][p1] <= MIN_PT_CUT || currentEventTracksPt[p2] <= MIN_PT_CUT) continue; // pT cut

                        if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;
        
                        double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                        hists->mix->Fill(qinv);
                        hists->mix_cor->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
                    }
                }
            };
        
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
                size_t chunk_size = n_tracks / thread_count;
                for (int t = 0; t < thread_count; t++) {
                    size_t start_i = t * chunk_size;
                    size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
            
                    threads.emplace_back(thread_task, start_i, end_i, nEv, local_hists[t].get());
                }
                for (auto& thread : threads) {
                    thread.join();
                }

                threads.clear();

                for (int t = 0; t < thread_count; ++t) {
                    h_qinv_mix_2l->Add(local_hists[t]->mix);
                    h_qinvCor_mix_2l->Add(local_hists[t]->mix_cor);
                    local_hists[t]->mix->Reset();
                    local_hists[t]->mix_cor->Reset();
                }
            }
        }
    }

    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);
    vectorEventTracksPt.push_back(currentEventTracksPt);
    
    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
        vectorEventTracksPt.clear();
    }
}

// Multithreaded version with pT cut and eta cut, only Corrected hists
void processMixQinv(
    Int_t n_tracks,
    int thread_count,
    size_t poolSize,
    Float_t vertexDistance,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    TH1D* h_qinvSSCor_mix_2l,
    TH1D* h_qinvOSCor_mix_2l)
{   
    if (pvZ >= vertexDistance) return;
    if (!currentEventTracks4V.empty()){
        if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()){
            if (thread_count == 0) thread_count = 1;
            struct LocalHistograms {
                TH1D* mix_SS_cor;
                TH1D* mix_OS_cor;
        
                LocalHistograms(const TH1D* base_mix_SS_cor, const TH1D* base_mix_OS_cor) {
                    mix_SS_cor = (TH1D*)base_mix_SS_cor->Clone(TString::Format("%s_clone_%p", base_mix_SS_cor->GetName(), this));
                    mix_OS_cor = (TH1D*)base_mix_OS_cor->Clone(TString::Format("%s_clone_%p", base_mix_OS_cor->GetName(), this));
                    mix_SS_cor->Reset();
                    mix_OS_cor->Reset();
                }

                ~LocalHistograms() {
                    delete mix_SS_cor;
                    delete mix_OS_cor;
                }
            };
        
            std::vector<std::thread> threads;
            std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        
            for (int i = 0; i < thread_count; ++i) {
                local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSSCor_mix_2l, h_qinvOSCor_mix_2l));
            }
        
            auto thread_task = [&](size_t start_p1, size_t end_p1, size_t nEv, LocalHistograms* hists) {
                for (size_t p2 = start_p1; p2 < end_p1; ++p2) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
                        if (vectorEventTracks4V[nEv][p1].Pt() <= MIN_PT_CUT 
                            || currentEventTracks4V[p2].Pt() <= MIN_PT_CUT) continue; // pT cut
                        if (std::abs(vectorEventTracks4V[nEv][p1].Eta()) >= ABS_ETA_CUT 
                            || std::abs(currentEventTracks4V[p2].Eta()) >= ABS_ETA_CUT) continue; // pT cut
                        if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;
        
                        double qinv = GetQ(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);

                        if (vectorEventTracksCharge[nEv][p1]*currentEventTracksCharge[p2] > 0){
                            hists->mix_SS_cor->Fill(qinv, vectorEventTracksWeight[nEv][p1]*currentEventTracksWeight[p2]);
                        } else {
                            hists->mix_OS_cor->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
                        }
                    }
                }
            };
        
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
                size_t chunk_size = n_tracks / thread_count;
                for (int t = 0; t < thread_count; t++) {
                    size_t start_i = t * chunk_size;
                    size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
            
                    threads.emplace_back(thread_task, start_i, end_i, nEv, local_hists[t].get());
                }
                for (auto& thread : threads) {
                    thread.join();
                }

                threads.clear();

                for (int t = 0; t < thread_count; ++t) {
                    h_qinvSSCor_mix_2l->Add(local_hists[t]->mix_SS_cor);
                    h_qinvOSCor_mix_2l->Add(local_hists[t]->mix_OS_cor);
                    local_hists[t]->mix_SS_cor->Reset();
                    local_hists[t]->mix_OS_cor->Reset();
                }
            }
        }
    }

    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);
    
    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
    }
}

// ===== Process Signal QLCMS =====
// ----- Signal -----
// Serialized version
void processSignalQLCMS(
    size_t i_n_tracks,
    size_t j_n_tracks,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSS_signal_2l,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOS_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)

{
    if (i_n_tracks <= 1) return; // check if event has 2 or more tracks

    for (size_t p1 = 0; p1 < i_n_tracks; p1++) {
        for (size_t p2 = p1+1; p2 < j_n_tracks; p2++) {
            // Checks
            if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue; // Check if weight is infinity
            // Calculate q_inv
            double qinv = GetQLCMS(currentEventTracks4V[p1], currentEventTracks4V[p2]);
            if (trkCharge[p1] * trkCharge[p2] > 0) { // Fills same charge particle pair histogram
                h_qinvSS_signal_2l->Fill(qinv);
                h_qinvSSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
            } else { // Fills opposite charge particle pair histogram
                h_qinvOS_signal_2l->Fill(qinv);
                h_qinvOSCor_signal_2l->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
            }
        }
    }
}

// Multithreaded version
void processSignalQLCMS(
    size_t n_tracks,
    int thread_count,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSS_signal_2l,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOS_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)

{
    if (n_tracks <= 1) return; // check if event has 2 or more tracks
    if (thread_count == 0) thread_count = 1;

    // A simple struct to hold the set of histograms for one thread
    struct LocalHistograms {
        TH1D* ss;
        TH1D* ss_cor;
        TH1D* os;
        TH1D* os_cor;

        // Constructor to clone the main histograms
        LocalHistograms(const TH1D* base_ss, const TH1D* base_ss_cor, const TH1D* base_os, const TH1D* base_os_cor) {
            // Give each clone a unique name to avoid ROOT conflicts
            ss      = (TH1D*)base_ss->Clone(TString::Format("%s_clone_%p", base_ss->GetName(), this));
            ss_cor  = (TH1D*)base_ss_cor->Clone(TString::Format("%s_clone_%p", base_ss_cor->GetName(), this));
            os      = (TH1D*)base_os->Clone(TString::Format("%s_clone_%p", base_os->GetName(), this));
            os_cor  = (TH1D*)base_os_cor->Clone(TString::Format("%s_clone_%p", base_os_cor->GetName(), this));

            // Reset clones to be empty
            ss->Reset();
            ss_cor->Reset();
            os->Reset();
            os_cor->Reset();
        }
        // Destructor to clean up cloned histograms
        ~LocalHistograms() {
            delete ss;
            delete ss_cor;
            delete os;
            delete os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    // Create a set of local histograms for each thread
    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l));
    }
    // The task for each thread
    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;

                double qinv = GetQLCMS(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    hists->ss->Fill(qinv);
                    hists->ss_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else {
                    hists->os->Fill(qinv);
                    hists->os_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
    };
    // Distribute the work among threads
    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;

        // Each thread gets a pointer to its own set of histograms
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());

    }
    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Merge the results from local histograms into the main ones
    for (int t = 0; t < thread_count; ++t) {
        h_qinvSS_signal_2l->Add(local_hists[t]->ss);
        h_qinvSSCor_signal_2l->Add(local_hists[t]->ss_cor);
        h_qinvOS_signal_2l->Add(local_hists[t]->os);
        h_qinvOSCor_signal_2l->Add(local_hists[t]->os_cor);

    }
}

// Multithreaded version with pT cut and eta cut, only Corrected hists
void processSignalQLCMS(
    size_t n_tracks,
    int thread_count,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH1D* h_qinvSSCor_signal_2l,
    TH1D* h_qinvOSCor_signal_2l)

{
    if (n_tracks <= 1) return; // check if event has 2 or more tracks
    if (thread_count == 0) thread_count = 1;

    // A simple struct to hold the set of histograms for one thread
    struct LocalHistograms {
        TH1D* ss_cor;
        TH1D* os_cor;

        // Constructor to clone the main histograms
        LocalHistograms(const TH1D* base_ss_cor, const TH1D* base_os_cor) {
            // Give each clone a unique name to avoid ROOT conflicts
            ss_cor  = (TH1D*)base_ss_cor->Clone(TString::Format("%s_clone_%p", base_ss_cor->GetName(), this));
            os_cor  = (TH1D*)base_os_cor->Clone(TString::Format("%s_clone_%p", base_os_cor->GetName(), this));

            // Reset clones to be empty
            ss_cor->Reset();
            os_cor->Reset();
        }
        // Destructor to clean up cloned histograms
        ~LocalHistograms() {
            delete ss_cor;
            delete os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    // Create a set of local histograms for each thread
    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSSCor_signal_2l, h_qinvOSCor_signal_2l));
    }
    // The task for each thread
    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;
                if (currentEventTracks4V[p1].Pt() <= MIN_PT_CUT 
                    || currentEventTracks4V[p2].Pt() <= MIN_PT_CUT) continue; // pT cut
                if (std::abs(currentEventTracks4V[p1].Eta()) >= ABS_ETA_CUT 
                    || std::abs(currentEventTracks4V[p2].Eta()) >= ABS_ETA_CUT) continue; // eta cut

                double qinv = GetQLCMS(currentEventTracks4V[p1], currentEventTracks4V[p2]);

                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    hists->ss_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                } else {
                    hists->os_cor->Fill(qinv, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
    };
    // Distribute the work among threads
    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;

        // Each thread gets a pointer to its own set of histograms
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());

    }
    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Merge the results from local histograms into the main ones
    for (int t = 0; t < thread_count; ++t) {
        h_qinvSSCor_signal_2l->Add(local_hists[t]->ss_cor);
        h_qinvOSCor_signal_2l->Add(local_hists[t]->os_cor);

    }
}

// ----- Mixing -----
// Serialized version
void processMixQLCMS(
    size_t poolSize,
    Float_t vertexDistance,
    Int_t pairChargeMult,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    TH1D* h_qinv_mix_2l,
    TH1D* h_qinvCor_mix_2l)
{
    if (currentEventTracks4V.empty()) return;
    if (vectorEventTracks4V.empty()) return;
    if (vectorEventTracks4V.size() >= poolSize) return;
    if (pvZ >= vertexDistance) return;

    for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
        for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
            for (size_t p2 = 0; p2 < currentEventTracks4V.size(); ++p2) {
                if (vectorEventTracksCharge[nEv][p1] * currentEventTracksCharge[p2] != pairChargeMult) continue;
                if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;

                double qinv = GetQLCMS(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                h_qinv_mix_2l->Fill(qinv);
                h_qinvCor_mix_2l->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
            }
        }
    }
    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);
    
    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
    }
}

// Multithreaded version
void processMixQLCMS(
    Int_t n_tracks,
    int thread_count,
    size_t poolSize,
    Float_t vertexDistance,
    Int_t pairChargeMult,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    TH1D* h_qinv_mix_2l,
    TH1D* h_qinvCor_mix_2l)
{   
    if (pvZ >= vertexDistance) return;
    if (!currentEventTracks4V.empty()){
        if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()){
            if (thread_count == 0) thread_count = 1;
            struct LocalHistograms {
                TH1D* mix;
                TH1D* mix_cor;
        
                LocalHistograms(const TH1D* base_mix, const TH1D* base_mix_cor) {
                    mix      = (TH1D*)base_mix->Clone(TString::Format("%s_clone_%p", base_mix->GetName(), this));
                    mix_cor  = (TH1D*)base_mix_cor->Clone(TString::Format("%s_clone_%p", base_mix_cor->GetName(), this));
                    mix->Reset();
                    mix_cor->Reset();
                }

                ~LocalHistograms() {
                    delete mix;
                    delete mix_cor;
                }
            };
        
            std::vector<std::thread> threads;
            std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        
            for (int i = 0; i < thread_count; ++i) {
                local_hists.push_back(std::make_unique<LocalHistograms>(h_qinv_mix_2l, h_qinvCor_mix_2l));
            }
        
            auto thread_task = [&](size_t start_p1, size_t end_p1, size_t nEv, LocalHistograms* hists) {
                for (size_t p2 = start_p1; p2 < end_p1; ++p2) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
                        if (vectorEventTracksCharge[nEv][p1] * currentEventTracksCharge[p2] != pairChargeMult) continue;
                        if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;
        
                        double qinv = GetQLCMS(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                        hists->mix->Fill(qinv);
                        hists->mix_cor->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
                    }
                }
            };
        
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
                size_t chunk_size = n_tracks / thread_count;
                for (int t = 0; t < thread_count; t++) {
                    size_t start_i = t * chunk_size;
                    size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
            
                    threads.emplace_back(thread_task, start_i, end_i, nEv, local_hists[t].get());
                }
                for (auto& thread : threads) {
                    thread.join();
                }

                threads.clear();

                for (int t = 0; t < thread_count; ++t) {
                    h_qinv_mix_2l->Add(local_hists[t]->mix);
                    h_qinvCor_mix_2l->Add(local_hists[t]->mix_cor);
                    local_hists[t]->mix->Reset();
                    local_hists[t]->mix_cor->Reset();
                }
            }
        }
    }

    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);

    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
    }
}

// Multithreaded version with pT cut
void processMixQLCMS(
    Int_t n_tracks,
    int thread_count,
    size_t poolSize,
    Float_t vertexDistance,
    Int_t pairChargeMult,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    std::vector<std::vector<Float_t>>& vectorEventTracksPt,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    const std::vector<Float_t>& currentEventTracksPt,
    TH1D* h_qinv_mix_2l,
    TH1D* h_qinvCor_mix_2l)
{   
    if (pvZ >= vertexDistance) return;
    if (!currentEventTracks4V.empty()){
        if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()){
            if (thread_count == 0) thread_count = 1;
            struct LocalHistograms {
                TH1D* mix;
                TH1D* mix_cor;
        
                LocalHistograms(const TH1D* base_mix, const TH1D* base_mix_cor) {
                    mix      = (TH1D*)base_mix->Clone(TString::Format("%s_clone_%p", base_mix->GetName(), this));
                    mix_cor  = (TH1D*)base_mix_cor->Clone(TString::Format("%s_clone_%p", base_mix_cor->GetName(), this));
                    mix->Reset();
                    mix_cor->Reset();
                }

                ~LocalHistograms() {
                    delete mix;
                    delete mix_cor;
                }
            };
        
            std::vector<std::thread> threads;
            std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        
            for (int i = 0; i < thread_count; ++i) {
                local_hists.push_back(std::make_unique<LocalHistograms>(h_qinv_mix_2l, h_qinvCor_mix_2l));
            }
        
            auto thread_task = [&](size_t start_p1, size_t end_p1, size_t nEv, LocalHistograms* hists) {
                for (size_t p2 = start_p1; p2 < end_p1; ++p2) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
                        if (vectorEventTracksCharge[nEv][p1] * currentEventTracksCharge[p2] != pairChargeMult) continue;
                        if (vectorEventTracksPt[nEv][p1] <= MIN_PT_CUT || currentEventTracksPt[p2] <= MIN_PT_CUT) continue; // pT cut
                        if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;
        
                        double qinv = GetQLCMS(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);
                        hists->mix->Fill(qinv);
                        hists->mix_cor->Fill(qinv, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
                    }
                }
            };
        
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
                size_t chunk_size = n_tracks / thread_count;
                for (int t = 0; t < thread_count; t++) {
                    size_t start_i = t * chunk_size;
                    size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
            
                    threads.emplace_back(thread_task, start_i, end_i, nEv, local_hists[t].get());
                }
                for (auto& thread : threads) {
                    thread.join();
                }

                threads.clear();

                for (int t = 0; t < thread_count; ++t) {
                    h_qinv_mix_2l->Add(local_hists[t]->mix);
                    h_qinvCor_mix_2l->Add(local_hists[t]->mix_cor);
                    local_hists[t]->mix->Reset();
                    local_hists[t]->mix_cor->Reset();
                }
            }
        }
    }

    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);
    vectorEventTracksPt.push_back(currentEventTracksPt);

    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
        vectorEventTracksPt.clear();
    }
}

// Multithreaded version with pT cut and eta cut, only Corrected hists
void processMixQLCMS(
    Int_t n_tracks,
    int thread_count,
    size_t poolSize,
    Float_t vertexDistance,
    Float_t pvZ,
    std::vector<std::vector<FourVector>>& vectorEventTracks4V,
    std::vector<std::vector<Int_t>>& vectorEventTracksCharge,
    std::vector<std::vector<Float_t>>& vectorEventTracksWeight,
    const std::vector<FourVector>& currentEventTracks4V,
    const std::vector<Int_t>& currentEventTracksCharge,
    const std::vector<Float_t>& currentEventTracksWeight,
    TH1D* h_qinvSSCor_mix_2l,
    TH1D* h_qinvOSCor_mix_2l)
{   
    if (pvZ >= vertexDistance) return;
    if (!currentEventTracks4V.empty()){
        if (vectorEventTracks4V.size() < poolSize && !vectorEventTracks4V.empty()){
            if (thread_count == 0) thread_count = 1;
            struct LocalHistograms {
                TH1D* mix_SS_cor;
                TH1D* mix_OS_cor;
        
                LocalHistograms(const TH1D* base_mix_SS_cor, const TH1D* base_mix_OS_cor) {
                    mix_SS_cor = (TH1D*)base_mix_SS_cor->Clone(TString::Format("%s_clone_%p", base_mix_SS_cor->GetName(), this));
                    mix_OS_cor = (TH1D*)base_mix_OS_cor->Clone(TString::Format("%s_clone_%p", base_mix_OS_cor->GetName(), this));
                    mix_SS_cor->Reset();
                    mix_OS_cor->Reset();
                }

                ~LocalHistograms() {
                    delete mix_SS_cor;
                    delete mix_OS_cor;
                }
            };
        
            std::vector<std::thread> threads;
            std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        
            for (int i = 0; i < thread_count; ++i) {
                local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSSCor_mix_2l, h_qinvOSCor_mix_2l));
            }
        
            auto thread_task = [&](size_t start_p1, size_t end_p1, size_t nEv, LocalHistograms* hists) {
                for (size_t p2 = start_p1; p2 < end_p1; ++p2) {
                    for (size_t p1 = 0; p1 < vectorEventTracks4V[nEv].size(); ++p1) {
                        if (vectorEventTracks4V[nEv][p1].Pt() <= MIN_PT_CUT 
                            || currentEventTracks4V[p2].Pt() <= MIN_PT_CUT) continue; // pT cut
                        if (std::abs(vectorEventTracks4V[nEv][p1].Eta()) >= ABS_ETA_CUT 
                            || std::abs(currentEventTracks4V[p2].Eta()) >= ABS_ETA_CUT) continue; // pT cut
                        if (std::isinf(vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2])) continue;
        
                        double qlcms = GetQLCMS(vectorEventTracks4V[nEv][p1], currentEventTracks4V[p2]);

                        if (vectorEventTracksCharge[nEv][p1]*currentEventTracksCharge[p2] > 0){
                            hists->mix_SS_cor->Fill(qlcms, vectorEventTracksWeight[nEv][p1]*currentEventTracksWeight[p2]);
                        } else {
                            hists->mix_OS_cor->Fill(qlcms, vectorEventTracksWeight[nEv][p1] * currentEventTracksWeight[p2]);
                        }
                    }
                }
            };
        
            for (size_t nEv = 0; nEv < vectorEventTracks4V.size(); ++nEv) {
                size_t chunk_size = n_tracks / thread_count;
                for (int t = 0; t < thread_count; t++) {
                    size_t start_i = t * chunk_size;
                    size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
            
                    threads.emplace_back(thread_task, start_i, end_i, nEv, local_hists[t].get());
                }
                for (auto& thread : threads) {
                    thread.join();
                }

                threads.clear();

                for (int t = 0; t < thread_count; ++t) {
                    h_qinvSSCor_mix_2l->Add(local_hists[t]->mix_SS_cor);
                    h_qinvOSCor_mix_2l->Add(local_hists[t]->mix_OS_cor);
                    local_hists[t]->mix_SS_cor->Reset();
                    local_hists[t]->mix_OS_cor->Reset();
                }
            }
        }
    }

    vectorEventTracks4V.push_back(currentEventTracks4V);
    vectorEventTracksCharge.push_back(currentEventTracksCharge);
    vectorEventTracksWeight.push_back(currentEventTracksWeight);

    if (vectorEventTracks4V.size() >= poolSize) {
        vectorEventTracks4V.clear();
        vectorEventTracksCharge.clear();
        vectorEventTracksWeight.clear();
    }
}

// Process 2D histograms for qt,qz, q0, |q|
// Multithreaded version with pT cut
void processSignalqtqzq0q(
    size_t n_tracks,
    int thread_count,
    const Float_t* trkWeight,
    const Int_t* trkCharge,
    const std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>& currentEventTracks4V,
    TH2D* h_qtqz_SS_cor,
    TH2D* h_qtqz_OS_cor,
    TH2D* h_q0q_SS_cor,
    TH2D* h_q0q_OS_cor)

{
    if (n_tracks <= 1) return; // check if event has 2 or more tracks
    if (thread_count == 0) thread_count = 1;

    // A simple struct to hold the set of histograms for one thread
    struct LocalHistograms {
        TH2D *qtqz_ss_cor, *qtqz_os_cor;
        TH2D *q0q_ss_cor, *q0q_os_cor;

        LocalHistograms(
            const TH2D* b2, const TH2D* b4,
            const TH2D* b6, const TH2D* b8)
        {
            qtqz_ss_cor = (TH2D*)b2->Clone(Form("%s_%p", b2->GetName(), this));
            qtqz_os_cor = (TH2D*)b4->Clone(Form("%s_%p", b4->GetName(), this));

            q0q_ss_cor  = (TH2D*)b6->Clone(Form("%s_%p", b6->GetName(), this));
            q0q_os_cor  = (TH2D*)b8->Clone(Form("%s_%p", b8->GetName(), this));

            qtqz_ss_cor->Reset();
            qtqz_os_cor->Reset();
            q0q_ss_cor->Reset();
            q0q_os_cor->Reset();
        }

        ~LocalHistograms() {
            delete qtqz_ss_cor;
            delete qtqz_os_cor;
            delete q0q_ss_cor;
            delete q0q_os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    // Create a set of local histograms for each thread
    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qtqz_SS_cor, h_qtqz_OS_cor, h_q0q_SS_cor, h_q0q_OS_cor));
    }
    // The task for each thread
    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (std::isinf(trkWeight[p1] * trkWeight[p2])) continue;
                if (currentEventTracks4V[p1].Pt() <= MIN_PT_CUT 
                    || currentEventTracks4V[p2].Pt() <= MIN_PT_CUT) continue; // pT cut
                if (std::abs(currentEventTracks4V[p1].Eta()) >= ABS_ETA_CUT 
                    || std::abs(currentEventTracks4V[p2].Eta()) >= ABS_ETA_CUT) continue; // eTa cut

                auto tp1 = currentEventTracks4V[p1], tp2 = currentEventTracks4V[p2];

                // ---- LCMS boost ----
                auto K = 0.5 * (tp1 + tp2);
                ROOT::Math::BoostZ boost(-K.Pz() / K.E());
                
                tp1 = boost(tp1);
                tp2 = boost(tp2);

                // ---- relative four-momentum ----
                auto q = tp1 - tp2;
                
                double q0 = q.E();
                double qz = q.Pz();
                double qt = std::sqrt(q.Px()*q.Px() + q.Py()*q.Py());
                double qabs = std::sqrt(qt*qt + qz*qz);

                if (trkCharge[p1] * trkCharge[p2] > 0) {
                    hists->qtqz_ss_cor->Fill(qt, qz, trkWeight[p1] * trkWeight[p2]);
                    hists->q0q_ss_cor->Fill(q0, qabs, trkWeight[p1] * trkWeight[p2]);
                } else {
                    hists->qtqz_os_cor->Fill(qt, qz, trkWeight[p1] * trkWeight[p2]);
                    hists->q0q_os_cor->Fill(q0, qabs, trkWeight[p1] * trkWeight[p2]);
                }
            }
        }
    };
    // Distribute the work among threads
    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;

        // Each thread gets a pointer to its own set of histograms
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());

    }
    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Merge the results from local histograms into the main ones
    for (int t = 0; t < thread_count; ++t) {
        h_qtqz_SS_cor->Add(local_hists[t]->qtqz_ss_cor);
        h_qtqz_OS_cor->Add(local_hists[t]->qtqz_os_cor);
        h_q0q_SS_cor->Add(local_hists[t]->q0q_ss_cor);
        h_q0q_OS_cor->Add(local_hists[t]->q0q_os_cor);
    }
}

// Process signal for Delta Eta Delta Phi
// --- Updated Signal Function ---
void processSignalDeltaEtaDeltaPhi(
    int thread_count,
    const EventData& currentEv,
    TH2D* h_DeltaEtaDeltaPhiSS_signal)
{   
    size_t n_tracks = currentEv.tracks.size();
    if (n_tracks <= 1) return;
    if (thread_count == 0) thread_count = 1;

    struct LocalHistograms {
        TH2D* hSigSS; 
        LocalHistograms(const TH2D* baseSSS) {
            hSigSS = (TH2D*)baseSSS->Clone(TString::Format("%s_clone_%p", baseSSS->GetName(), this));
            hSigSS->Reset();
        }
        ~LocalHistograms() { delete hSigSS; }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_DeltaEtaDeltaPhiSS_signal));
    }

    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; p1++) {
            for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                if (currentEv.charges[p1] * currentEv.charges[p2] > 0) { // SS
                    const double dEta = std::abs(currentEv.tracks[p1].Eta() - currentEv.tracks[p2].Eta());
                    const double dPhi = std::abs(TVector2::Phi_mpi_pi(currentEv.tracks[p1].Phi() - currentEv.tracks[p2].Phi()));
                    
                    // Study range check
                    if (dEta < 0.15 && dPhi < 0.15) {
                        hists->hSigSS->Fill(dEta, dPhi);
                    }
                }
            }
        }
    };

    size_t chunk_size = n_tracks / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk_size;
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());
    }
    for (auto& thread : threads) thread.join();
    for (int t = 0; t < thread_count; ++t) h_DeltaEtaDeltaPhiSS_signal->Add(local_hists[t]->hSigSS);
}

// --- Updated Mixing Function ---
void processMixDeltaEtaDeltaPhi(
    int thread_count,
    const EventData& currentEv,
    const std::vector<EventData>& zBucketPool,
    TH2D* h_DeltaEtaDeltaPhiSS_mix)
{   
    if (currentEv.tracks.empty() || zBucketPool.empty()) return;
    if (thread_count == 0) thread_count = 1;

    struct LocalHistograms {
        TH2D* mixSS;
        LocalHistograms(const TH2D* baseSS) {
            mixSS = (TH2D*)baseSS->Clone(TString::Format("%s_clone_%p", baseSS->GetName(), this));
            mixSS->Reset();
        }
        ~LocalHistograms() { delete mixSS; }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_DeltaEtaDeltaPhiSS_mix));
    }

    auto thread_task = [&](size_t start_p1, size_t end_p1, LocalHistograms* hists) {
        for (size_t p1 = start_p1; p1 < end_p1; ++p1) {
            for (const auto& pastEv : zBucketPool) {
                for (size_t p2 = 0; p2 < pastEv.tracks.size(); ++p2) {
                    if (currentEv.charges[p1] * pastEv.charges[p2] > 0) { // SS
                        const double dEta = std::abs(currentEv.tracks[p1].Eta() - pastEv.tracks[p2].Eta());
                        const double dPhi = std::abs(TVector2::Phi_mpi_pi(currentEv.tracks[p1].Phi() - pastEv.tracks[p2].Phi()));
                        
                        if (dEta < 0.15 && dPhi < 0.15) {
                            hists->mixSS->Fill(dEta, dPhi);
                        }
                    }
                }
            }
        }
    };

    size_t chunk_size = currentEv.tracks.size() / thread_count;
    for (int t = 0; t < thread_count; t++) {
        size_t start_i = t * chunk_size;
        size_t end_i = (t == thread_count - 1) ? currentEv.tracks.size() : (t + 1) * chunk_size;
        threads.emplace_back(thread_task, start_i, end_i, local_hists[t].get());
    }
    for (auto& thread : threads) thread.join();
    for (int t = 0; t < thread_count; ++t) h_DeltaEtaDeltaPhiSS_mix->Add(local_hists[t]->mixSS);
}

#endif