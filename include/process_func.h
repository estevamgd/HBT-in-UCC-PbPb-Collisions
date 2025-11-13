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
#include "TMath.h"
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>
#include <chrono>
#include <string>
#include "data_func.h"
#include <TSystem.h>


void processSignal(
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

void processSignal(
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

            ss->Reset();
            ss_cor->Reset();
            os->Reset();
            os_cor->Reset();
        }
        ~LocalHistograms() {
            delete ss;
            delete ss_cor;
            delete os;
            delete os_cor;
        }
    };

    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<LocalHistograms>> local_hists;

    for (int i = 0; i < thread_count; ++i) {
        local_hists.push_back(std::make_unique<LocalHistograms>(h_qinvSS_signal_2l, h_qinvSSCor_signal_2l, h_qinvOS_signal_2l, h_qinvOSCor_signal_2l));
    }
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


void processMix(
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


void processMix(
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
#endif