#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <thread>
#include <string>
#include <memory>
#include <algorithm>

// ----- ROOT Headers (Must be included before my_func.h) -----
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"

// ----- Custom Headers -----
#include "../include/my_func.h"

// ----- Struct for Input Configuration -----
struct InputParams {
    std::string fileInput;
    std::string treeInput;
};

// ----- Struct for Analysis Parameters and Toggles -----
struct AnalysisParams {
    // Event Selection
    double selVarMoreeq = 0.0;
    double selVarLess = 100.0;
    ControlVar selectionVarType = ControlVar::CENT;

    // Event Cuts (Pile-up)
    bool usePileUpCut = false;
    double pileUpXStart = 675.0;
    double pileUpYStart = 600000.0;
    double pileUpXEnd   = 3900.0;
    double pileUpYEnd   = 0.0;
    double pileUpSlope = 0.0;     // Calculated in constructor
    double pileUpIntercept = 0.0; // Calculated in constructor

    // Processing Limits (Toggles)
    bool useTestLimitSig = false;
    int testLimitSig = -1;
    bool useTestLimitMix = false;
    int testLimitMix = -1;

    // Track Filters (Toggles + Values)
    bool usePtMin = true;
    float ptMin = 0.5; 
    bool usePtMax = false;
    float ptMax = 5.0;
    bool useEtaCut = true;
    float etaCut = 0.95;
    bool usePixHit = true;
    float pixHit = 1.0;
    bool useNvalidHits = false;
    float nvalidHits = 0.0;
    
    // Pair Filters (Split Tracks)
    bool useSplitCut = false;
    double splitCosCut = 0.99996; // Reject if cos(a) > this
    double splitDeltaPtCut = 0.04; // AND deltaPt < this

    // Vertex and Mixing Config
    float vertexDistance = 15.0; 
    int poolSizeInt = 10;
    float zBinWidth = 2.0;

    // Constructor to calculate PileUp parameters automatically
    AnalysisParams() {
        pileUpSlope = (pileUpYEnd - pileUpYStart) / (pileUpXEnd - pileUpXStart);
        pileUpIntercept = pileUpYEnd - pileUpSlope * pileUpXEnd;
    }
};

// ----- Process Class for Analysis Management -----
class Process {
private:
    InputParams inputs;
    AnalysisParams cfg;

    // --- Helper Math Functions (Inline) ---
    double GetQ(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& p1, 
                const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& p2) {
        auto q = p1 - p2;
        return std::sqrt(std::abs(q.M2())); 
    }

    double GetQLCMS(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& p1, 
                    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& p2) {
        auto K = 0.5 * (p1 + p2);
        ROOT::Math::BoostZ boost(-K.Pz() / K.E());
        auto p1_lcms = boost(p1);
        auto p2_lcms = boost(p2);
        auto q = p1_lcms - p2_lcms;
        return std::sqrt(q.Px()*q.Px() + q.Py()*q.Py() + q.Pz()*q.Pz()); 
    }

    // --- Cut Functions ---
    
    // Returns TRUE if event is good (passes pile-up cut)
    bool IsGoodEventPileUp(float zdcSum, float hfSum) {
        if (!cfg.usePileUpCut) return true;
        // Formula: zdcSum < (slope * HF + intercept)
        return zdcSum < (cfg.pileUpSlope * hfSum + cfg.pileUpIntercept);
    }

    // Returns TRUE if pair is good (NOT a split track)
    bool IsGoodPairSplit(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& vec1, 
                         const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& vec2) {
        if (!cfg.useSplitCut) return true;

        Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
        Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());

        // Reject if very collinear (high cos) AND similar Pt (low deltaPt)
        if (cosa > cfg.splitCosCut && deltapt < cfg.splitDeltaPtCut) return false;
        
        return true;
    }

    bool PassesTrackCuts(float pt, float eta, float nPixLayers, float nValidHits, float weight) const {
        if (cfg.usePtMin && (pt < cfg.ptMin)) return false;
        if (cfg.usePtMax && (pt > cfg.ptMax)) return false;
        if (cfg.useEtaCut && (std::abs(eta) > cfg.etaCut)) return false;
        if (cfg.usePixHit && (nPixLayers < cfg.pixHit)) return false;
        if (cfg.useNvalidHits && (nValidHits <= cfg.nvalidHits)) return false;
        if (std::isinf(weight)) return false;
        return true;
    }

    // =================================================================================
    // Internal Process Functions 
    // =================================================================================

    // --- DeltaEta-DeltaPhi Processes ---
    void processSignalDeltaEtaDeltaPhi(int thread_count, const EventData& currentEv, TH2D* hSigSS) {
        size_t n_tracks = currentEv.tracks.size();
        if (n_tracks <= 1) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH2D* h;
            LocalHistograms(const TH2D* base) {
                h = (TH2D*)base->Clone(TString::Format("%s_clone_%p", base->GetName(), this));
                h->Reset();
            }
            ~LocalHistograms() { delete h; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hSigSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; p1++) {
                for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                    if (currentEv.charges[p1] * currentEv.charges[p2] > 0) { // SS Only
                        
                        // --- Split Cut Check ---
                        if (cfg.useSplitCut) {
                            if (!IsGoodPairSplit(currentEv.tracks[p1], currentEv.tracks[p2])) continue;
                        }

                        const double dEta = std::abs(currentEv.tracks[p1].Eta() - currentEv.tracks[p2].Eta());
                        const double dPhi = std::abs(TVector2::Phi_mpi_pi(currentEv.tracks[p1].Phi() - currentEv.tracks[p2].Phi()));
                        if (dEta < 0.1 && dPhi < 0.1) {
                            lh->h->Fill(dEta, dPhi, currentEv.weights[p1] * currentEv.weights[p2]);
                        }
                    }
                }
            }
        };

        size_t chunk = n_tracks / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) hSigSS->Add(lh->h);
    }

    void processMixDeltaEtaDeltaPhi(int thread_count, const EventData& currentEv, const std::vector<EventData>& pool, TH2D* hMixSS) {
        if (currentEv.tracks.empty() || pool.empty()) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH2D* h;
            LocalHistograms(const TH2D* base) {
                h = (TH2D*)base->Clone(TString::Format("%s_clone_%p", base->GetName(), this));
                h->Reset();
            }
            ~LocalHistograms() { delete h; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hMixSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; ++p1) {
                for (const auto& pastEv : pool) {
                    for (size_t p2 = 0; p2 < pastEv.tracks.size(); ++p2) {
                        if (currentEv.charges[p1] * pastEv.charges[p2] > 0) { // SS Only
                            
                            // --- Split Cut Check ---
                            if (cfg.useSplitCut) {
                                if (!IsGoodPairSplit(currentEv.tracks[p1], pastEv.tracks[p2])) continue;
                            }

                            const double dEta = std::abs(currentEv.tracks[p1].Eta() - pastEv.tracks[p2].Eta());
                            const double dPhi = std::abs(TVector2::Phi_mpi_pi(currentEv.tracks[p1].Phi() - pastEv.tracks[p2].Phi()));
                            if (dEta < 0.1 && dPhi < 0.1) {
                                lh->h->Fill(dEta, dPhi, currentEv.weights[p1] * pastEv.weights[p2]);
                            }
                        }
                    }
                }
            }
        };

        size_t chunk = currentEv.tracks.size() / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? currentEv.tracks.size() : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) hMixSS->Add(lh->h);
    }

    // --- Qinv Processes ---
    void processSignalQinv(int thread_count, const EventData& currentEv, TH1D* hSigSS) {
        size_t n_tracks = currentEv.tracks.size();
        if (n_tracks <= 1) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH1D* h;
            LocalHistograms(const TH1D* base) {
                h = (TH1D*)base->Clone(TString::Format("%s_clone_%p", base->GetName(), this));
                h->Reset();
            }
            ~LocalHistograms() { delete h; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hSigSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; p1++) {
                for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                    if (currentEv.charges[p1] * currentEv.charges[p2] > 0) { // SS Only
                        
                        if (cfg.useSplitCut) {
                            if (!IsGoodPairSplit(currentEv.tracks[p1], currentEv.tracks[p2])) continue;
                        }

                        double q = GetQ(currentEv.tracks[p1], currentEv.tracks[p2]);
                        lh->h->Fill(q, currentEv.weights[p1] * currentEv.weights[p2]);
                    }
                }
            }
        };

        size_t chunk = n_tracks / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) hSigSS->Add(lh->h);
    }

    void processMixQinv(int thread_count, const EventData& currentEv, const std::vector<EventData>& pool, TH1D* hMixSS) {
        if (currentEv.tracks.empty() || pool.empty()) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH1D* h;
            LocalHistograms(const TH1D* base) {
                h = (TH1D*)base->Clone(TString::Format("%s_clone_%p", base->GetName(), this));
                h->Reset();
            }
            ~LocalHistograms() { delete h; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hMixSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; ++p1) {
                for (const auto& pastEv : pool) {
                    for (size_t p2 = 0; p2 < pastEv.tracks.size(); ++p2) {
                        if (currentEv.charges[p1] * pastEv.charges[p2] > 0) { // SS Only
                            
                            if (cfg.useSplitCut) {
                                if (!IsGoodPairSplit(currentEv.tracks[p1], pastEv.tracks[p2])) continue;
                            }

                            double q = GetQ(currentEv.tracks[p1], pastEv.tracks[p2]);
                            lh->h->Fill(q, currentEv.weights[p1] * pastEv.weights[p2]);
                        }
                    }
                }
            }
        };

        size_t chunk = currentEv.tracks.size() / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? currentEv.tracks.size() : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) hMixSS->Add(lh->h);
    }

    // --- QLCMS Processes ---
    void processSignalQLCMS(int thread_count, const EventData& currentEv, TH1D* hSigSS) {
        size_t n_tracks = currentEv.tracks.size();
        if (n_tracks <= 1) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH1D* h;
            LocalHistograms(const TH1D* base) {
                h = (TH1D*)base->Clone(TString::Format("%s_clone_%p", base->GetName(), this));
                h->Reset();
            }
            ~LocalHistograms() { delete h; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hSigSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; p1++) {
                for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                    if (currentEv.charges[p1] * currentEv.charges[p2] > 0) { // SS Only
                        
                        if (cfg.useSplitCut) {
                            if (!IsGoodPairSplit(currentEv.tracks[p1], currentEv.tracks[p2])) continue;
                        }

                        double q = GetQLCMS(currentEv.tracks[p1], currentEv.tracks[p2]);
                        lh->h->Fill(q, currentEv.weights[p1] * currentEv.weights[p2]);
                    }
                }
            }
        };

        size_t chunk = n_tracks / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) hSigSS->Add(lh->h);
    }

    void processMixQLCMS(int thread_count, const EventData& currentEv, const std::vector<EventData>& pool, TH1D* hMixSS) {
        if (currentEv.tracks.empty() || pool.empty()) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH1D* h;
            LocalHistograms(const TH1D* base) {
                h = (TH1D*)base->Clone(TString::Format("%s_clone_%p", base->GetName(), this));
                h->Reset();
            }
            ~LocalHistograms() { delete h; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hMixSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; ++p1) {
                for (const auto& pastEv : pool) {
                    for (size_t p2 = 0; p2 < pastEv.tracks.size(); ++p2) {
                        if (currentEv.charges[p1] * pastEv.charges[p2] > 0) { // SS Only
                            
                            if (cfg.useSplitCut) {
                                if (!IsGoodPairSplit(currentEv.tracks[p1], pastEv.tracks[p2])) continue;
                            }

                            double q = GetQLCMS(currentEv.tracks[p1], pastEv.tracks[p2]);
                            lh->h->Fill(q, currentEv.weights[p1] * pastEv.weights[p2]);
                        }
                    }
                }
            }
        };

        size_t chunk = currentEv.tracks.size() / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? currentEv.tracks.size() : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) hMixSS->Add(lh->h);
    }

    
    // --- QtQzQ0Q Processes (2D) ---
    void processSignalQtQzQ0Q(int thread_count, const EventData& currentEv, TH2D* hQtQzSS, TH2D* hQ0QSS) {
        size_t n_tracks = currentEv.tracks.size();
        if (n_tracks <= 1) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH2D *qtqz, *q0q;
            LocalHistograms(const TH2D* base1, const TH2D* base2) {
                qtqz = (TH2D*)base1->Clone(TString::Format("%s_clone_%p", base1->GetName(), this));
                q0q = (TH2D*)base2->Clone(TString::Format("%s_clone_%p", base2->GetName(), this));
                qtqz->Reset(); q0q->Reset();
            }
            ~LocalHistograms() { delete qtqz; delete q0q; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hQtQzSS, hQ0QSS));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; p1++) {
                for (size_t p2 = p1 + 1; p2 < n_tracks; p2++) {
                    if (currentEv.charges[p1] * currentEv.charges[p2] > 0) { // SS Only
                        
                        if (cfg.useSplitCut) {
                            if (!IsGoodPairSplit(currentEv.tracks[p1], currentEv.tracks[p2])) continue;
                        }

                        auto tp1 = currentEv.tracks[p1];
                        auto tp2 = currentEv.tracks[p2];
                        auto K = 0.5 * (tp1 + tp2);
                        ROOT::Math::BoostZ boost(-K.Pz() / K.E());
                        auto tp1_lcms = boost(tp1);
                        auto tp2_lcms = boost(tp2);
                        auto q = tp1_lcms - tp2_lcms;

                        double q0 = q.E();
                        double qz = q.Pz();
                        double qt = std::sqrt(q.Px()*q.Px() + q.Py()*q.Py());
                        double q0_2 = q0 * q0;
                        double qz_2 = qz * qz;
                        double qt_2 = qt * qt;
                        double qabs_2 = qt_2 + qz_2;

                        lh->qtqz->Fill(qt_2, qz_2, currentEv.weights[p1] * currentEv.weights[p2]);
                        lh->q0q->Fill(q0_2, qabs_2, currentEv.weights[p1] * currentEv.weights[p2]);
                    }
                }
            }
        };

        size_t chunk = n_tracks / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) {
            hQtQzSS->Add(lh->qtqz);
            hQ0QSS->Add(lh->q0q);
        }
    }

    void processMixQtQzQ0Q(int thread_count, const EventData& currentEv, const std::vector<EventData>& pool, TH2D* hQtQzMix, TH2D* hQ0QMix) {
        if (currentEv.tracks.empty() || pool.empty()) return;
        if (thread_count == 0) thread_count = 1;

        struct LocalHistograms {
            TH2D *qtqz, *q0q;
            LocalHistograms(const TH2D* base1, const TH2D* base2) {
                qtqz = (TH2D*)base1->Clone(TString::Format("%s_clone_%p", base1->GetName(), this));
                q0q = (TH2D*)base2->Clone(TString::Format("%s_clone_%p", base2->GetName(), this));
                qtqz->Reset(); q0q->Reset();
            }
            ~LocalHistograms() { delete qtqz; delete q0q; }
        };

        std::vector<std::thread> threads;
        std::vector<std::unique_ptr<LocalHistograms>> local_hists;
        for (int i = 0; i < thread_count; ++i) local_hists.push_back(std::make_unique<LocalHistograms>(hQtQzMix, hQ0QMix));

        auto thread_task = [&](size_t start, size_t end, LocalHistograms* lh) {
            for (size_t p1 = start; p1 < end; ++p1) {
                for (const auto& pastEv : pool) {
                    for (size_t p2 = 0; p2 < pastEv.tracks.size(); ++p2) {
                        if (currentEv.charges[p1] * pastEv.charges[p2] > 0) { // SS Only

                            if (cfg.useSplitCut) {
                                if (!IsGoodPairSplit(currentEv.tracks[p1], pastEv.tracks[p2])) continue;
                            }

                            auto tp1 = currentEv.tracks[p1];
                            auto tp2 = pastEv.tracks[p2];
                            auto K = 0.5 * (tp1 + tp2);
                            ROOT::Math::BoostZ boost(-K.Pz() / K.E());
                            auto tp1_lcms = boost(tp1);
                            auto tp2_lcms = boost(tp2);
                            auto q = tp1_lcms - tp2_lcms;

                            double q0 = q.E();
                            double qz = q.Pz();
                            double qt = std::sqrt(q.Px()*q.Px() + q.Py()*q.Py());
                            double q0_2 = q0 * q0;
                            double qz_2 = qz * qz;
                            double qt_2 = qt * qt;
                            double qabs_2 = qt_2 + qz_2;

                            double pairWeight = currentEv.weights[p1] * pastEv.weights[p2];
                            lh->qtqz->Fill(qt_2, qz_2, pairWeight);
                            lh->q0q->Fill(q0_2, qabs_2, pairWeight);
                        }
                    }
                }
            }
        };

        size_t chunk = currentEv.tracks.size() / thread_count;
        for (int t = 0; t < thread_count; t++) {
            size_t start = t * chunk;
            size_t end = (t == thread_count - 1) ? currentEv.tracks.size() : (t + 1) * chunk;
            threads.emplace_back(thread_task, start, end, local_hists[t].get());
        }
        for (auto& t : threads) t.join();
        for (auto& lh : local_hists) {
            hQtQzMix->Add(lh->qtqz);
            hQ0QMix->Add(lh->q0q);
        }
    }

public:
    Process(InputParams in, AnalysisParams out) : inputs(in), cfg(out) {}

    // ----- DeltaEta-DeltaPhi Analysis (2D) -----
    void buildDeltaEtaDeltaPhi() 
    {
        // ----- INITIALIZATION & DATA LOADING -----
        TFile *fr = nullptr; TTree *t = nullptr;
        getFileTree(inputs.fileInput.c_str(), inputs.treeInput.c_str(), fr, t);

        Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
        Float_t pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize],
                trkNpixLayers[maxSize], trkNvalidHits[maxSize], pionMass = 0.13957039;
        Float_t HFsumET, zdcSum; 

        t->SetBranchAddress("HFsumET", &HFsumET);
        t->SetBranchAddress("zdcSum", &zdcSum);
        t->SetBranchAddress("Ntrk", &Ntrk);
        t->SetBranchAddress("hiBin", &hiBin);
        t->SetBranchAddress("pvZ", &pvZ); 
        t->SetBranchAddress("trkNpixLayers", trkNpixLayers);
        t->SetBranchAddress("trkCharge", trkCharge);
        t->SetBranchAddress("trkWeight", trkWeight);
        t->SetBranchAddress("trkPt", trkPt);
        t->SetBranchAddress("trkEta", trkEta);
        t->SetBranchAddress("trkPhi", trkPhi);
        t->SetBranchAddress("trkNvalidHits", trkNvalidHits);

        // ----- SELECTION VARIABLE SETUP -----
        double hiBinProxy, NtrkProxy, HFsumETProxy;
        double* selectionVar;
        double displayMoreeq = cfg.selVarMoreeq, displayLess = cfg.selVarLess;
        const char* selectionVarName;

        if (cfg.selectionVarType == ControlVar::CENT) {
            selectionVar = &hiBinProxy; cfg.selVarMoreeq *= 2; cfg.selVarLess *= 2; selectionVarName = "CENT";
        } else if (cfg.selectionVarType == ControlVar::MULT) {
            selectionVar = &NtrkProxy; selectionVarName = "MULT";
        } else {
            selectionVar = &HFsumETProxy; selectionVarName = "CENTHF";
        }

        // ----- BINS & POOL SETUP -----
        int nZBins = static_cast<int>(std::ceil((cfg.vertexDistance * 2.0) / cfg.zBinWidth));
        std::vector<std::vector<EventData>> pool(nZBins);
        size_t poolSize = static_cast<size_t>(cfg.poolSizeInt);

        // ----- HISTOGRAMS & TIMING SETUP -----
        double nnscale = 100, x0 = 0., xt = 0.1;
        double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;
        
        TH2D* hSigSS = new TH2D("hSigSS", ";#Delta#eta;#Delta#varphi;#Pairs", nnscale, x0, xt, nnscale, x0, xt);
        TH2D* hMixSS = new TH2D("hMixSS", ";#Delta#eta;#Delta#varphi;#Pairs", nnscale, x0, xt, nnscale, x0, xt);
        hSigSS->Sumw2(); hMixSS->Sumw2();

        Long64_t nentries = t->GetEntries();
        int processedEventsSig = 0, processedEventsMix = 0, thread_count = std::thread::hardware_concurrency();

        // ----- MAIN EVENT LOOP -----
        auto start_full = std::chrono::high_resolution_clock::now();
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            if (i % 500 == 0) printProgressBar(i + 1, nentries);

            hiBinProxy = (double)hiBin; NtrkProxy = (double)Ntrk; HFsumETProxy = (double)HFsumET;

            // ----- Processing Limit Check -----
            bool sigLimitReached = cfg.useTestLimitSig && (processedEventsSig >= cfg.testLimitSig);
            bool mixLimitReached = cfg.useTestLimitMix && (processedEventsMix >= cfg.testLimitMix);
            if (sigLimitReached || mixLimitReached) break;

            // ----- Event-level filters -----
            if (!(*selectionVar >= cfg.selVarMoreeq && *selectionVar < cfg.selVarLess)) continue;
            if (std::abs(pvZ) >= cfg.vertexDistance) continue;
            
            // --- Pile-Up Cut ---
            if (!IsGoodEventPileUp(zdcSum, HFsumET)) continue;

            // ----- Z-Bin calculation -----
            int zBin = (int)((pvZ + cfg.vertexDistance) / cfg.zBinWidth);
            if (zBin < 0) zBin = 0; if (zBin >= nZBins) zBin = nZBins - 1;

            EventData currentEv;
            for (int j = 0; j < Ntrk; j++) {
                if (!PassesTrackCuts(trkPt[j], trkEta[j], trkNpixLayers[j], trkNvalidHits[j], trkWeight[j])) continue;

                currentEv.tracks.push_back(ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass));
                currentEv.weights.push_back(trkWeight[j]);
                currentEv.charges.push_back(trkCharge[j]);
            }

            if (currentEv.tracks.size() < 2) continue;
            processedEventsSig++;

            // ----- Mixing processing -----
            auto start_mix_lap = std::chrono::high_resolution_clock::now();
            if (!pool[zBin].empty()) {
                processMixDeltaEtaDeltaPhi(thread_count, currentEv, pool[zBin], hMixSS);
                processedEventsMix++;
            }
            auto end_mix_lap = std::chrono::high_resolution_clock::now();
            duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();

            // ----- Signal processing -----
            auto start_signal_lap = std::chrono::high_resolution_clock::now();
            processSignalDeltaEtaDeltaPhi(thread_count, currentEv, hSigSS);
            auto end_signal_lap = std::chrono::high_resolution_clock::now();
            duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();

            // ----- Pool Update -----
            pool[zBin].push_back(currentEv);
            if (pool[zBin].size() > poolSize) pool[zBin].erase(pool[zBin].begin());
        }
        auto end_full = std::chrono::high_resolution_clock::now();
        duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

        // ----- CORRELATION CALCULATION -----
        TH2D* hRatio = (TH2D*)hSigSS->Clone("hRatioSS");
        hRatio->Divide(hMixSS);
        
        double nSig = hSigSS->GetEntries();
        double nMix = hMixSS->GetEntries();
        double scale = (nSig > 0) ? (nMix / nSig) : 1.0;
        hRatio->Scale(scale);
        hRatio->SetTitle("Correlation Function C(#Delta#eta, #Delta#varphi);#Delta#eta;#Delta#varphi;C");

        // ----- SAVING -----
        char prefix[256];
        sprintf(prefix, "DeltaEtaDeltaPhi_eta-abs%g_%s_%g-%g", cfg.etaCut, selectionVarName, displayMoreeq, displayLess);
        
        TCanvas *cCorr = new TCanvas("cCorr", "Correlation", 1200, 800);
        cCorr->cd(); 
        hRatio->Draw("COLZ");
        drawCMSHeaders("#bf{CMS} #it{Work in Progress}","PbPb 2.76 TeV | Correlation SS");

        save_benchmark_chrono({duration_full, duration_signal, duration_mix}, {"Total", "Signal", "Mix"}, "benchmarks", prefix, processedEventsSig, processedEventsMix);
        
        TH2D* hists[] = {hSigSS, hMixSS, hRatio};
        save_histograms2d(hists, 3, "./data/signal_mix/", prefix, displayMoreeq, displayLess);
        
        save_canvas_images(new TCanvas*[1]{cCorr}, 1, "./imgs/test/signal_mix/", prefix, "png");
        save_canvas_images(new TCanvas*[1]{cCorr}, 1, "./imgs/test/signal_mix/", prefix, "pdf");

        AnalysisLog::instance().save("./logs", "sig_DeltaEtaDeltaPhi");
        
        delete cCorr;
        delete hSigSS;
        delete hMixSS;
        delete hRatio;
        if(fr) { fr->Close(); delete fr; }
    }

    // ----- Q-Invariant Analysis (1D) -----
    void buildQinv() 
    {
        // ----- INITIALIZATION & DATA LOADING -----
        TFile *fr = nullptr; TTree *t = nullptr;
        getFileTree(inputs.fileInput.c_str(), inputs.treeInput.c_str(), fr, t);

        Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
        Float_t pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize],
                trkNpixLayers[maxSize], trkNvalidHits[maxSize], pionMass = 0.13957039;
        Float_t HFsumET, zdcSum; 

        t->SetBranchAddress("HFsumET", &HFsumET);
        t->SetBranchAddress("zdcSum", &zdcSum);
        t->SetBranchAddress("Ntrk", &Ntrk);
        t->SetBranchAddress("hiBin", &hiBin);
        t->SetBranchAddress("pvZ", &pvZ); 
        t->SetBranchAddress("trkNpixLayers", trkNpixLayers);
        t->SetBranchAddress("trkCharge", trkCharge);
        t->SetBranchAddress("trkWeight", trkWeight);
        t->SetBranchAddress("trkPt", trkPt);
        t->SetBranchAddress("trkEta", trkEta);
        t->SetBranchAddress("trkPhi", trkPhi);
        t->SetBranchAddress("trkNvalidHits", trkNvalidHits);

        // ----- SELECTION VARIABLE SETUP -----
        double hiBinProxy, NtrkProxy, HFsumETProxy;
        double* selectionVar;
        double displayMoreeq = cfg.selVarMoreeq, displayLess = cfg.selVarLess;
        const char* selectionVarName;

        if (cfg.selectionVarType == ControlVar::CENT) {
            selectionVar = &hiBinProxy; cfg.selVarMoreeq *= 2; cfg.selVarLess *= 2; selectionVarName = "CENT";
        } else if (cfg.selectionVarType == ControlVar::MULT) {
            selectionVar = &NtrkProxy; selectionVarName = "MULT";
        } else {
            selectionVar = &HFsumETProxy; selectionVarName = "CENTHF";
        }

        // ----- BINS & POOL SETUP -----
        int nZBins = static_cast<int>(std::ceil((cfg.vertexDistance * 2.0) / cfg.zBinWidth));
        std::vector<std::vector<EventData>> pool(nZBins);
        size_t poolSize = static_cast<size_t>(cfg.poolSizeInt);

        // ----- HISTOGRAMS & TIMING SETUP -----
        double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;
        double nnscale = 10000, x0 = 0., xt = 10.;
        
        TH1D* hSigSS = new TH1D("hSigSS", ";q_{inv} [GeV];#Pairs", nnscale, x0, xt);
        TH1D* hMixSS = new TH1D("hMixSS", ";q_{inv} [GeV];#Pairs", nnscale, x0, xt);
        hSigSS->Sumw2(); hMixSS->Sumw2();

        Long64_t nentries = t->GetEntries();
        int processedEventsSig = 0, processedEventsMix = 0, thread_count = std::thread::hardware_concurrency();

        // ----- MAIN EVENT LOOP -----
        auto start_full = std::chrono::high_resolution_clock::now();
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            if (i % 500 == 0) printProgressBar(i + 1, nentries);

            hiBinProxy = (double)hiBin; NtrkProxy = (double)Ntrk; HFsumETProxy = (double)HFsumET;

            // ----- Processing Limit Check -----
            bool sigLimitReached = cfg.useTestLimitSig && (processedEventsSig >= cfg.testLimitSig);
            bool mixLimitReached = cfg.useTestLimitMix && (processedEventsMix >= cfg.testLimitMix);
            if (sigLimitReached || mixLimitReached) break;

            // ----- Event-level filters -----
            if (!(*selectionVar >= cfg.selVarMoreeq && *selectionVar < cfg.selVarLess)) continue;
            if (std::abs(pvZ) >= cfg.vertexDistance) continue;
            
            // --- Pile-Up Cut ---
            if (!IsGoodEventPileUp(zdcSum, HFsumET)) continue;

            // ----- Z-Bin calculation -----
            int zBin = (int)((pvZ + cfg.vertexDistance) / cfg.zBinWidth);
            if (zBin < 0) zBin = 0; if (zBin >= nZBins) zBin = nZBins - 1;

            EventData currentEv;
            for (int j = 0; j < Ntrk; j++) {
                if (!PassesTrackCuts(trkPt[j], trkEta[j], trkNpixLayers[j], trkNvalidHits[j], trkWeight[j])) continue;

                currentEv.tracks.push_back(ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass));
                currentEv.weights.push_back(trkWeight[j]);
                currentEv.charges.push_back(trkCharge[j]);
            }

            if (currentEv.tracks.size() < 2) continue;
            processedEventsSig++;

            // ----- Mixing processing -----
            auto start_mix_lap = std::chrono::high_resolution_clock::now();
            if (!pool[zBin].empty()) {
                processMixQinv(thread_count, currentEv, pool[zBin], hMixSS);
                processedEventsMix++;
            }
            auto end_mix_lap = std::chrono::high_resolution_clock::now();
            duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();

            // ----- Signal processing -----
            auto start_signal_lap = std::chrono::high_resolution_clock::now();
            processSignalQinv(thread_count, currentEv, hSigSS);
            auto end_signal_lap = std::chrono::high_resolution_clock::now();
            duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();

            // ----- Pool Update -----
            pool[zBin].push_back(currentEv);
            if (pool[zBin].size() > poolSize) pool[zBin].erase(pool[zBin].begin());
        }
        auto end_full = std::chrono::high_resolution_clock::now();
        duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

        // ----- SAVING -----
        char prefix[256];
        sprintf(prefix, "qinv_eta-abs%g_%s_%g-%g", cfg.etaCut, selectionVarName, displayMoreeq, displayLess);
        TH1D *hists[] = {hSigSS, hMixSS};
        save_benchmark_chrono({duration_full, duration_signal, duration_mix}, {"Total", "Signal", "Mix"}, "benchmarks", prefix, processedEventsSig, processedEventsMix);
        save_histograms(hists, 2, "./data/signal_mix/", prefix, displayMoreeq, displayLess);
        AnalysisLog::instance().save("./logs", "sig_qinv");
        close_program(nullptr, 0, hists, 2, nullptr, 0, fr);
    }

    // ----- Q-LCMS Analysis (1D) -----
    void buildQlcms() 
    {
        // ----- INITIALIZATION & DATA LOADING -----
        TFile *fr = nullptr; TTree *t = nullptr;
        getFileTree(inputs.fileInput.c_str(), inputs.treeInput.c_str(), fr, t);

        Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
        Float_t pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize],
                trkNpixLayers[maxSize], trkNvalidHits[maxSize], pionMass = 0.13957039;
        Float_t HFsumET, zdcSum; 

        t->SetBranchAddress("HFsumET", &HFsumET);
        t->SetBranchAddress("zdcSum", &zdcSum);
        t->SetBranchAddress("Ntrk", &Ntrk);
        t->SetBranchAddress("hiBin", &hiBin);
        t->SetBranchAddress("pvZ", &pvZ); 
        t->SetBranchAddress("trkNpixLayers", trkNpixLayers);
        t->SetBranchAddress("trkCharge", trkCharge);
        t->SetBranchAddress("trkWeight", trkWeight);
        t->SetBranchAddress("trkPt", trkPt);
        t->SetBranchAddress("trkEta", trkEta);
        t->SetBranchAddress("trkPhi", trkPhi);
        t->SetBranchAddress("trkNvalidHits", trkNvalidHits);

        // ----- SELECTION VARIABLE SETUP -----
        double hiBinProxy, NtrkProxy, HFsumETProxy;
        double* selectionVar;
        double displayMoreeq = cfg.selVarMoreeq, displayLess = cfg.selVarLess;
        const char* selectionVarName;

        if (cfg.selectionVarType == ControlVar::CENT) {
            selectionVar = &hiBinProxy; cfg.selVarMoreeq *= 2; cfg.selVarLess *= 2; selectionVarName = "CENT";
        } else if (cfg.selectionVarType == ControlVar::MULT) {
            selectionVar = &NtrkProxy; selectionVarName = "MULT";
        } else {
            selectionVar = &HFsumETProxy; selectionVarName = "CENTHF";
        }

        // ----- BINS & POOL SETUP -----
        int nZBins = static_cast<int>(std::ceil((cfg.vertexDistance * 2.0) / cfg.zBinWidth));
        std::vector<std::vector<EventData>> pool(nZBins);
        size_t poolSize = static_cast<size_t>(cfg.poolSizeInt);

        // ----- HISTOGRAMS & TIMING SETUP -----
        double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;
        double nnscale = 10000, x0 = 0., xt = 10.;

        TH1D* hSigSS = new TH1D("hSigSS", ";q_{LCMS} [GeV];#Pairs", nnscale, x0, xt);
        TH1D* hMixSS = new TH1D("hMixSS", ";q_{LCMS} [GeV];#Pairs", nnscale, x0, xt);
        hSigSS->Sumw2(); hMixSS->Sumw2();

        Long64_t nentries = t->GetEntries();
        int processedEventsSig = 0, processedEventsMix = 0, thread_count = std::thread::hardware_concurrency();

        // ----- MAIN EVENT LOOP -----
        auto start_full = std::chrono::high_resolution_clock::now();
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            if (i % 500 == 0) printProgressBar(i + 1, nentries);

            hiBinProxy = (double)hiBin; NtrkProxy = (double)Ntrk; HFsumETProxy = (double)HFsumET;

            // ----- Processing Limit Check -----
            bool sigLimitReached = cfg.useTestLimitSig && (processedEventsSig >= cfg.testLimitSig);
            bool mixLimitReached = cfg.useTestLimitMix && (processedEventsMix >= cfg.testLimitMix);
            if (sigLimitReached || mixLimitReached) break;

            // ----- Event-level filters -----
            if (!(*selectionVar >= cfg.selVarMoreeq && *selectionVar < cfg.selVarLess)) continue;
            if (std::abs(pvZ) >= cfg.vertexDistance) continue;
            
            // --- Pile-Up Cut ---
            if (!IsGoodEventPileUp(zdcSum, HFsumET)) continue;

            // ----- Z-Bin calculation -----
            int zBin = (int)((pvZ + cfg.vertexDistance) / cfg.zBinWidth);
            if (zBin < 0) zBin = 0; if (zBin >= nZBins) zBin = nZBins - 1;

            EventData currentEv;
            for (int j = 0; j < Ntrk; j++) {
                if (!PassesTrackCuts(trkPt[j], trkEta[j], trkNpixLayers[j], trkNvalidHits[j], trkWeight[j])) continue;

                currentEv.tracks.push_back(ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass));
                currentEv.weights.push_back(trkWeight[j]);
                currentEv.charges.push_back(trkCharge[j]);
            }

            if (currentEv.tracks.size() < 2) continue;
            processedEventsSig++;

            // ----- Mixing processing -----
            auto start_mix_lap = std::chrono::high_resolution_clock::now();
            if (!pool[zBin].empty()) {
                processMixQLCMS(thread_count, currentEv, pool[zBin], hMixSS);
                processedEventsMix++;
            }
            auto end_mix_lap = std::chrono::high_resolution_clock::now();
            duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();

            // ----- Signal processing -----
            auto start_signal_lap = std::chrono::high_resolution_clock::now();
            processSignalQLCMS(thread_count, currentEv, hSigSS);
            auto end_signal_lap = std::chrono::high_resolution_clock::now();
            duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();

            // ----- Pool Update -----
            pool[zBin].push_back(currentEv);
            if (pool[zBin].size() > poolSize) pool[zBin].erase(pool[zBin].begin());
        }
        auto end_full = std::chrono::high_resolution_clock::now();
        duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

        // ----- SAVING -----
        char prefix[256];
        sprintf(prefix, "qlcms_eta-abs%g_%s_%g-%g", cfg.etaCut, selectionVarName, displayMoreeq, displayLess);
        TH1D *hists[] = {hSigSS, hMixSS};
        save_histograms(hists, 2, "./data/signal_mix/", prefix, displayMoreeq, displayLess);
        AnalysisLog::instance().save("./logs", "sig_qlcms");
    }

    // ----- qtqz and q0q Analysis (2D) -----
    void buildQtqzq0q() 
    {
        // ----- INITIALIZATION & DATA LOADING -----
        TFile *fr = nullptr; TTree *t = nullptr;
        getFileTree(inputs.fileInput.c_str(), inputs.treeInput.c_str(), fr, t);

        Int_t maxSize = 50000, Ntrk, hiBin, trkCharge[maxSize]; 
        Float_t pvZ, trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize],
                trkNpixLayers[maxSize], trkNvalidHits[maxSize], pionMass = 0.13957039;
        Float_t HFsumET, zdcSum; 

        t->SetBranchAddress("HFsumET", &HFsumET);
        t->SetBranchAddress("zdcSum", &zdcSum);
        t->SetBranchAddress("Ntrk", &Ntrk);
        t->SetBranchAddress("hiBin", &hiBin);
        t->SetBranchAddress("pvZ", &pvZ); 
        t->SetBranchAddress("trkNpixLayers", trkNpixLayers);
        t->SetBranchAddress("trkCharge", trkCharge);
        t->SetBranchAddress("trkWeight", trkWeight);
        t->SetBranchAddress("trkPt", trkPt);
        t->SetBranchAddress("trkEta", trkEta);
        t->SetBranchAddress("trkPhi", trkPhi);
        t->SetBranchAddress("trkNvalidHits", trkNvalidHits);

        // ----- SELECTION VARIABLE SETUP -----
        double hiBinProxy, NtrkProxy, HFsumETProxy;
        double* selectionVar;
        double displayMoreeq = cfg.selVarMoreeq, displayLess = cfg.selVarLess;
        const char* selectionVarName;

        if (cfg.selectionVarType == ControlVar::CENT) {
            selectionVar = &hiBinProxy; cfg.selVarMoreeq *= 2; cfg.selVarLess *= 2; selectionVarName = "CENT";
        } else if (cfg.selectionVarType == ControlVar::MULT) {
            selectionVar = &NtrkProxy; selectionVarName = "MULT";
        } else {
            selectionVar = &HFsumETProxy; selectionVarName = "CENTHF";
        }

        // ----- BINS & POOL SETUP -----
        int nZBins = static_cast<int>(std::ceil((cfg.vertexDistance * 2.0) / cfg.zBinWidth));
        std::vector<std::vector<EventData>> pool(nZBins);
        size_t poolSize = static_cast<size_t>(cfg.poolSizeInt);

        // ----- HISTOGRAMS & TIMING SETUP -----
        double duration_full = 0.0, duration_mix = 0.0, duration_signal = 0.0;
        // Paper-driven range for q^2 maps and coarser binning for stable 2D rendering.
        int nnscale = 150 / 1;
        double x0 = 0., xt = 0.015;

        TH2D* hSigSS_q0q = new TH2D("hSigSS_q0q", ";q_{0}^{2} [GeV^{2}];|#vec{q}|^{2} [GeV^{2}];#Pairs", nnscale, x0, xt, nnscale, x0, xt);
        TH2D* hSigSS_qtqz = new TH2D("hSigSS_qtqz", ";q_{t}^{2} [GeV^{2}];|q_{z}^{2} [GeV^{2}];#Pairs", nnscale, x0, xt, nnscale, x0, xt);
        TH2D* hMixSS_qtqz = new TH2D("hMixSS_qtqz", ";q_{t}^{2} [GeV^{2}];|q_{z}^{2} [GeV^{2}];#Pairs", nnscale, x0, xt, nnscale, x0, xt);
        TH2D* hMixSS_q0q = new TH2D("hMixSS_q0q", ";q_{0}^{2} [GeV^{2}];|#vec{q}|^{2} [GeV^{2}];#Pairs", nnscale, x0, xt, nnscale, x0, xt);

        hSigSS_q0q->Sumw2();
        hSigSS_qtqz->Sumw2();
        hMixSS_qtqz->Sumw2();
        hMixSS_q0q->Sumw2();

        Long64_t nentries = t->GetEntries();
        int processedEventsSig = 0, processedEventsMix = 0, thread_count = std::thread::hardware_concurrency();

        // ----- MAIN EVENT LOOP -----
        auto start_full = std::chrono::high_resolution_clock::now();
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            if (i % 500 == 0) printProgressBar(i + 1, nentries);

            hiBinProxy = (double)hiBin; NtrkProxy = (double)Ntrk; HFsumETProxy = (double)HFsumET;

            bool sigLimitReached = cfg.useTestLimitSig && (processedEventsSig >= cfg.testLimitSig);
            bool mixLimitReached = cfg.useTestLimitMix && (processedEventsMix >= cfg.testLimitMix);
            if (sigLimitReached || mixLimitReached) break;

            if (!(*selectionVar >= cfg.selVarMoreeq && *selectionVar < cfg.selVarLess)) continue;
            if (std::abs(pvZ) >= cfg.vertexDistance) continue;
            
            // --- Pile-Up Cut ---
            if (!IsGoodEventPileUp(zdcSum, HFsumET)) continue;

            // ----- Z-Bin calculation -----
            int zBin = (int)((pvZ + cfg.vertexDistance) / cfg.zBinWidth);
            if (zBin < 0) zBin = 0; if (zBin >= nZBins) zBin = nZBins - 1;

            EventData currentEv;
            for (int j = 0; j < Ntrk; j++) {
                // ----- TRACK CUTS LOGIC -----
                if (!PassesTrackCuts(trkPt[j], trkEta[j], trkNpixLayers[j], trkNvalidHits[j], trkWeight[j])) continue;

                currentEv.tracks.push_back(ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(trkPt[j], trkEta[j], trkPhi[j], pionMass));
                currentEv.weights.push_back(trkWeight[j]);
                currentEv.charges.push_back(trkCharge[j]);
            }

            if (currentEv.tracks.size() < 2) continue;
            processedEventsSig++;

            // ----- Mixing processing -----
            auto start_mix_lap = std::chrono::high_resolution_clock::now();
            if (!pool[zBin].empty()) {
                processMixQtQzQ0Q(thread_count, currentEv, pool[zBin], hMixSS_qtqz, hMixSS_q0q);
                processedEventsMix++;
            }
            auto end_mix_lap = std::chrono::high_resolution_clock::now();
            duration_mix += std::chrono::duration_cast<std::chrono::duration<double>>(end_mix_lap - start_mix_lap).count();

            auto start_signal_lap = std::chrono::high_resolution_clock::now();
            processSignalQtQzQ0Q(thread_count, currentEv, hSigSS_qtqz, hSigSS_q0q);
            auto end_signal_lap = std::chrono::high_resolution_clock::now();
            duration_signal += std::chrono::duration_cast<std::chrono::duration<double>>(end_signal_lap - start_signal_lap).count();

            // ----- Pool Update -----
            pool[zBin].push_back(currentEv);
            if (pool[zBin].size() > poolSize) pool[zBin].erase(pool[zBin].begin());
        }
        auto end_full = std::chrono::high_resolution_clock::now();
        duration_full = std::chrono::duration_cast<std::chrono::duration<double>>(end_full - start_full).count();

        // ----- SAVING -----
        char prefix[256];
        sprintf(prefix, "qtqzq0q_eta-abs%g_%s_%g-%g", cfg.etaCut, selectionVarName, displayMoreeq, displayLess);
        TH2D *hists[] = {hSigSS_q0q, hMixSS_q0q, hSigSS_qtqz, hMixSS_qtqz};
        save_benchmark_chrono({duration_full, duration_signal, duration_mix}, {"Total", "Signal", "Mix"}, "benchmarks", prefix, processedEventsSig, processedEventsMix);
        save_histograms2d(hists, 4, "./data/signal_mix/", prefix, displayMoreeq, displayLess);
        AnalysisLog::instance().save("./logs", "sig_qtqzq0q");
        if (fr) { fr->Close(); delete fr; }
    }
};

// ----- Simplified Testing Main -----
int main() {
// to compile use:
// g++ -std=c++17 -pthread makeSignal.cpp -o makeSignal `root-config --cflags --libs` -lGenVector
// to run use:
// ./makeSignal
    ROOT::EnableImplicitMT();

    // ----- Configuration -----
    InputParams inputs;
    inputs.fileInput = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";
    inputs.treeInput = "demo/TreeMBUCC";

    std::vector<double> centralityBins = {3200, 3300}; 
    
    // Base parameters setup
    AnalysisParams baseCfg;
    baseCfg.selectionVarType = ControlVar::CENTHF;
    baseCfg.ptMin = 0.5;
    baseCfg.usePtMax = false; 
    baseCfg.etaCut = 0.95;
    baseCfg.pixHit = 1.0;
    baseCfg.nvalidHits = 6.0;
    baseCfg.vertexDistance = 15.0;
    baseCfg.poolSizeInt = 10;
    baseCfg.zBinWidth = 2.0;
    
    baseCfg.useTestLimitSig = false;
    baseCfg.useTestLimitMix = false;
    baseCfg.usePixHit = true;
    baseCfg.useNvalidHits = false;
    
    // --- New Pile-Up Cut Toggles ---
    baseCfg.usePileUpCut = true; // Activated by default based on requirement
    
    // --- New Split Cut Toggles ---
    baseCfg.useSplitCut = true;  // Activated by default based on requirement

    // ----- Processing Loop -----
    std::cout << "Starting analysis for " << centralityBins.size() - 1 << " bin(s)." << std::endl;
    
    for (size_t i = 0; i < centralityBins.size() - 1; ++i) {
        AnalysisParams currentCfg = baseCfg; // Copy base settings
        
        currentCfg.selVarMoreeq = centralityBins[i];
        currentCfg.selVarLess = centralityBins[i+1];

        std::cout << "\n--- Running for bin: " << currentCfg.selVarMoreeq << " - " << currentCfg.selVarLess << " ---" << std::endl;
        
        Process analysis(inputs, currentCfg);
        
        // Running all analysis methods
        //analysis.buildDeltaEtaDeltaPhi();
        //analysis.buildQinv();
        analysis.buildQlcms();
        //analysis.buildQtqzq0q();
    }

    std::cout << "\nAll analysis bins finished." << std::endl;
    return 0;
}
