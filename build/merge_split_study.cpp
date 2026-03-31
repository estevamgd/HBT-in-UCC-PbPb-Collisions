#include <iostream>
#include <vector>
#include <deque>
#include <cmath>
#include <string>
#include <thread>
#include <vector>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector2.h"
#include "TLatex.h"
#include "TLine.h" 
#include "TColor.h" 

// --- Configuração da Região de Interesse (ROI) ---
const double ROI_DELTA_ETA_MAX = 0.01; 
const double ROI_DELTA_PHI_MIN = 0.00; 
const double ROI_DELTA_PHI_MAX = 0.10; 

// --- Configuração do Mixing e Vertex ---
const int MIXING_POOL_SIZE = 10; 
const float VERTEX_Z_MAX = 15.0; 
const float Z_BIN_WIDTH = 2.0;   

// --- Cortes Básicos de Análise ---
const double PT_MIN = 0.3;       // Changed to 0.3
const double PT_MAX = 1.0;       // Added upper limit 1.0
const double ETA_MIN_CUT = 0.0;  // Lower limit for |eta|
const double ETA_MAX_CUT = 0.95;  // Upper limit for |eta|
const double PIX_HIT_MIN = 1.0;
const double HF_MIN = 3200.0; 
const double HF_MAX = 3300.0; 

// Estrutura leve para guardar tracks no buffer de mixing
struct MiniTrack {
    float pt, eta, phi, weight;
    int charge;
    float npix, ptRes, dzSig, dxySig; 
};

// Estrutura para agrupar histogramas de sinal (facilita clonagem em threads)
struct SignalHistos {
    TH2D* h_sig;
    TH1D *h_pt, *h_eta, *h_phi, *h_npix, *h_ptRes, *h_dzSig, *h_dxySig;

    SignalHistos(const SignalHistos* base) {
        h_sig   = (TH2D*)base->h_sig->Clone(); h_sig->Reset();
        h_pt    = (TH1D*)base->h_pt->Clone(); h_pt->Reset();
        h_eta   = (TH1D*)base->h_eta->Clone(); h_eta->Reset();
        h_phi   = (TH1D*)base->h_phi->Clone(); h_phi->Reset();
        h_npix  = (TH1D*)base->h_npix->Clone(); h_npix->Reset();
        h_ptRes = (TH1D*)base->h_ptRes->Clone(); h_ptRes->Reset();
        h_dzSig = (TH1D*)base->h_dzSig->Clone(); h_dzSig->Reset();
        h_dxySig= (TH1D*)base->h_dxySig->Clone(); h_dxySig->Reset();
    }
    
    SignalHistos(TH2D* s, TH1D* p, TH1D* e, TH1D* ph, TH1D* n, TH1D* pr, TH1D* dz, TH1D* dxy) 
        : h_sig(s), h_pt(p), h_eta(e), h_phi(ph), h_npix(n), h_ptRes(pr), h_dzSig(dz), h_dxySig(dxy) {}

    void Delete() {
        delete h_sig; delete h_pt; delete h_eta; delete h_phi;
        delete h_npix; delete h_ptRes; delete h_dzSig; delete h_dxySig;
    }
};

// Função para desenhar legenda
void draw_info_legend() {
    TLegend* leg = new TLegend(0.40, 0.60, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->SetHeader("#bf{Selections & ROI}");
    leg->AddEntry((TObject*)0, Form("Event: %.0f < HF < %.0f", HF_MIN, HF_MAX), "");
    leg->AddEntry((TObject*)0, Form("Event: |v_{z}| < %.1f cm", VERTEX_Z_MAX), "");
    // Updated Legend for pT and Eta ranges
    leg->AddEntry((TObject*)0, Form("Track: %.1f < p_{T} < %.1f, %.1f < |#eta| < %.1f", PT_MIN, PT_MAX, ETA_MIN_CUT, ETA_MAX_CUT), "");
    leg->AddEntry((TObject*)0, "Pair: Same-Sign Only", "");
    leg->AddEntry((TObject*)0, Form("Correlation ROI:"), "");
    leg->AddEntry((TObject*)0, Form("#Delta#eta < %.3f", ROI_DELTA_ETA_MAX), "");
    leg->AddEntry((TObject*)0, Form("%.2f < #Delta#phi < %.2f", ROI_DELTA_PHI_MIN, ROI_DELTA_PHI_MAX), "");
    leg->Draw();
}

// --- Funções de Processamento Paralelo ---

void processSignalParallel(const std::vector<MiniTrack>& tracks, SignalHistos* mainHistos, int thread_count) {
    size_t n_tracks = tracks.size();
    if (n_tracks <= 1) return;
    if (thread_count < 1) thread_count = 1;

    std::vector<std::thread> threads;
    std::vector<SignalHistos*> localHistos(thread_count);

    for(int i=0; i<thread_count; ++i) {
        localHistos[i] = new SignalHistos(mainHistos);
    }

    auto task = [&](size_t start_idx, size_t end_idx, SignalHistos* h) {
        for (size_t it = start_idx; it < end_idx; ++it) {
            for (size_t jt = it + 1; jt < n_tracks; ++jt) {
                const auto& t1 = tracks[it];
                const auto& t2 = tracks[jt];

                if (t1.charge * t2.charge <= 0) continue; // Same-Sign Check

                double dEta = std::abs(t1.eta - t2.eta);
                double dPhi = std::abs(TVector2::Phi_mpi_pi(t1.phi - t2.phi));

                if (dEta < ROI_DELTA_ETA_MAX && dPhi >= ROI_DELTA_PHI_MIN && dPhi < ROI_DELTA_PHI_MAX) {
                    h->h_sig->Fill(dEta, dPhi, t1.weight * t2.weight);

                    h->h_pt->Fill(t1.pt); h->h_pt->Fill(t2.pt);
                    h->h_eta->Fill(t1.eta); h->h_eta->Fill(t2.eta);
                    h->h_phi->Fill(t1.phi); h->h_phi->Fill(t2.phi);
                    h->h_npix->Fill(t1.npix); h->h_npix->Fill(t2.npix);
                    h->h_ptRes->Fill(t1.ptRes); h->h_ptRes->Fill(t2.ptRes);
                    h->h_dzSig->Fill(t1.dzSig); h->h_dzSig->Fill(t2.dzSig);
                    h->h_dxySig->Fill(t1.dxySig); h->h_dxySig->Fill(t2.dxySig);
                }
            }
        }
    };

    size_t chunk = n_tracks / thread_count;
    for(int t=0; t<thread_count; ++t) {
        size_t start = t * chunk;
        size_t end = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk;
        threads.emplace_back(task, start, end, localHistos[t]);
    }

    for(auto& t : threads) t.join();
    for(int i=0; i<thread_count; ++i) {
        mainHistos->h_sig->Add(localHistos[i]->h_sig);
        mainHistos->h_pt->Add(localHistos[i]->h_pt);
        mainHistos->h_eta->Add(localHistos[i]->h_eta);
        mainHistos->h_phi->Add(localHistos[i]->h_phi);
        mainHistos->h_npix->Add(localHistos[i]->h_npix);
        mainHistos->h_ptRes->Add(localHistos[i]->h_ptRes);
        mainHistos->h_dzSig->Add(localHistos[i]->h_dzSig);
        mainHistos->h_dxySig->Add(localHistos[i]->h_dxySig);
        localHistos[i]->Delete(); delete localHistos[i];
    }
}

void processMixParallel(const std::vector<MiniTrack>& current_tracks, const std::vector<std::vector<MiniTrack>>& pool, TH2D* h_mix, int thread_count) {
    if (current_tracks.empty() || pool.empty()) return;
    if (thread_count < 1) thread_count = 1;

    std::vector<std::thread> threads;
    std::vector<TH2D*> localHistos(thread_count);

    for(int i=0; i<thread_count; ++i) {
        localHistos[i] = (TH2D*)h_mix->Clone();
        localHistos[i]->Reset();
    }

    auto task = [&](size_t start_idx, size_t end_idx, TH2D* h) {
        for (size_t it = start_idx; it < end_idx; ++it) {
            const auto& t1 = current_tracks[it];
            for (const auto& prev_event_tracks : pool) {
                for (const auto& t2 : prev_event_tracks) {
                    if (t1.charge * t2.charge <= 0) continue; // Same-Sign Check

                    double dEta = std::abs(t1.eta - t2.eta);
                    double dPhi = std::abs(TVector2::Phi_mpi_pi(t1.phi - t2.phi));

                    if (dEta < ROI_DELTA_ETA_MAX && dPhi >= ROI_DELTA_PHI_MIN && dPhi < ROI_DELTA_PHI_MAX) {
                        h->Fill(dEta, dPhi, t1.weight * t2.weight);
                    }
                }
            }
        }
    };

    size_t n_tracks = current_tracks.size();
    size_t chunk = n_tracks / thread_count;
    for(int t=0; t<thread_count; ++t) {
        size_t start = t * chunk;
        size_t end = (t == thread_count - 1) ? n_tracks : (t + 1) * chunk;
        threads.emplace_back(task, start, end, localHistos[t]);
    }

    for(auto& t : threads) t.join();
    for(int i=0; i<thread_count; ++i) {
        h_mix->Add(localHistos[i]);
        delete localHistos[i];
    }
}


void study_correlation_region(const char* inputFile, const char* treeName) {
    
    TFile* f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) { std::cerr << "Erro no arquivo" << std::endl; return; }
    TTree* t = (TTree*)f->Get(treeName);
    if (!t) { std::cerr << "Erro na Tree" << std::endl; return; }

    const int MAX_TRK = 50000;
    Int_t Ntrk;
    Int_t trkCharge[MAX_TRK]; 
    Float_t trkPt[MAX_TRK], trkEta[MAX_TRK], trkPhi[MAX_TRK], trkWeight[MAX_TRK];
    Float_t trkPtRes[MAX_TRK], trkDzSig[MAX_TRK], trkDxySig[MAX_TRK], trkNpixLayers[MAX_TRK]; 
    Float_t HFsumET, pvZ;

    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("HFsumET", &HFsumET); 
    t->SetBranchAddress("pvZ", &pvZ);         
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    t->SetBranchAddress("trkCharge", trkCharge); 
    t->SetBranchAddress("trkWeight", trkWeight); 
    
    if(t->GetBranch("trkPtRes")) t->SetBranchAddress("trkPtRes", trkPtRes);
    if(t->GetBranch("trkDzSig")) t->SetBranchAddress("trkDzSig", trkDzSig);
    if(t->GetBranch("trkDxySig")) t->SetBranchAddress("trkDxySig", trkDxySig);
    if(t->GetBranch("trkNpixLayers")) t->SetBranchAddress("trkNpixLayers", trkNpixLayers);

    // --- Histogramas ---
    int nBinsEta = 50;
    int nBinsPhi = 50;
    TH2D* h_sig = new TH2D("h_sig", "Signal", nBinsEta, 0, ROI_DELTA_ETA_MAX, nBinsPhi, ROI_DELTA_PHI_MIN, ROI_DELTA_PHI_MAX);
    TH2D* h_mix = new TH2D("h_mix", "Mixed Event", nBinsEta, 0, ROI_DELTA_ETA_MAX, nBinsPhi, ROI_DELTA_PHI_MIN, ROI_DELTA_PHI_MAX);
    
    h_sig->Sumw2();
    h_mix->Sumw2();

    TH1D* h_pt      = new TH1D("h_pt", "p_{T} of ROI Tracks;p_{T} [GeV];#Events", 100, 0, 2.0); // Range 0-2 for 0.3-1.0 cut
    TH1D* h_eta     = new TH1D("h_eta", "#eta of ROI Tracks;#eta;#Events", 100, -3, 3);
    TH1D* h_phi     = new TH1D("h_phi", "#phi of ROI Tracks;#phi;#Events", 100, -3.2, 3.2);
    TH1D* h_npix    = new TH1D("h_npix", "Pixel Layers;N Pixel Layers;#Events", 10, -0.5, 9.5);
    TH1D* h_ptRes   = new TH1D("h_ptRes", "p_{T} Resolution;#sigma_{pT}/p_{T};#Events", 100, 0, 0.2);
    TH1D* h_dzSig   = new TH1D("h_dzSig", "z-Vertex Sig;dzSig;#Events", 100, -10, 10);
    TH1D* h_dxySig  = new TH1D("h_dxySig", "xy-Vertex Sig;dxySig;#Events", 100, -10, 10);

    SignalHistos mainSigHistos(h_sig, h_pt, h_eta, h_phi, h_npix, h_ptRes, h_dzSig, h_dxySig);

    int nZBins = static_cast<int>(std::ceil((VERTEX_Z_MAX * 2.0) / Z_BIN_WIDTH));
    std::vector<std::deque<std::vector<MiniTrack>>> mixing_pools(nZBins);

    int n_threads = std::thread::hardware_concurrency();
    
    Long64_t nentries = t->GetEntries();
    std::cout << "Processando " << nentries << " eventos com " << n_threads << " threads e " << nZBins << " Z-bins..." << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        if (i % 1000 == 0) {
            float prog = (float)(i)/(float)nentries;
            std::cout << int(prog*100) << " %\r" << std::flush;
        }

        if (!(HFsumET >= HF_MIN && HFsumET < HF_MAX)) continue;
        if (std::abs(pvZ) >= VERTEX_Z_MAX) continue;

        int zBin = (int)((pvZ + VERTEX_Z_MAX) / Z_BIN_WIDTH);
        if (zBin < 0) zBin = 0; 
        if (zBin >= nZBins) zBin = nZBins - 1;

        std::vector<MiniTrack> current_tracks;
        for (int j = 0; j < Ntrk; j++) {
            // --- Cuts Updated: 0.3 < pT < 1.0, 0.1 < |eta| < 1.0 ---
            if (trkPt[j] < PT_MIN || trkPt[j] > PT_MAX) continue;
            
            double absEta = std::abs(trkEta[j]);
            if (absEta <= ETA_MIN_CUT || absEta >= ETA_MAX_CUT) continue;

            if (trkNpixLayers[j] < PIX_HIT_MIN) continue;
            if (std::isinf(trkWeight[j])) continue;

            MiniTrack m;
            m.pt = trkPt[j]; m.eta = trkEta[j]; m.phi = trkPhi[j];
            m.charge = trkCharge[j]; m.weight = trkWeight[j];
            m.npix = trkNpixLayers[j]; m.ptRes = trkPtRes[j]; 
            m.dzSig = trkDzSig[j]; m.dxySig = trkDxySig[j];
            
            current_tracks.push_back(m);
        }

        processSignalParallel(current_tracks, &mainSigHistos, n_threads);

        if (!mixing_pools[zBin].empty()) {
            std::vector<std::vector<MiniTrack>> pool_vec(mixing_pools[zBin].begin(), mixing_pools[zBin].end());
            processMixParallel(current_tracks, pool_vec, h_mix, n_threads);
        }

        mixing_pools[zBin].push_back(current_tracks);
        if (mixing_pools[zBin].size() > MIXING_POOL_SIZE) {
            mixing_pools[zBin].pop_front();
        }
    }
    std::cout << "100% Done." << std::endl;

    TH2D* h_corr = (TH2D*)h_sig->Clone("h_corr");
    h_corr->SetTitle("Correlation C(#Delta#eta, #Delta#varphi) in ROI;|#Delta#eta|;|#Delta#phi|;C");
    
    double normNum = h_sig->Integral();
    double normDen = h_mix->Integral();
    
    h_corr->Divide(h_mix); 
    if (normNum > 0 && normDen > 0) {
        h_corr->Scale(normDen / normNum); 
    }

    TCanvas* c = new TCanvas("c_corr_study", "Correlation Study in Split Region", 1200, 800);
    c->Divide(4, 2);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow); 

    for(int i=1; i<=8; ++i) {
        c->cd(i);
        gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.13); gPad->SetBottomMargin(0.15);
    }

    auto style_h = [](TH1* h) {
        h->GetYaxis()->SetTitleOffset(1.6); h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
    };
    style_h(h_corr); style_h(h_pt); style_h(h_eta); style_h(h_phi);
    style_h(h_npix); style_h(h_ptRes); style_h(h_dzSig); style_h(h_dxySig);

    c->cd(1); 
    h_corr->GetZaxis()->SetRangeUser(0.6, 1.4); 
    h_corr->Draw("COLZ"); draw_info_legend();
    c->cd(2); h_pt->Draw(); gPad->SetLogy();
    c->cd(3); h_eta->Draw(); gPad->SetLogy();
    c->cd(4); h_phi->Draw();
    
    c->cd(5); h_npix->Draw(); 
    c->cd(6); h_ptRes->Draw();
    c->cd(7); h_dzSig->Draw();
    c->cd(8); h_dxySig->Draw(); gPad->SetLogy();

    c->SaveAs("imgs/test/control_plots/investigation_Correlation_ROI.pdf");
    c->SaveAs("imgs/test/control_plots/investigation_Correlation_ROI.png");

    delete h_sig; delete h_mix;
}

int main() {
    const char* file = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";
    const char* tree = "demo/TreeMBUCC";
    study_correlation_region(file, tree);
    return 0;
}