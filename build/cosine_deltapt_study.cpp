#include <iostream>
#include <vector>
#include <cmath>
#include <string>
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

// --- Configuração do Mixing e Vertex ---
const float VERTEX_Z_MAX = 15.0; 

// --- Cortes Básicos de Análise ---
const double PT_MIN = 0.3;       // Alterado para 0.3 GeV
const double ETA_MIN_CUT = 0.0;  
const double ETA_MAX_CUT = 0.95;  
const double PIX_HIT_MIN = 1.0;
const double HF_MIN = 3200.0; 
const double HF_MAX = 3300.0; 

// Cortes Históricos para Referência Visual (Linhas no gráfico)
const double HIST_COS_CUT = 0.99996;
const double HIST_DPT_CUT = 0.04;

// Estrutura auxiliar para cálculo vetorial rápido
struct VectorInfo {
    double px, py, pz, p, pt, eta, phi;
    int charge;
    float weight;
    
    VectorInfo(float pt_in, float eta_in, float phi_in, int q_in, float w_in) {
        pt = pt_in; eta = eta_in; phi = phi_in;
        charge = q_in; weight = w_in;
        
        px = pt * std::cos(phi);
        py = pt * std::sin(phi);
        pz = pt * std::sinh(eta);
        p  = pt * std::cosh(eta);
    }
};

void draw_info_legend() {
    TLegend* leg = new TLegend(0.40, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);    
    leg->SetHeader("#bf{Selections (No #Delta#eta/#Delta#phi cut)}");
    leg->AddEntry((TObject*)0, Form("Event: %.0f < HF < %.0f", HF_MIN, HF_MAX), "");
    // Atualizado pT min para 0.3
    leg->AddEntry((TObject*)0, Form("Track: p_{T} > %.1f, %.2f < |#eta| < %.2f", PT_MIN, ETA_MIN_CUT, ETA_MAX_CUT), "");
    leg->AddEntry((TObject*)0, "Pair: Same-Sign Only", "");
    leg->Draw();
}

void analyze_vectors(const char* inputFile, const char* treeName) {
    
    TFile* f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) { std::cerr << "Erro no arquivo" << std::endl; return; }
    TTree* t = (TTree*)f->Get(treeName);
    if (!t) { std::cerr << "Erro na Tree" << std::endl; return; }

    const int MAX_TRK = 50000;
    Int_t Ntrk;
    Int_t trkCharge[MAX_TRK]; 
    Float_t trkPt[MAX_TRK], trkEta[MAX_TRK], trkPhi[MAX_TRK], trkWeight[MAX_TRK];
    Float_t trkNpixLayers[MAX_TRK]; 
    Float_t HFsumET, pvZ;

    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("HFsumET", &HFsumET); 
    t->SetBranchAddress("pvZ", &pvZ);         
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    t->SetBranchAddress("trkCharge", trkCharge); 
    t->SetBranchAddress("trkWeight", trkWeight); 
    
    if(t->GetBranch("trkNpixLayers")) t->SetBranchAddress("trkNpixLayers", trkNpixLayers);

    // --- Histogramas ---
    
    // 1. Zoom Extremo (Direita) - Região do Split
    TH2D* h_cosa_deltapt = new TH2D("h_cosa_deltapt", 
        "Splitting Region (Ultra Zoom);cos(#theta);#Delta p_{T} [GeV];#Pairs", 
        200, 0.9995, 1.000005, // X-Axis: Ultra Zoom
        100, 0, 0.2);          // Y-Axis: Zoom

    // 2. Zoom Médio (Esquerda) - Região Solicitada
    // Ranges: 0.9 < cos(theta) < 1.0 e 0 < dPt < 0.1
    TH2D* h_cosa_deltapt_med = new TH2D("h_cosa_deltapt_med", 
        "Medium Zoom;cos(#theta);#Delta p_{T} [GeV];#Pairs", 
        200, 0.9, 1.0,  // X-Axis ajustado
        200, 0, 0.1);   // Y-Axis ajustado

    h_cosa_deltapt->Sumw2();
    h_cosa_deltapt_med->Sumw2();

    Long64_t nentries = t->GetEntries();
    std::cout << "Processando " << nentries << " eventos..." << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        if (i % 1000 == 0) {
            float prog = (float)(i)/(float)nentries;
            std::cout << int(prog*100) << " %\r" << std::flush;
        }

        if (!(HFsumET >= HF_MIN && HFsumET < HF_MAX)) continue;
        if (std::abs(pvZ) >= VERTEX_Z_MAX) continue;

        std::vector<VectorInfo> good_tracks;
        good_tracks.reserve(Ntrk);

        for (int j = 0; j < Ntrk; j++) {
            if (trkPt[j] < PT_MIN) continue; // PT > 0.3
            
            double absEta = std::abs(trkEta[j]);
            if (absEta <= ETA_MIN_CUT || absEta >= ETA_MAX_CUT) continue;
            if (trkNpixLayers[j] < PIX_HIT_MIN) continue;
            if (std::isinf(trkWeight[j])) continue;

            good_tracks.emplace_back(trkPt[j], trkEta[j], trkPhi[j], trkCharge[j], trkWeight[j]);
        }

        for (size_t it = 0; it < good_tracks.size(); ++it) {
            for (size_t jt = it + 1; jt < good_tracks.size(); ++jt) {
                const auto& vec1 = good_tracks[it];
                const auto& vec2 = good_tracks[jt];

                // Filtro Same-Sign Apenas
                if (vec1.charge * vec2.charge <= 0) continue;

                double dot_product = vec1.px * vec2.px + vec1.py * vec2.py + vec1.pz * vec2.pz;
                double mag_prod = vec1.p * vec2.p;
                
                Double_t cosa = TMath::Abs(dot_product) / mag_prod;
                Double_t deltapt = TMath::Abs(vec1.pt - vec2.pt);

                h_cosa_deltapt->Fill(cosa, deltapt, vec1.weight * vec2.weight);
                h_cosa_deltapt_med->Fill(cosa, deltapt, vec1.weight * vec2.weight);
            }
        }
    }
    std::cout << "100% Done." << std::endl;

    // --- PLOTTING ---
    TCanvas* c = new TCanvas("c_vectors", "Vector Analysis", 1200, 600);
    c->Divide(2, 1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    // Pad 1: Medium Zoom (0.9 < cos < 1.0, dPt < 0.1)
    c->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetTheta(30); 
    gPad->SetPhi(30);
    
    h_cosa_deltapt_med->GetXaxis()->SetTitleOffset(1.8);
    h_cosa_deltapt_med->GetYaxis()->SetTitleOffset(1.8);
    h_cosa_deltapt_med->GetZaxis()->SetTitleOffset(1.2);
    h_cosa_deltapt_med->Draw("LEGO2 0");
    draw_info_legend();
    
    // Pad 2: Ultra Zoomed Split Region (0.9995 < cos < 1.0)
    c->cd(2);
    gPad->SetRightMargin(0.10);
    gPad->SetLeftMargin(0.15);
    gPad->SetTheta(30); 
    gPad->SetPhi(30);
    
    h_cosa_deltapt->GetXaxis()->SetTitleOffset(1.8);
    h_cosa_deltapt->GetYaxis()->SetTitleOffset(1.8);
    h_cosa_deltapt->GetZaxis()->SetTitleOffset(1.2);
    h_cosa_deltapt->Draw("LEGO2 0"); 
    
    TLatex lat; 
    lat.SetNDC(); lat.SetTextSize(0.035);
    lat.DrawLatex(0.1, 0.92, "#bf{Ultra Zoom: Split Region}");
    lat.DrawLatex(0.1, 0.87, Form("Hist Cuts: cos(#theta) > %.5f, #Delta p_{T} < %.2f", HIST_COS_CUT, HIST_DPT_CUT));

    c->SaveAs("imgs/test/control_plots/investigation_CosTheta_DeltaPt_LEGO.pdf");
    c->SaveAs("imgs/test/control_plots/investigation_CosTheta_DeltaPt_LEGO.png");

    delete h_cosa_deltapt;
    delete h_cosa_deltapt_med;
}

int main() {
    const char* file = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";
    const char* tree = "demo/TreeMBUCC";
    analyze_vectors(file, tree);
    return 0;
}