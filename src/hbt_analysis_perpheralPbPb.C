#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TH1D.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <TStyle.h>
#include "TLegend.h"
#include <TText.h>
#include <TBenchmark.h>
#include <TProfile.h>
#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"

void hbt_analysis_perpheralPbPb() {
    // Keeping track of time
    time_t ttime1 = time(NULL);
    struct tm date1 = *localtime(&ttime1);
    std::cout << "Start: " << date1.tm_mday << "/" << date1.tm_mon + 1 << "/" << date1.tm_year + 1900 <<  " " 
        << date1.tm_hour << ":" << date1.tm_min << ":" << date1.tm_sec << std::endl;

    TFile *fr = nullptr;
    TTree *t = nullptr;

    getFileTree("data/HiForestAOD_HIMinBiasUPC_Skim_rereco_2760PbPbMB_pixeltracks_0000_300kevts_part1.root", "demo/TreeMBUCC", fr, t);

    // Variables
    Int_t maxSize = 17100;
    Int_t Ntrk = maxSize, trkCharge[Ntrk];
    Float_t HFsumET, HFsumETPlus, HFsumETMinus;
    Float_t trkPt[Ntrk], trkEta[Ntrk], trkPhi[Ntrk], trkPtRes[Ntrk], 
            trkDxySig[Ntrk], trkNpixLayers[Ntrk], trkDzSig[Ntrk], trkWeight[Ntrk];

    // Arrays of variables
    void* variables[] = {
        &HFsumET, &HFsumETPlus, &HFsumETMinus, &Ntrk, trkPt, trkEta, trkPhi, trkPtRes, trkDxySig, 
        trkNpixLayers, trkDzSig, trkCharge, trkWeight
    };

    // Branch names
    const char* branchNames[] = {
        "HFsumET", "HFsumETPlus", "HFsumETMinus", "Ntrk", "trkPt", "trkEta", "trkPhi", "trkPtRes", 
        "trkDxySig", "trkNpixLayers", "trkDzSig", "trkCharge", "trkWeight"
    };

    // Setting addresses
    int numBranches = sizeof(variables) / sizeof(variables[0]);
    for (int i = 0; i < numBranches; i++) {
        t->SetBranchAddress(branchNames[i], variables[i]);
    }

    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    double d_nentries = t->GetEntries();

    // Setting canvases
    TCanvas *c1 = new TCanvas("c1", "Histograms", 7680, 4320);
    c1->Divide(3, 2);

    TCanvas *c2 = new TCanvas("c2", "Histograms", 7680, 4320);
    c2->Divide(3, 2);

    TCanvas *c3 = new TCanvas("c3", "Histograms", 7680, 4320);
    c3->Divide(3,2);

    TCanvas *canvases[] = {c1, c2, c3};
    int numCanvases = 3;

    // Setting histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;

    TH1D* hHFsumET = cHist("hHFsumET", "HFsumET[GeV]", "#Events/bin", 0, 4000, ninterval, nlength, nscale);
    TH1D* hHFsumETPlus = cHist("hHFsumETPlus", "HFsumETPlus[GeV]", "#Events/bin", 0, 2500, ninterval, nlength, nscale);
    TH1D* hHFsumETMinus = cHist("hHFsumETMinus", "HFsumETMinus[GeV]", "#Events/bin", 0, 2500, ninterval, nlength, nscale);
    TH1D* hNtrk = cHist("hNtrk", "Ntrk", "#Events/bin", 0, 4000, ninterval, nlength, nscale);
    TH1D* htrkCharge = cHist("htrkCharge", "trkCharge", "#Tracks/bin", -2, 2, ninterval, nlength, nscale);
    TH1D* htrkWeight = cHist("htrkWeight", "trkWeight", "#Tracks/bin", 0, 20, ninterval, nlength, nscale);
    TH1D* htrkPt = cHist("htrkPt", "pT[GeV]", "#Tracks/bin", 0, 45, ninterval, nlength, nscale);
    TH1D* htrkEta = cHist("htrkEta", "trkEta", "#Tracks/bin", -3, 3, ninterval, nlength, nscale);
    TH1D* htrkPhi = cHist("htrkPhi", "trkPhi", "#Tracks/bin", -3.4, 3.4, ninterval, nlength, nscale);
    TH1D* htrkPtRes = cHist("htrkPtRes", "trkPtRes", "#Tracks/bin", 0, 1, ninterval, nlength, nscale);
    TH1D* htrkDxySig = cHist("htrkDxySig", "trkDxySig", "#Tracks/bin", -3, 3, ninterval, nlength, nscale);
    TH1D* htrkDzSig = cHist("htrkDzSig", "trkDzSig", "#Tracks/bin", -3.5, 3.5, ninterval, nlength, nscale);
    TH1D* htrkNpixLayers = cHist("htrkNpixLayers", "trkNpixLayers", "#Tracks/bin", 0, 6, ninterval, nlength, nscale);
    
    TH1D* htrkPtCor = cHist("htrkPtCor", "pT[GeV]", "#Tracks/bin", 0, 45, ninterval, nlength, nscale);
    TH1D* htrkEtaCor = cHist("htrkEtaCor", "trkEta", "#Tracks/bin", -3, 3, ninterval, nlength, nscale);
    TH1D* htrkPhiCor = cHist("htrkPhiCor", "trkPhi", "#Tracks/bin", -3.4, 3.4, ninterval, nlength, nscale);
    
    TH1D* histograms[] = {hHFsumET, hHFsumETPlus, hHFsumETMinus, hNtrk, htrkCharge, htrkWeight, htrkPt, 
                            htrkEta, htrkPhi, htrkPtRes, htrkDxySig, htrkDzSig, htrkNpixLayers, htrkPtCor, htrkEtaCor, htrkPhiCor};
    int numHistograms = 16;

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        hHFsumET->Fill(HFsumET);
        hHFsumETPlus->Fill(HFsumETPlus);
        hHFsumETMinus->Fill(HFsumETMinus);
        hNtrk->Fill(Ntrk);

        for (int j = 0; j < Ntrk; j++) {
            if (!TMath::IsNaN(trkWeight[j]) && TMath::Finite(trkWeight[j])) {
                htrkPtCor->Fill(trkPt[j], trkWeight[j]);
                htrkEtaCor->Fill(trkEta[j], trkWeight[j]);
                htrkPhiCor->Fill(trkPhi[j], trkWeight[j]);
            } 
            htrkCharge->Fill(trkCharge[j]);
            htrkWeight->Fill(trkWeight[j]);
            htrkPt->Fill(trkPt[j]);
            htrkEta->Fill(trkEta[j]);
            htrkPhi->Fill(trkPhi[j]);
            htrkPtRes->Fill(trkPtRes[j]);
            htrkDxySig->Fill(trkDxySig[j]);
            htrkDzSig->Fill(trkDzSig[j]);
            htrkNpixLayers->Fill(trkNpixLayers[j]);
        }
        
        if (i % 100 == 0) {
            float progress = (float)(i + 1) / nentries * 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
        }
        
    }
    std::cout<<std::endl;

    TH1D *htrkPt_norm = (TH1D *)htrkPt->Clone("htrkPt_norm");
    TH1D *htrkPt_norm_cor = (TH1D *)htrkPtCor->Clone("htrkPt_norm_cor");
    TH1D *htrkEta_norm = (TH1D *)htrkEta->Clone("htrkEta_norm");
    TH1D *htrkEta_norm_cor = (TH1D *)htrkEtaCor->Clone("htrkEta_norm_cor");
    
    htrkPt_norm->Scale(1.0/htrkPt_norm->Integral());
    htrkPt_norm_cor->Scale(1.0/htrkPt_norm_cor->Integral());
    htrkEta_norm->Scale(1.0/htrkEta_norm->Integral());
    htrkEta_norm_cor->Scale(1.0/htrkEta_norm_cor->Integral());
    
    // Drawing histograms
    c1->cd(1); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hHFsumET->Draw();
    c1->cd(2); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hHFsumETPlus->Draw();
    c1->cd(3); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hHFsumETMinus->Draw();
    c1->cd(4); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hNtrk->Draw();
    c1->cd(5); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkCharge->Draw();
    c1->cd(6); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkWeight->Draw();
    
    c2->cd(1); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkPt->Draw("HIST");
    htrkPt_norm->SetLineColor(kRed); htrkPt_norm->SetLineWidth(5); htrkPt_norm->Draw("HIST SAME");
    htrkPt_norm_cor->SetLineColor(kBlue); htrkPt_norm_cor->SetLineWidth(5); htrkPt_norm_cor->Draw("HIST SAME");

    c2->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkEta->Draw("HIST");
    htrkEta_norm->SetLineColor(kRed); htrkEta_norm->SetLineWidth(5); htrkEta_norm->Draw("HIST SAME");
    htrkEta_norm_cor->SetLineColor(kBlue); htrkEta_norm_cor->SetLineWidth(5); htrkEta_norm_cor->Draw("HIST SAME");

    c2->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkPhi->Draw();

    c2->cd(4); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkPtCor->Draw("HIST");
    htrkPt_norm->SetLineColor(kRed); htrkPt_norm->SetLineWidth(5); htrkPt_norm->Draw("HIST SAME");
    htrkPt_norm_cor->SetLineColor(kBlue); htrkPt_norm_cor->SetLineWidth(5); htrkPt_norm_cor->Draw("HIST SAME");

    c2->cd(5); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkEtaCor->Draw("HIST");
    htrkEta_norm->SetLineColor(kRed); htrkEta_norm->SetLineWidth(5); htrkEta_norm->Draw("HIST SAME");
    htrkEta_norm_cor->SetLineColor(kBlue); htrkEta_norm_cor->SetLineWidth(5); htrkEta_norm_cor->Draw("HIST SAME");

    c2->cd(6); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkPhiCor->Draw("HIST");

    c3->cd(1); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkPtRes->Draw();
    c3->cd(2); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkDxySig->Draw();
    c3->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkDzSig->Draw();
    c3->cd(5); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkNpixLayers->Draw();

    // Save canvas images
    const char *path = "./imgs/test/";
    const char *prefix = "all-histograms";
    const char *file_type = "png";

    save_canvas_images(canvases, numCanvases, path, prefix, file_type);

    // close file
    close_program(canvases, numCanvases, histograms, numHistograms, fr);
}