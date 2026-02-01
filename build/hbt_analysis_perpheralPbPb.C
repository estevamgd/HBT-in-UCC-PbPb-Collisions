#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TF1.h>
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

    // Open file
    getFileTree("data/merged_2760PbPbMB_pixeltracks_UCC_skim.root", "demo/TreeMBUCC", fr, t);

    // --- Updated Variable Declarations ---
    const Int_t maxSize = 50000;
    Int_t Ntrk, hiBin;
    Int_t trkCharge[maxSize]; 
    
    Float_t HFsumET, pvZ;
    Float_t trkPt[maxSize], trkEta[maxSize], trkPhi[maxSize], trkWeight[maxSize];
    Float_t trkDxySig[maxSize], trkNpixLayers[maxSize], trkDzSig[maxSize];
    
    Float_t pionMass = 0.13957039; // Pion mass [GeV] from PDG

    // --- Setting Branch Addresses Explicitly ---
    t->SetBranchAddress("HFsumET", &HFsumET);
    t->SetBranchAddress("Ntrk", &Ntrk);
    t->SetBranchAddress("hiBin", &hiBin);
    t->SetBranchAddress("pvZ", &pvZ); 
    t->SetBranchAddress("trkCharge", trkCharge);
    t->SetBranchAddress("trkWeight", trkWeight);
    t->SetBranchAddress("trkPt", trkPt);
    t->SetBranchAddress("trkEta", trkEta);
    t->SetBranchAddress("trkPhi", trkPhi);
    t->SetBranchAddress("trkDxySig", trkDxySig);
    t->SetBranchAddress("trkNpixLayers", trkNpixLayers);
    t->SetBranchAddress("trkDzSig", trkDzSig);

    // Getting how many entries
    Long64_t nentries = t->GetEntries();

    // Setting canvases
    TCanvas *c1 = new TCanvas("c1", "Event Info", 7680, 4320);
    c1->Divide(3, 2);

    TCanvas *c2 = new TCanvas("c2", "Kinematics", 7680, 4320);
    c2->Divide(3, 2);

    TCanvas *c3 = new TCanvas("c3", "Track Quality", 7680, 4320);
    c3->Divide(3, 2);

    TCanvas *canvases[] = {c1, c2, c3};
    int numCanvases = 3;

    // Setting histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;

    // Event Level Histograms (C1)
    TH1D* hHFsumET = cHist("hHFsumET", "HFsumET[GeV]", "#Events/bin", 0, 4000, ninterval, nlength, nscale);
    TH1D* hNtrk = cHist("hNtrk", "Ntrk", "#Events/bin", 0, 4000, ninterval, nlength, nscale);
    TH1D* hHiBin = cHist("hHiBin", "hiBin", "#Events/bin", 0, 200, ninterval, nlength, nscale);
    TH1D* hPvZ = cHist("hPvZ", "pvZ [cm]", "#Events/bin", -30, 30, ninterval, nlength, nscale);
    TH1D* htrkCharge = cHist("htrkCharge", "trkCharge", "#Tracks/bin", -2, 2, ninterval, nlength, nscale);
    TH1D* htrkWeight = cHist("htrkWeight", "trkWeight", "#Tracks/bin", 0, 20, ninterval, nlength, nscale);

    // Track Kinematics Histograms (C2)
    TH1D* htrkPt = cHist("htrkPt", "pT[GeV]", "#Tracks/bin", 0, 45, ninterval, nlength, nscale);
    TH1D* htrkEta = cHist("htrkEta", "trkEta", "#Tracks/bin", -3, 3, ninterval, nlength, nscale);
    TH1D* htrkPhi = cHist("htrkPhi", "trkPhi", "#Tracks/bin", -3.4, 3.4, ninterval, nlength, nscale);
    
    // Corrected Histograms (C2)
    TH1D* htrkPtCor = cHist("htrkPtCor", "pT[GeV]", "#Tracks/bin", 0, 45, ninterval, nlength, nscale);
    TH1D* htrkEtaCor = cHist("htrkEtaCor", "trkEta", "#Tracks/bin", -3, 3, ninterval, nlength, nscale);
    TH1D* htrkPhiCor = cHist("htrkPhiCor", "trkPhi", "#Tracks/bin", -3.4, 3.4, ninterval, nlength, nscale);

    // Track Quality Histograms (C3)
    TH1D* htrkDxySig = cHist("htrkDxySig", "trkDxySig", "#Tracks/bin", -3, 3, ninterval, nlength, nscale);
    TH1D* htrkDzSig = cHist("htrkDzSig", "trkDzSig", "#Tracks/bin", -3.5, 3.5, ninterval, nlength, nscale);
    TH1D* htrkNpixLayers = cHist("htrkNpixLayers", "trkNpixLayers", "#Tracks/bin", 0, 6, ninterval, nlength, nscale);

    TH1D* histograms[] = {
        hHFsumET, hNtrk, hHiBin, hPvZ, htrkCharge, htrkWeight, 
        htrkPt, htrkEta, htrkPhi, htrkPtCor, htrkEtaCor, htrkPhiCor,
        htrkDxySig, htrkDzSig, htrkNpixLayers
    };
    int numHistograms = 15;

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        // Fill Event level variables
        //hHFsumET->Fill(HFsumET);
        //hNtrk->Fill(Ntrk);
        //hHiBin->Fill(hiBin);
        //hPvZ->Fill(pvZ);

        // Loop over tracks
        for (int j = 0; j < Ntrk; j++) {
            // Check weight validity
            if (!TMath::IsNaN(trkWeight[j]) && TMath::Finite(trkWeight[j])) {
                //htrkPtCor->Fill(trkPt[j], trkWeight[j]);
                //htrkEtaCor->Fill(trkEta[j], trkWeight[j]);
                //htrkPhiCor->Fill(trkPhi[j], trkWeight[j]);
            } 
            
            //htrkCharge->Fill(trkCharge[j]);
            //htrkWeight->Fill(trkWeight[j]);
            //htrkPt->Fill(trkPt[j]);
            //htrkEta->Fill(trkEta[j]);
            //htrkPhi->Fill(trkPhi[j]);
            //
            //htrkDxySig->Fill(trkDxySig[j]);
            //htrkDzSig->Fill(trkDzSig[j]);
            htrkNpixLayers->Fill(trkNpixLayers[j]);
        }
        
        if (i % 100 == 0) {
            float progress = (float)(i + 1) / nentries * 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
        }
    }
    std::cout<<std::endl;

    // Normalization clones
    TH1D *htrkPt_norm = (TH1D *)htrkPt->Clone("htrkPt_norm");
    TH1D *htrkPt_norm_cor = (TH1D *)htrkPtCor->Clone("htrkPt_norm_cor");
    TH1D *htrkEta_norm = (TH1D *)htrkEta->Clone("htrkEta_norm");
    TH1D *htrkEta_norm_cor = (TH1D *)htrkEtaCor->Clone("htrkEta_norm_cor");
    
    if (htrkPt_norm->Integral() > 0) htrkPt_norm->Scale(1.0/htrkPt_norm->Integral());
    if (htrkPt_norm_cor->Integral() > 0) htrkPt_norm_cor->Scale(1.0/htrkPt_norm_cor->Integral());
    if (htrkEta_norm->Integral() > 0) htrkEta_norm->Scale(1.0/htrkEta_norm->Integral());
    if (htrkEta_norm_cor->Integral() > 0) htrkEta_norm_cor->Scale(1.0/htrkEta_norm_cor->Integral());
    
    // Drawing histograms
    // Canvas 1: Event Variables
    c1->cd(1); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hHFsumET->Draw();
    c1->cd(2); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hHiBin->Draw();
    c1->cd(3); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hPvZ->Draw();
    c1->cd(4); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); hNtrk->Draw();
    c1->cd(5); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkCharge->Draw();
    c1->cd(6); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkWeight->Draw();
    
    // Canvas 2: Kinematics (Raw vs Weighted)
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

    // Canvas 3: Track Quality
    c3->cd(1); gPad->SetGrid(); gPad->SetLogy(); gPad->SetLeftMargin(0.15); htrkDxySig->Draw();
    c3->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkDzSig->Draw();
    c3->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.15); htrkNpixLayers->Draw();

    // Save canvas images
    const char *path = "./imgs/test/";
    const char *prefix = "all-histograms";
    const char *file_type = "png";

    TLegend *legends[] = {};
    int numLegends = 0;

    save_canvas_images(canvases, numCanvases, path, prefix, file_type);

    // close file
    close_program(canvases, numCanvases, histograms, numHistograms, legends, numLegends, fr);
}