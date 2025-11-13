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
#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"

void coulomb_correction() {
    // Setting up, copy from HERE to...
    // Getting TTree
    TFile *fr = nullptr;
    TTree *t = nullptr;

    getFileTree("data/HiForestAOD_UPC.root", "demo/HBT", fr, t);
    
    // Variables
    Int_t maxSize = 17100;
    Int_t Ntrk = maxSize, NSSpair = maxSize, NOSpair = maxSize;
    Float_t HFsumET;
    Float_t coulombWOS[NOSpair], coulombWSS[NSSpair], qinvSigSS[NSSpair], qinvSigOS[NOSpair], 
            trkPt[Ntrk], trkEta[Ntrk], trkPhi[Ntrk], trkPtRes[Ntrk], 
            trkDxySig[Ntrk], trkNpixLayers[Ntrk], trkDzSig[Ntrk];

    // Arrays of variables
    void* variables[] = {
        &HFsumET, &Ntrk, &NSSpair, &NOSpair, coulombWOS, coulombWSS, qinvSigOS, qinvSigSS, 
        trkPt, trkEta, trkPhi, trkPtRes, trkDxySig, trkDzSig, trkNpixLayers
    };

    // Branch names
    const char* branchNames[] = {
        "HFsumET", "Ntrk", "NSSpair", "NOSpair", "coulombWOS", "coulombWSS", 
        "qinvSigOS", "qinvSigSS", "trkPt", "trkEta", "trkPhi", "trkPtRes", 
        "trkDxySig", "trkDzSig", "trkNpixLayers"
    };

    // Setting addresses
    int numBranches = sizeof(variables) / sizeof(variables[0]);
    for (int i = 0; i < numBranches; i++) {
        t->SetBranchAddress(branchNames[i], variables[i]);
    }

    // ...HERE
    // Getting how many entries
    Long64_t nentries = t->GetEntries();
    double d_nentries = t->GetEntries();

    // Setting canvases
    TCanvas *c1 = new TCanvas("c1", "", 7680, 4320);
    TCanvas *c2 = new TCanvas("c2", "", 7680, 4320);
    TCanvas *c3 = new TCanvas("c3", "", 7680, 4320);

    TCanvas *canvases[] = {c1, c2, c3};
    int numCanvases = 3;

    // Setting histograms
    double ninterval = 1., nlength = 0.02, nscale = 1./1.;

    double nnscale = numBins(ninterval, nlength, nscale);
    TH1D* h11 = new TH1D("h11", "", nnscale, 0., 1.);
    TH1D* h22 = new TH1D("h22", "", nnscale, 0., 1.);
    TH1D* h33 = new TH1D("h33", "", nnscale, 0., 1.);
    TH1D* h44 = new TH1D("h44", "", nnscale, 0., 1.);
    /*
    TH1D *h1 = cHist("h1", "qinv[GeV]", "#Pairs/bin", 0., 1., ninterval, nlength, nscale);
    TH1D *h2 = cHist("h2", "qinv[GeV]", "#Pairs/bin", 0., 1., ninterval, nlength, nscale);
    TH1D *h3 = cHist("h3", "", "", 0., 1., ninterval, nlength, nscale);
    TH1D *h4 = cHist("h4", "", "", 0., 1., ninterval, nlength, nscale);
    */
    TH1D *tonormhist[] = {h11, h22, h33, h44};
    int numTonorm = 4;

    // Filling histograms
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        if (HFsumET > 100 && HFsumET < 375) { // This selects a centrality
            for (int k = 0; k < NOSpair; k++) {
                h33->Fill(qinvSigOS[k], coulombWOS[k]);
                h11->Fill(qinvSigOS[k]);
            }
            for (int l = 0; l < NSSpair; l++) {
                h44->Fill(qinvSigSS[l], coulombWSS[l]);
                h22->Fill(qinvSigSS[l]);
            }    
        }

        // Display progress
        float progress = (float)(i + 1) / nentries * 100;
        std::cout << "\rFilling Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
    }

    // Normalize histograms
    Double_t scale = 1;
    normalizer(tonormhist, numTonorm, scale);
    
    // Dividing SS/OS2
    TH1D *sr = (TH1D *)h22->Clone("sr");
    TH1D *sr_cor = (TH1D *)h44->Clone("sr_cor");
    sr->Divide(h11);
    sr_cor->Divide(h33);

    // Adding Title
    sr->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV");
    h11->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV");
    h22->SetTitle("CMS Open Data 2011 - PbPb 2.76 TeV");

    // Adjusting fill colors
    h11->SetFillColorAlpha(kRed, 0.5);
    h22->SetFillColorAlpha(kRed, 0.5);
    h33->SetFillColorAlpha(kBlue, 0.5);
    h44->SetFillColorAlpha(kBlue, 0.5);
    sr->SetFillColorAlpha(kRed, 0.5);
    sr_cor->SetFillColorAlpha(kBlue, 0.5);

    // Adjusting fill colors
    h11->SetLineColor(kRed);
    h22->SetLineColor(kRed);
    h33->SetLineColor(kBlue);
    h44->SetLineColor(kBlue);
    sr->SetLineWidth(4);
    sr->SetLineColor(kRed);
    sr_cor->SetLineWidth(4);
    sr_cor->SetLineColor(kBlue);

    // Set titles and axis labels
    h11->GetXaxis()->SetTitle("qinv[GeV]");
    h11->GetYaxis()->SetTitle("#Pairs/bin");

    h22->GetXaxis()->SetTitle("qinv[GeV]");
    h22->GetYaxis()->SetTitle("#Pairs/bin");

    h33->GetXaxis()->SetTitle("qinv[GeV]");
    h33->GetYaxis()->SetTitle("#Pairs/bin");

    h44->GetXaxis()->SetTitle("qinv[GeV]");
    h44->GetYaxis()->SetTitle("#Pairs/bin");

    sr->GetXaxis()->SetTitle("qinv[GeV]");
    sr->GetYaxis()->SetTitle("Single Ratio SS/OS");

    sr_cor->GetXaxis()->SetTitle("qinv[GeV]");
    sr_cor->GetYaxis()->SetTitle("Single Ratio SS/OS");
    
    // Set label sizes
    h11->GetXaxis()->SetLabelSize(0.04);
    h22->GetYaxis()->SetLabelSize(0.04);
    h33->GetYaxis()->SetLabelSize(0.04);
    h44->GetYaxis()->SetLabelSize(0.04);
    sr->GetYaxis()->SetLabelSize(0.04);
    sr_cor->GetYaxis()->SetLabelSize(0.04);

    // Adjust title font size and offset
    /*sr->GetXaxis()->SetTitleSize(0.05);
    sr->GetYaxis()->SetTitleSize(0.05);
    */
    // Single Ratio Legend
    TLegend *h1_legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    h1_legend->AddEntry((TObject*)0, "Invariant 4-momentum", "");
    h1_legend->AddEntry((TObject*)0, "Opposite charge pairs", "");
    h1_legend->AddEntry((TObject*)0, "Centrality: 50-70%", "");
    h1_legend->AddEntry(h11, "Without Coulomb Correction", "l");
    h1_legend->AddEntry(h33, "With Coulomb Correction", "l");

    // Single Ratio Legend
    TLegend *h2_legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    h2_legend->AddEntry((TObject*)0, "Invariant 4-momentum", "");
    h2_legend->AddEntry((TObject*)0, "Same charge pairs", "");
    h2_legend->AddEntry((TObject*)0, "Centrality: 50-70%", "");
    h2_legend->AddEntry(h11, "Without Coulomb Correction", "l");
    h2_legend->AddEntry(h33, "With Coulomb Correction", "l");

    // Single Ratio Legend
    TLegend *sr_legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    sr_legend->AddEntry((TObject*)0, "Single Ratio SS/OS", "");
    sr_legend->AddEntry((TObject*)0, "Centrality: 50-70%", "");
    sr_legend->AddEntry(sr, "Without Coulomb Correction", "l");
    sr_legend->AddEntry(sr_cor, "With Coulomb Correction", "l");
    
    // Setting y range in single ratio to 0.95<y<1.6
    sr->GetYaxis()->SetRangeUser(0.95, 1.4);
    sr_cor->GetYaxis()->SetRangeUser(0.95, 1.4);
    /*
    // Setting x range in qinv to qinv<0.1
    h1->GetXaxis()->SetRangeUser(0,0.1);
    h2->GetXaxis()->SetRangeUser(0,0.1);
    h3->GetXaxis()->SetRangeUser(0,0.1);
    h4->GetXaxis()->SetRangeUser(0,0.1);
    */
    /*
    h11->SetFillStyle(1001);
    h22->SetFillStyle(1001);
    h33->SetFillStyle(1001);
    h44->SetFillStyle(1001);
    sr->SetFillStyle(1001);
    sr_cor->SetFillStyle(1001);
    */
    

    // Removing statistics box
    TH1D *histograms[] = {h11, h22, h33, h44, sr, sr_cor};
    int numHistograms = 6;

    no_statbox(histograms, numHistograms);

    // sr->SetMarkerStyle(4);
    // sr_cor->SetMarkerStyle(4);
   
    // Drawing
    c1->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); h33->Draw("HIST"); h11->Draw("HIST SAME"); h1_legend->Draw();
    c2->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); h44->Draw("HIST"); h22->Draw("HIST SAME"); h2_legend->Draw();
    c3->cd(); gPad->SetGrid(); gPad->SetLeftMargin(0.15); sr_cor->Draw("HIST"); sr->Draw("HIST SAME"); sr_legend->Draw();
    
    // Saving image
    const char *path = "./imgs/final/";
    const char *prefix = "final-coulomb-correction";
    const char *file_type = "png";
    save_canvas_images(canvases, numCanvases, path, prefix, file_type);

    // Closing program
    close_program(canvases, numCanvases, histograms, numHistograms, fr);
}
