#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TF1.h>
#include <TMath.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TH1D.h"
#include <TStyle.h>
#include "TLegend.h"
#include <TText.h>
#include <TBenchmark.h>
#include <TLine.h>
#include <TFitResult.h>
#include <TSystem.h> // Include for gSystem if you need globbing
#include <TString.h> 

#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"
#include "../include/data_func.h"

TString findFile(const TString& pattern) {
    TString command = TString::Format("ls -1 %s 2>/dev/null | head -n 1", pattern.Data());

    return TString(gSystem->GetFromPipe(command)).Strip(TString::kBoth, '\n');
}

auto safeSumw2(TH1* h) {
    if (h && h->GetSumw2N() == 0) h->Sumw2();
}

void hist_comparison() {
    std::vector<int> selectionVars = {0, 10, 20, 60, 100, 140};
    int numCentBins = selectionVars.size() - 1;
    
    std::vector<std::string> histTypes = {"SS", "OS", "SSCor", "OSCor"};
    int numTypes = histTypes.size();

    TCanvas *c1 = new TCanvas("c1", "q_{inv} Comparison", 3840, 2160);
    TCanvas *canvases[] = {c1};
    int numCanvases = 1;
    
    c1->Divide(5, 4);

    for (int i = 0; i < numTypes; i++){
        for (int j = 0; j < numCentBins; j++){
            int cent_start = selectionVars[j];
            int cent_end = selectionVars[j+1];
            
            int cent_start_cor = selectionVars[j] / 2;
            int cent_end_cor = selectionVars[j+1] / 2;
            
            // Construct the expected filenames
            TString sl_pattern = TString::Format("data/Sig_mix/sigsl_cent_%d-%d*.root", cent_start, cent_end);
            TString sl_utn_pattern = TString::Format("data/Sig_mix/sig_sl_utn_cent_%d-%d*.root", cent_start, cent_end);
            TString dl_pattern = TString::Format("data/Sig_mix/sig_dl_cent_%d-%d*.root", cent_start, cent_end);
            
            TString sl_filename = findFile(sl_pattern);
            TString sl_utn_filename = findFile(sl_utn_pattern);
            TString dl_filename = findFile(dl_pattern);

            // Open the ROOT files and check if they exist
            TFile *sl_f = TFile::Open(sl_filename, "READ");
            TFile *sl_utn_f = TFile::Open(sl_utn_filename, "READ");
            TFile *dl_f = TFile::Open(dl_filename, "READ");

            // Check if files were opened successfully
            if (!sl_f || !sl_utn_f || !dl_f) {
                std::cerr << "Error: One or more files not found for centrality " 
                          << cent_start << "-" << cent_end << std::endl;
                std::cerr << "Attempted to open:" << std::endl;
                if (!sl_f) std::cerr << " - " << sl_filename << std::endl;
                if (!sl_utn_f) std::cerr << " - " << sl_utn_filename << std::endl;
                if (!dl_f) std::cerr << " - " << dl_filename << std::endl;
                
                // Close any files that might have been opened before returning
                if (sl_f) sl_f->Close();
                if (sl_utn_f) sl_utn_f->Close();
                if (dl_f) dl_f->Close();
                continue; // Skip to the next centrality bin
            }
            
            TH1D *sl_h_qinv_signal_1l = (TH1D *)sl_f->Get(TString::Format("h_qinv%s_signal_1l", histTypes[i].c_str()));
            TH1D *sl_utn_h_qinv_signal_1l = (TH1D *)sl_utn_f->Get(TString::Format("h_qinv%s_signal_1l", histTypes[i].c_str()));
            TH1D *dl_h_qinv_signal_1l = (TH1D *)dl_f->Get(TString::Format("h_qinv%s_signal_2l", histTypes[i].c_str()));

            safeSumw2(sl_h_qinv_signal_1l);
            safeSumw2(sl_utn_h_qinv_signal_1l);
            safeSumw2(dl_h_qinv_signal_1l);

            // Also check if histograms are retrieved successfully
            if (!sl_h_qinv_signal_1l || !sl_utn_h_qinv_signal_1l || !dl_h_qinv_signal_1l) {
                std::cerr << "Error: One or more histograms not found in files for centrality " 
                          << cent_start << "-" << cent_end << " and histType " << histTypes[i] << std::endl;
                // Close files before continuing to next iteration
                sl_f->Close();
                sl_utn_f->Close();
                dl_f->Close();
                continue; // Skip to the next iteration (next centrality bin or histType)
            }

            TH1D *h1 = (TH1D *)sl_h_qinv_signal_1l->Clone(TString::Format("h1_cent%d-%d_%s", cent_start, cent_end, histTypes[i].c_str()));    
            TH1D *h2 = (TH1D *)sl_h_qinv_signal_1l->Clone(TString::Format("h2_cent%d-%d_%s", cent_start, cent_end, histTypes[i].c_str()));    
            TH1D *h3 = (TH1D *)sl_utn_h_qinv_signal_1l->Clone(TString::Format("h3_cent%d-%d_%s", cent_start, cent_end, histTypes[i].c_str()));    

            // Check if division can be performed (denominators are not null and not empty)
            if (sl_utn_h_qinv_signal_1l->GetEntries() > 0) {
                h1->Divide(sl_utn_h_qinv_signal_1l);
            } else {
                std::cerr << "Warning: Division by zero for h1 (sl_utn_h_qinv_signal_1l is empty) in centrality " 
                << cent_start << "-" << cent_end << " and histType " << histTypes[i] << std::endl;
            }
            
            if (dl_h_qinv_signal_1l->GetEntries() > 0) {
                h2->Divide(dl_h_qinv_signal_1l);
                h3->Divide(dl_h_qinv_signal_1l);
            } else {
                std::cerr << "Warning: Division by zero for h2 or h3 (dl_h_qinv_signal_1l is empty) in centrality " 
                << cent_start << "-" << cent_end << " and histType " << histTypes[i] << std::endl;
            }
            
            TLegend *legend_c1 = new TLegend(0.2, 0.6, 0.9, 0.9);

            legend_c1->AddEntry(h1, "Single Loop / Single Loop UTN", "l");
            legend_c1->AddEntry(h2, "Single Loop / Double Loop", "l");
            legend_c1->AddEntry(h3, "Single Loop UTN / Double Loop", "l");

            h1->GetYaxis()->SetRangeUser(0.9, 1.1);
            h2->GetYaxis()->SetRangeUser(0.9, 1.1);
            h3->GetYaxis()->SetRangeUser(0.9, 1.1);
            
            Float_t tittle_size=0.09;
            Float_t marg = 0.16;
            
            h1->SetTitle("");
            h2->SetTitle("");
            h3->SetTitle("");
            
            if (i == 0 && j < 5) {
                h1->SetTitle(TString::Format("%d-%d%%", cent_start_cor, cent_end_cor));
                gStyle->SetTitleFontSize(tittle_size);
            }

            h1->SetYTitle("");
            h2->SetYTitle("");  
            h3->SetYTitle("");

            if (j==0) {
                h1->SetYTitle(histTypes[i].c_str());
                h1->SetTitleSize(tittle_size, "Y");
            }

            c1->cd(i*5 + j + 1  ); 
            gStyle->SetOptStat(0);

            h1->SetLineStyle(1);
            h1->Draw("HIST PLC");
            
            h2->SetLineStyle(2);
            h2->Draw("HIST PLC SAME");
            
            h3->SetLineStyle(3);
            h3->Draw("HIST PLC SAME");
            
            if (i == 0 && j == 0){
                legend_c1->Draw();
            }
            gStyle->SetPalette(kLake);
            if (j==0) {
                gPad->SetLeftMargin(marg);
            }
        }
    }

    // Saving image
    const char *path = "./imgs/test/hist_comparison";
    const char *prefix = "hist_comparison";
    const char *file_type = "png";
    const char *file_type2 = "pdf";
    save_canvas_images(canvases, numCanvases, path, prefix, file_type);
    save_canvas_images(canvases, numCanvases, path, prefix, file_type2);

    // Clean up canvases to avoid memory leaks
    delete c1;
}