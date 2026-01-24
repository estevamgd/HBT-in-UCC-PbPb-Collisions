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
#include "TLine.h"
#include "../include/my_func.h"
#include "../include/normalizer.h"
#include "../include/analyze_tools.h"

void qinvqlcmsPlot(double plotXMin, double plotXMax, double plotYMin, double plotYMax,
                   TH1D* hQinv_data = nullptr, TH1D* hQLCMS_data = nullptr) {
    
    TString searchPatternQinv, searchPatternQLCMS;
    const char* nameHist = "hRatioClone";
    if (hQinv_data == nullptr){
        searchPatternQinv = TString::Format("./data/correlation_ratios/sr_cor_qinv*.root");
        hQinv_data = getHistogram(searchPatternQinv.Data(), nameHist);
        if (hQinv_data == nullptr){
            std::cerr << "Error: Histogram " << nameHist << " not found in files matching " 
                      << searchPatternQinv.Data() << std::endl;
            return;
        }
    }
    if (hQLCMS_data == nullptr){
        searchPatternQLCMS = TString::Format("./data/correlation_ratios/sr_cor_qlcms*.root");
        hQLCMS_data = getHistogram(searchPatternQLCMS.Data(), nameHist);
        if (hQLCMS_data == nullptr){
            std::cerr << "Error: Histogram " << nameHist << " not found in files matching " 
                      << searchPatternQLCMS.Data() << std::endl;
            return;
        }
    }

    TCanvas *cComp = new TCanvas("cComp", "Comparison qinv x qlcm", 1200, 800);
    gStyle->SetOptStat(0);
    cComp->cd();
    cComp->SetLeftMargin(0.12);
    cComp->SetBottomMargin(0.12);

    hQinv_data->SetMarkerStyle(24);
    hQinv_data->SetMarkerColor(colors[0]);
    hQinv_data->SetLineColor(colors[0]-7);
    hQinv_data->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
    hQinv_data->GetYaxis()->SetRangeUser(plotYMin, plotYMax);
    hQinv_data->GetYaxis()->SetTitleOffset(1.2);
    hQinv_data->SetTitle("; q [GeV]; C_{2}(q) = Data/Fit");
    hQinv_data->Draw("E1");

    hQLCMS_data->SetMarkerStyle(25);
    hQLCMS_data->SetMarkerColor(colors[1]);
    hQLCMS_data->SetLineColor(kOrange);
    hQLCMS_data->Draw("E1 SAME");

    TLine *line = new TLine(plotXMin, 1.0, plotXMax, 1.0); 

    line->SetLineColor(kGray + 2);
    line->SetLineStyle(kDashed);
    line->SetLineWidth(2);
    line->Draw("SAME");

    TLegend *legendComp = new TLegend(0.55, 0.60, 0.83, 0.80);
    legendComp->AddEntry((TObject*)nullptr, "3200 < HF #Sigma E_{T} < 3300", "");
    legendComp->AddEntry((TObject*)nullptr, "pT > 0.5 GeV", "");
    legendComp->AddEntry(hQinv_data, "Double Ratio (q_{inv})", "lep");
    legendComp->AddEntry(hQLCMS_data, "Double Ratio (q_{LCMS})", "lep");
    legendComp->SetFillStyle(0);
    legendComp->SetBorderSize(0);
    legendComp->Draw();

    drawCMSHeaders("#bf{CMS} #it{Work in Progress}", "PbPb 2.76 TeV | q inv x q LCMS Comparison", 1.0);

    const char* imagePath = "./imgs/test/correlation_ratios/";
    TCanvas *canvasesToSave[] = { cComp };
    save_canvas_images(canvasesToSave, 1, imagePath, "qinv_qlcms_comparison", "png");
    save_canvas_images(canvasesToSave, 1, imagePath, "qinv_qlcms_comparison", "pdf");
    std::cout << "Comparison plot saved to " << imagePath << std::endl;
    
    AnalysisLog::instance().save("./logs", "qinvqlcmsPlot");

    delete legendComp;
    delete cComp;
}

int main() {
    // to compile use:
    // g++ -std=c++17 -pthread somePlots.cpp -o qinvqlcmsPlot `root-config --cflags --libs`
    // to run use:
    // ./qinvqlcmsPlot

    qinvqlcmsPlot(0.0, 0.3, 0.9, 2.1);

    return 0;
}