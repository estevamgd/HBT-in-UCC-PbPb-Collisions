#ifndef FITS_H
#define FITS_H
#include <iostream>
#include <string>
#include <vector>
#include <filesystem> 
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h" 
#include "TSystem.h" 
#include "Math/MinimizerOptions.h"
#include "../include/my_func.h" 
#include "../include/data_func.h"


struct FitOutput {
    FitFunctionType type;
    TF1* function;
    TFitResultPtr result;
    std::vector<FitParamInfo> params;
    std::string displayName;
};


FitModelConfig getFitModelConfig(FitFunctionType type)
{
    switch (type) {

    case FitFunctionType::EXPONENTIAL:
        return {
            "fitExp", "Exponential Fit", FitExp, 4,
            {"Const","#lambda","R (fm)","#epsilon"},
            {{0.,0.},{0.,1.},{0.,0.},{0.,0.}},
            {{1, "#lambda", ""},{2, "R", " fm"}}
        };

    case FitFunctionType::GAUSSIAN:
        return {
            "fitGauss", "Gaussian Fit", FitGauss, 4,
            {"Const","#lambda","R (fm)","#epsilon"},
            {{0.,0.},{0.,1.},{0.,0.},{0.,0.}},
            {{1, "#lambda", ""},{2, "R", " fm"}}
        };

    case FitFunctionType::LEVY:
        return {
            "fitLevy", "Levy Fit", FitLevy, 5,
            {"Const","#lambda","R (fm)","#epsilon","#alpha"},
            {{0.,0.},{0.,1.},{0.,0.},{0.,0.},{1.,2.}},
            {{1, "#lambda", ""},{2, "R", " fm"},{4, "#alpha", ""}}
        };

    case FitFunctionType::LEVY2:
        return {
            "fitLevy2", "Levy2 Fit", FitLevy2, 3,
            {"#lambda","R (fm)","#alpha"},
            {{0.,1.},{0.,0.},{1.,2.}},
            {{0, "#lambda", ""},{1, "R", " fm"},{2, "#alpha", ""}}
        };

    case FitFunctionType::DOUBLE_LEVY:
        return {
            "fitLevyDR", "Double Ratio Levy Fit", FitLevyDR, 5,
            {"#lambda","R","#alpha","#epsilon","Norm"},
            {{0.,1.},{0.,0.},{1.,2.},{0.,0.},{0.,0.}},
            {{0, "#lambda", ""},{1, "R", " fm"},{2, "#alpha", ""}}
        };
    
    case FitFunctionType::BACKGROUND:
        return {
            "fitBG", "Background Fit", FitBG, 5,
            {"Norm","a","b","c","d"},
            {{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}},
            {}
        };
    }

    

    throw std::runtime_error("Unknown FitFunctionType");
}


TF1* fitHistogram(
    TH1D* hist,
    FitFunctionType type,
    const FitInit& init,
    double fitMin,
    double fitMax,
    TFitResultPtr* fitResult
) {
    auto cfg = getFitModelConfig(type);
    
    if ((int)init.values.size() != cfg.nPar)
        throw std::runtime_error("Wrong number of initial parameters");

    TF1* f = new TF1(cfg.name, cfg.func, fitMin, fitMax, cfg.nPar);

    for (int i = 0; i < cfg.nPar; ++i) {
        f->SetParameter(i, init.values[i]);
        f->SetParName(i, cfg.parNames[i].c_str());

        if (i < (int)cfg.parLimits.size() &&
            cfg.parLimits[i].first != cfg.parLimits[i].second){
                f->SetParLimits(i,
                    cfg.parLimits[i].first,
                    cfg.parLimits[i].second);
                std::cout << "Set limits for parameter " << i << ": "
                          << cfg.parLimits[i].first << " to "
                          << cfg.parLimits[i].second << std::endl;
            }
    }

    TFitResultPtr res = hist->Fit(f, "S R E M");
    if (fitResult) *fitResult = res;

    return f;
}

std::vector<FitOutput> fitHistogramMultiple(
    TH1D* hist,
    const std::vector<std::pair<FitFunctionType, FitInit>>& models,
    double fitMin,
    double fitMax
) {
    std::vector<FitOutput> outputs;

    for (const auto& [type, init] : models) {
        TFitResultPtr res;
        TF1* f = fitHistogram(hist, type, init, fitMin, fitMax, &res);
        auto cfg = getFitModelConfig(type);
        outputs.push_back({type, f, res, cfg.legendParams, cfg.displayName});
    }

    return outputs;
}

void drawAndSaveFits(
    TH1D* hist,
    const std::vector<FitOutput>& fits,
    const char* canvasName,
    TString headerLabel,
    const char* outputPrefix,
    qMode mode = qMode::QINV,
    double plotXMin = 0.0, double plotXMax = 10.0,
    double plotYMin = 0.9, double plotYMax = 2.1) {
    TH1D* histClone = (TH1D*)(hist->Clone("histClone"));
    gStyle->SetOptStat(0);     
    gStyle->SetPalette(kBlueRedYellow);
    TCanvas* c = new TCanvas(canvasName, "Correlation Fits", 1200, 900);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    histClone->SetMarkerStyle(20);
    histClone->SetMarkerColor(kBlack);
    histClone->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
    histClone->GetYaxis()->SetRangeUser(plotYMin, plotYMax);
    histClone->SetXTitle(Form("q_{%s} [GeV]", 
        (mode == qMode::QLCMS) ? "LCMS" : "inv"));
    histClone->SetYTitle(Form("C(q_{%s}) = %s", 
        (mode == qMode::QLCMS) ? "LCMS" : "inv",
        (headerLabel.Contains("Double Ratio")) ? "Data/Fit" : "SS/OS"));
    histClone->Draw("E1 P");
    
    TLegend* leg = new TLegend(0.55, 0.40, 0.88, 0.88);
    leg->SetTextFont(42);
    leg->SetTextSize(0.022);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->SetEntrySeparation(0.12);
    leg->SetMargin(0.18);
    
    leg->AddEntry(histClone, "Data", "lep");
    leg->AddEntry((TObject*)nullptr, " ", "");
    
    Int_t fitIndex = 0;
    for (const auto& fit : fits) {
        fit.function->Draw("SAME");
        fit.function->SetLineColor(colors[fitIndex]);
         fit.function->SetLineStyle(lineStyles[1]);
        fit.function->SetLineWidth(2);

        leg->AddEntry(fit.function, fit.displayName.c_str(), "l");
        
        for (const auto& p : fit.params) {
        leg->AddEntry(
                (TObject*)nullptr,
                Form(" %s = %.2f #pm %.2f%s",
                     p.label.c_str(),
                     fit.function->GetParameter(p.index),
                     fit.function->GetParError(p.index),
                     p.unit.c_str()),
                ""
            );
        }
        
        leg->AddEntry(
            (TObject*)nullptr,
            Form("#chi^{2}/NDF = %.1f / %d",
                 fit.result->Chi2(),
                 fit.result->Ndf()),
            ""
        );

        leg->AddEntry(
            (TObject*)nullptr,
            Form("p-value = %.5f", fit.result->Prob()),
            ""
        );

        leg->AddEntry((TObject*)nullptr, " ", "");
        fitIndex++;
    }
    
    leg->Draw();
    double sizeFactor = 1.0;
    if (headerLabel.Length() > 40) sizeFactor = 0.8;
    if (headerLabel.Length() > 60) sizeFactor = 0.7;
    drawCMSHeaders("#bf{CMS} #it{Preliminary}", headerLabel, sizeFactor);

    TCanvas* canvases[] = {c};
    
    save_histograms(&histClone, 1, "./data/fit_correlation/", outputPrefix);
    save_canvas_images(canvases, 1, "./imgs/test/fit_correlation/", outputPrefix, "png");
    save_canvas_images(canvases, 1, "./imgs/test/fit_correlation/", outputPrefix, "pdf");

    delete leg;
    delete c;
}

#endif