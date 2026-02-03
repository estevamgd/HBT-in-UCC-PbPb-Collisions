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
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include <vector>
#include <string>

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


void plot_zdc_vs_hf(std::string input_file, std::string filename_out = "imgs/test/ZDC_vs_HF_UCC.pdf") {
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("demo/TreeMBUCC", input_file);

    auto df_filtered = df.Filter("HLT_HIUCC010_v2 == 1", "Trigger UCC");

    std::string nome_HF = "HFsumET";
    std::string nome_ZDC = "zdcSum";

    auto hist2d = df_filtered.Histo2D(
        {"zdc_vs_hf", ";HF #Sigma E_{T} [GeV];ZDC #Sigma E_{T} [GeV];#Events", 
         500, 0, 7000, 500, 0, 600000}, 
        nome_HF, 
        nome_ZDC
    );

    TCanvas *canvas = new TCanvas("c_zdc_hf", "", 1200, 800);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15); 
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.13);
    canvas->SetLogz();
    canvas->SetTickx(1);
    canvas->SetTicky(1);

    hist2d->SetStats(0);
    hist2d->GetXaxis()->SetTitleSize(0.05);
    hist2d->GetYaxis()->SetTitleSize(0.05);
    hist2d->GetZaxis()->SetTitleSize(0.05);
    hist2d->GetXaxis()->SetLabelSize(0.04);
    hist2d->GetYaxis()->SetLabelSize(0.04);
    hist2d->GetZaxis()->SetLabelSize(0.04);
    hist2d->GetXaxis()->SetTitleOffset(1.1);
    hist2d->GetYaxis()->SetTitleOffset(1.4);
    
    hist2d->Draw("COLZ");

    // --- Linha de Rejeição ---
    double x_start = 1450.0;
    double y_start = 600000.0;
    double x_end   = 3900.0;
    double y_end   = 0.0;
    double m = (y_end - y_start) / (x_end - x_start);
    double b = y_end - m * x_end;

    TLine *line = new TLine(x_start, y_start, x_end, y_end);
    line->SetLineColor(kRed);
    line->SetLineStyle(2); 
    line->SetLineWidth(4);
    line->Draw();

    // --- Legenda para a Linha ---
    // Posicionada no canto superior direito, abaixo do texto do sistema de colisão
    TLegend *leg = new TLegend(0.45, 0.75, 0.62, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // Transparente para não tapar o COLZ
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    // Adicionamos a linha manualmente na legenda usando o ponteiro 'line'
    leg->AddEntry(line, "Pile-Up Rejection Cut", "l");
    leg->AddEntry((TObject*)0, Form("Inclination (m) = %.2f", m), "");
    leg->AddEntry((TObject*)0, Form("Cut: ZDC < %.1f * HF + %.1f", m, b), "");
    leg->Draw();

    // --- Textos CMS ---
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(61);
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.18, 0.88, "CMS");
    
    latex.SetTextFont(52);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.28, 0.88, "Work in Progress");

    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.85, 0.94, "PbPb 2.76 TeV");

    canvas->SaveAs(filename_out.c_str());
    
    std::cout << "✅ Gráfico ZDC vs HF com legenda salvo em " << filename_out << std::endl;
}

void plot_hf_comparison_cms_style(std::string input_file, std::string filename_out = "imgs/test/HF_CMS_Style.pdf") {
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("demo/TreeMBUCC", input_file);

    // 1. Filtros e Definição da Reta
    auto df_triggered = df.Filter("HLT_HIUCC010_v2 == 1", "Trigger UCC");
    std::string cut_formula = "zdcSum < (-244.898 * HFsumET + 955102.0)";
    auto h_raw = df_triggered.Histo1D({"h_raw", "", 100, 0, 7000}, "HFsumET");
    auto h_clean = df_triggered.Filter(cut_formula).Histo1D({"h_clean", "", 100, 0, 7000}, "HFsumET");

    // 2. Configuração do Canvas (Estilo CMS)
    TCanvas *c = new TCanvas("c", "", 1200, 800);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.13);
    c->SetLogy();
    c->SetTickx(1);
    c->SetTicky(1);

    gStyle->SetOptStat(0);

    // 3. Estilização dos Histogramas (Pontos com Barras de Erro)
    auto style_histo = [](auto h, Color_t color, int marker) {
        h->SetMarkerStyle(marker);
        h->SetMarkerSize(1.2);
        h->SetMarkerColor(color);
        h->SetLineColor(color);
        h->SetLineWidth(1);
        h->SetStats(0);
        
        // Fontes e Títulos (Helvetica/Arial padrão CMS)
        h->GetXaxis()->SetTitle("HF #Sigma E_{T} [GeV]");
        h->GetYaxis()->SetTitle("#Events");
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.4);
    };

    style_histo(h_raw, kMagenta+2, 24); // Semelhante ao Double Ratio (q_inv)
    style_histo(h_clean, kOrange+1, 25); // Quadrado vazado, similar ao estilo da imagem

    h_raw->SetMinimum(0.5);
    h_raw->Draw("PE"); // P = Pontos, E = Erros estatísticos
    h_clean->Draw("PE SAME");

    // 4. Linha de Referência (Opcional, se quiser a linha tracejada em y=1 como na imagem)
    // TLine *line = new TLine(0, 1, 7000, 1);
    // line->SetLineStyle(2);
    // line->Draw();

    // 5. Legenda (Estilo Limpo)
    TLegend *leg = new TLegend(0.75, 0.70, 0.90, 0.85);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_raw.GetPtr(), "Raw Data", "pe");
    leg->AddEntry(h_clean.GetPtr(), "Pile-Up Cut", "pe");
    leg->Draw();

    // 6. Textos CMS (Lado Superior Esquerdo e Direito)
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(61); // Negrito
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.18, 0.88, "CMS");
    
    latex.SetTextFont(52); // Itálico
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.28, 0.88, "Work in Progress");
    latex.SetTextFont(42);
    latex.SetTextAlign(31); // Alinhado à direita
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.85, 0.94, "PbPb 2.76 TeV");


    c->SaveAs(filename_out.c_str());
}

int main() {
    // to compile use:
    // g++ -std=c++17 -pthread somePlots.cpp -o somePlots `root-config --cflags --libs`
    // to run use:
    // ./somePlots
    std::string input_file = "data/merged_2760PbPbMB_pixeltracks_UCC_skim.root";
    plot_zdc_vs_hf(input_file);
    plot_hf_comparison_cms_style(input_file);

    //qinvqlcmsPlot(0.0, 0.1, 0.9, 2.1);

    return 0;
}