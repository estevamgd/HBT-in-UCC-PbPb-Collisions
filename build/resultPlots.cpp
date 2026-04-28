#include "../include/ratiosAndFits.h"
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TH1F.h>


void resultPlots(double bins[], size_t numBins){
    ControlVar selectedControlVar = ControlVar::CENTHF;
    qMode modeLCMS = qMode::QLCMS;
    qMode modeQinv = qMode::QINV;
    
    Double_t etaMin = 0.95;
    Double_t etaMax = -0.95;
    Double_t ptMin  = 0.50;

    std::vector<std::pair<FitFunctionType, FitInit>> models = {
        {FitFunctionType::EXPONENTIAL, {{0.6, 4.0, 0.0, 0.0, 1.0}}},
        {FitFunctionType::GAUSSIAN, {{0.6, 4.0, 0.0, 0.0, 1.0}}},
        {FitFunctionType::DOUBLE_LEVY, {{0.6, 4.0, 1.5, 0.0, 1.0}}}
    };
    
    double plotXMin = 0.0;
    double plotXMax = 0.5;

    double plotYMin = 0.9;
    double plotYMax = 2.1;
    
    double fitMin = 0.04;
    double fitMax = 0.2;
    double fitMinBg = 0.2;

    // Normalization qinv range
    Double_t q1 = 6.82;
    Double_t q2 = 8.4;

    std::vector<std::pair<double, double>> sr_y_values_exp, sr_y_values_gaus, sr_y_values_levy;
    std::vector<std::pair<double, double>> dr_y_values_exp, dr_y_values_gaus, dr_y_values_levy;
    
    std::vector<std::pair<double, double>> sr_y_erros_exp, sr_y_erros_gaus, sr_y_erros_levy;
    std::vector<std::pair<double, double>> dr_y_erros_exp, dr_y_erros_gaus, dr_y_erros_levy;

    std::vector<std::string> bin_labels;

    for (size_t i = 0; i + 1 < numBins; i++) {
        double bin_low = bins[i];
        double bin_high = bins[i+1];

        std::string binLabel = std::to_string(static_cast<int>(bin_low)) + "-" + std::to_string(static_cast<int>(bin_high));
        
        MultipleFitResults fitResults;
        fitResults = doubleRatioMixFit(
            selectedControlVar,
            models,
            bin_low, bin_high,
            etaMin, etaMax, ptMin,
            q1, q2,
            modeLCMS,
            fitMin, fitMax, fitMinBg,
            plotXMin, plotXMax, plotYMin, plotYMax
        );

        int j = 0;
        for (const auto& res : fitResults.results) {
            double lambda = res->Parameter(0);
            double R = res->Parameter(1);

            double lambdaErr = res->ParError(0);
            double RErr = res->ParError(1);

            if (fitResults.modelNames[j].find("sr") != std::string::npos) {
                if (fitResults.modelNames[j].find("Exponential") != std::string::npos) {
                    sr_y_values_exp.push_back({R, lambda});
                    sr_y_erros_exp.push_back({RErr, lambdaErr});
                } else if (fitResults.modelNames[j].find("Gaussian") != std::string::npos) {
                    sr_y_values_gaus.push_back({R, lambda});
                    sr_y_erros_gaus.push_back({RErr, lambdaErr});
                } else if (fitResults.modelNames[j].find("Levy") != std::string::npos) {
                    sr_y_values_levy.push_back({R, lambda});
                    sr_y_erros_levy.push_back({RErr, lambdaErr});
                }
            } else if (fitResults.modelNames[j].find("dr") != std::string::npos) {
                if (fitResults.modelNames[j].find("Exponential") != std::string::npos) {
                    dr_y_values_exp.push_back({R, lambda});
                    dr_y_erros_exp.push_back({RErr, lambdaErr});
                } else if (fitResults.modelNames[j].find("Gaussian") != std::string::npos) {
                    dr_y_values_gaus.push_back({R, lambda});
                    dr_y_erros_gaus.push_back({RErr, lambdaErr});
                } else if (fitResults.modelNames[j].find("Levy") != std::string::npos) {
                    dr_y_values_levy.push_back({R, lambda});
                    dr_y_erros_levy.push_back({RErr, lambdaErr});
                }
            }
            j++;
        }
        bin_labels.push_back(binLabel);
    }

    int n = static_cast<int>(bin_labels.size());
    std::vector<double> x(n), xErr(n, 0.0);
    for (int i = 0; i < n; i++) {
        x[i] = i + 1.0;
    }

    TCanvas* c1_R = new TCanvas("c1_R", "Single Ratio", 1200, 800);
    TCanvas* c1_lambda = new TCanvas("c1_lambda", "Single Lambda", 1200, 800);
    TCanvas* c2_R = new TCanvas("c2_R", "Double Ratio", 1200, 800);
    TCanvas* c2_lambda = new TCanvas("c2_lambda", "Double Lambda", 1200, 800);
    c1_R->SetTopMargin(0.10);
    c1_lambda->SetTopMargin(0.10);
    c2_R->SetTopMargin(0.10);
    c2_lambda->SetTopMargin(0.10);
    c1_R->SetBottomMargin(0.12);
    c1_lambda->SetBottomMargin(0.12);
    c2_R->SetBottomMargin(0.12);
    c2_lambda->SetBottomMargin(0.12);

    std::vector<double> sr_R_exp, sr_R_gaus, sr_R_levy;
    std::vector<double> sr_lambda_exp, sr_lambda_gaus, sr_lambda_levy;
    std::vector<double> dr_R_exp, dr_R_gaus, dr_R_levy;
    std::vector<double> dr_lambda_exp, dr_lambda_gaus, dr_lambda_levy;

    std::vector<double> sr_R_err_exp, sr_R_err_gaus, sr_R_err_levy;
    std::vector<double> sr_lambda_err_exp, sr_lambda_err_gaus, sr_lambda_err_levy;
    std::vector<double> dr_R_err_exp, dr_R_err_gaus, dr_R_err_levy;
    std::vector<double> dr_lambda_err_exp, dr_lambda_err_gaus, dr_lambda_err_levy;

    for (int i = 0; i < n; i++) {
        sr_R_exp.push_back(sr_y_values_exp[i].first);
        sr_R_gaus.push_back(sr_y_values_gaus[i].first);
        sr_R_levy.push_back(sr_y_values_levy[i].first);
        sr_lambda_exp.push_back(sr_y_values_exp[i].second);
        sr_lambda_gaus.push_back(sr_y_values_gaus[i].second);
        sr_lambda_levy.push_back(sr_y_values_levy[i].second);

        dr_R_exp.push_back(dr_y_values_exp[i].first);
        dr_R_gaus.push_back(dr_y_values_gaus[i].first);
        dr_R_levy.push_back(dr_y_values_levy[i].first);
        dr_lambda_exp.push_back(dr_y_values_exp[i].second);
        dr_lambda_gaus.push_back(dr_y_values_gaus[i].second);
        dr_lambda_levy.push_back(dr_y_values_levy[i].second);

        sr_R_err_exp.push_back(sr_y_erros_exp[i].first);
        sr_R_err_gaus.push_back(sr_y_erros_gaus[i].first);
        sr_R_err_levy.push_back(sr_y_erros_levy[i].first);
        sr_lambda_err_exp.push_back(sr_y_erros_exp[i].second);
        sr_lambda_err_gaus.push_back(sr_y_erros_gaus[i].second);
        sr_lambda_err_levy.push_back(sr_y_erros_levy[i].second);

        dr_R_err_exp.push_back(dr_y_erros_exp[i].first);
        dr_R_err_gaus.push_back(dr_y_erros_gaus[i].first);
        dr_R_err_levy.push_back(dr_y_erros_levy[i].first);
        dr_lambda_err_exp.push_back(dr_y_erros_exp[i].second);
        dr_lambda_err_gaus.push_back(dr_y_erros_gaus[i].second);
        dr_lambda_err_levy.push_back(dr_y_erros_levy[i].second);
    }

    TGraphErrors* gr_sr_R_exp = new TGraphErrors(n, x.data(), sr_R_exp.data(), xErr.data(), sr_R_err_exp.data());
    TGraphErrors* gr_sr_R_gaus = new TGraphErrors(n, x.data(), sr_R_gaus.data(), xErr.data(), sr_R_err_gaus.data());
    TGraphErrors* gr_sr_R_levy = new TGraphErrors(n, x.data(), sr_R_levy.data(), xErr.data(), sr_R_err_levy.data());
    
    TGraphErrors* gr_sr_lambda_exp = new TGraphErrors(n, x.data(), sr_lambda_exp.data(), xErr.data(), sr_lambda_err_exp.data());
    TGraphErrors* gr_sr_lambda_gaus = new TGraphErrors(n, x.data(), sr_lambda_gaus.data(), xErr.data(), sr_lambda_err_gaus.data());
    TGraphErrors* gr_sr_lambda_levy = new TGraphErrors(n, x.data(), sr_lambda_levy.data(), xErr.data(), sr_lambda_err_levy.data());
    
    TGraphErrors* gr_dr_R_exp = new TGraphErrors(n, x.data(), dr_R_exp.data(), xErr.data(), dr_R_err_exp.data());
    TGraphErrors* gr_dr_R_gaus = new TGraphErrors(n, x.data(), dr_R_gaus.data(), xErr.data(), dr_R_err_gaus.data());
    TGraphErrors* gr_dr_R_levy = new TGraphErrors(n, x.data(), dr_R_levy.data(), xErr.data(), dr_R_err_levy.data());
    
    TGraphErrors* gr_dr_lambda_exp = new TGraphErrors(n, x.data(), dr_lambda_exp.data(), xErr.data(), dr_lambda_err_exp.data());
    TGraphErrors* gr_dr_lambda_gaus = new TGraphErrors(n, x.data(), dr_lambda_gaus.data(), xErr.data(), dr_lambda_err_gaus.data());
    TGraphErrors* gr_dr_lambda_levy = new TGraphErrors(n, x.data(), dr_lambda_levy.data(), xErr.data(), dr_lambda_err_levy.data());

    gr_sr_R_exp->SetTitle("R - Single Ratio");
    gr_sr_lambda_exp->SetTitle("Lambda - Single Ratio");
    gr_dr_R_exp->SetTitle("R - Double Ratio");
    gr_dr_lambda_exp->SetTitle("Lambda - Double Ratio");

    gr_sr_R_exp->SetMarkerStyle(20);
    gr_sr_R_gaus->SetMarkerStyle(21);
    gr_sr_R_levy->SetMarkerStyle(22);
    gr_sr_lambda_exp->SetMarkerStyle(20);
    gr_sr_lambda_gaus->SetMarkerStyle(21);
    gr_sr_lambda_levy->SetMarkerStyle(22);
    gr_dr_R_exp->SetMarkerStyle(20);
    gr_dr_R_gaus->SetMarkerStyle(21);
    gr_dr_R_levy->SetMarkerStyle(22);
    gr_dr_lambda_exp->SetMarkerStyle(20);
    gr_dr_lambda_gaus->SetMarkerStyle(21);
    gr_dr_lambda_levy->SetMarkerStyle(22);

    gr_sr_R_exp->SetMarkerColor(kRed + 1);
    gr_sr_R_gaus->SetMarkerColor(kBlue + 1);
    gr_sr_R_levy->SetMarkerColor(kGreen + 2);
    gr_sr_lambda_exp->SetMarkerColor(kRed + 1);
    gr_sr_lambda_gaus->SetMarkerColor(kBlue + 1);
    gr_sr_lambda_levy->SetMarkerColor(kGreen + 2);
    gr_dr_R_exp->SetMarkerColor(kRed + 1);
    gr_dr_R_gaus->SetMarkerColor(kBlue + 1);
    gr_dr_R_levy->SetMarkerColor(kGreen + 2);
    gr_dr_lambda_exp->SetMarkerColor(kRed + 1);
    gr_dr_lambda_gaus->SetMarkerColor(kBlue + 1);
    gr_dr_lambda_levy->SetMarkerColor(kGreen + 2);

    gr_sr_R_exp->SetLineColor(kRed + 1);
    gr_sr_R_gaus->SetLineColor(kBlue + 1);
    gr_sr_R_levy->SetLineColor(kGreen + 2);
    gr_sr_lambda_exp->SetLineColor(kRed + 1);
    gr_sr_lambda_gaus->SetLineColor(kBlue + 1);
    gr_sr_lambda_levy->SetLineColor(kGreen + 2);
    gr_dr_R_exp->SetLineColor(kRed + 1);
    gr_dr_R_gaus->SetLineColor(kBlue + 1);
    gr_dr_R_levy->SetLineColor(kGreen + 2);
    gr_dr_lambda_exp->SetLineColor(kRed + 1);
    gr_dr_lambda_gaus->SetLineColor(kBlue + 1);
    gr_dr_lambda_levy->SetLineColor(kGreen + 2);

    gr_sr_R_exp->SetLineWidth(2);
    gr_sr_R_gaus->SetLineWidth(2);
    gr_sr_R_levy->SetLineWidth(2);
    gr_sr_lambda_exp->SetLineWidth(2);
    gr_sr_lambda_gaus->SetLineWidth(2);
    gr_sr_lambda_levy->SetLineWidth(2);
    gr_dr_R_exp->SetLineWidth(2);
    gr_dr_R_gaus->SetLineWidth(2);
    gr_dr_R_levy->SetLineWidth(2);
    gr_dr_lambda_exp->SetLineWidth(2);
    gr_dr_lambda_gaus->SetLineWidth(2);
    gr_dr_lambda_levy->SetLineWidth(2);

    auto setup_frame = [&](const char* name, const char* yTitle, double yMin, double yMax) {
        TH1F* frame = new TH1F(name, "", n, 0.5, n + 0.5);
        frame->SetStats(0);
        frame->GetXaxis()->SetTitle("HFsumET bin");
        frame->GetYaxis()->SetTitle(yTitle);
        frame->GetYaxis()->SetRangeUser(yMin, yMax);
        for (int i = 0; i < n; i++) {
            frame->GetXaxis()->SetBinLabel(i + 1, bin_labels[i].c_str());
        }
        return frame;
    };

    auto draw_three_graphs = [&](TCanvas* canvas,
                                 TH1F* frame,
                                 const char* headerTitle,
                                 TGraphErrors* g1,
                                 TGraphErrors* g2,
                                 TGraphErrors* g3) {
        canvas->cd();
        frame->Draw();
        g1->Draw("LP SAME");
        g2->Draw("LP SAME");
        g3->Draw("LP SAME");

        TLegend* legend = new TLegend(0.70, 0.72, 0.90, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(g1, "Exponential", "lp");
        legend->AddEntry(g2, "Gaussian", "lp");
        legend->AddEntry(g3, "Levy", "lp");
        legend->Draw();
        drawCMSHeaders("#bf{CMS} #it{Work in Progress}", headerTitle, 3.0);
    };

    TH1F* frame_sr_R = setup_frame("frame_sr_R", "R [fm]", 0.0, 10.0);
    TH1F* frame_sr_lambda = setup_frame("frame_sr_lambda", "#lambda", 0.0, 1.2);
    TH1F* frame_dr_R = setup_frame("frame_dr_R", "R [fm]", 0.0, 10.0);
    TH1F* frame_dr_lambda = setup_frame("frame_dr_lambda", "#lambda", 0.0, 1.2);

    draw_three_graphs(c1_R, frame_sr_R, "PbPb 2.76 TeV | Single Ratio", gr_sr_R_exp, gr_sr_R_gaus, gr_sr_R_levy);
    draw_three_graphs(c1_lambda, frame_sr_lambda, "PbPb 2.76 TeV | Single Ratio", gr_sr_lambda_exp, gr_sr_lambda_gaus, gr_sr_lambda_levy);
    draw_three_graphs(c2_R, frame_dr_R, "PbPb 2.76 TeV | Double Ratio", gr_dr_R_exp, gr_dr_R_gaus, gr_dr_R_levy);
    draw_three_graphs(c2_lambda, frame_dr_lambda, "PbPb 2.76 TeV | Double Ratio", gr_dr_lambda_exp, gr_dr_lambda_gaus, gr_dr_lambda_levy);

    TCanvas* canvases[] = {c1_R, c1_lambda, c2_R, c2_lambda};
    const char* outputPath = "./imgs/test/resultPlots/";
    const char* outputPrefix = "resultPlots";

    save_canvas_images(canvases, 4, outputPath, outputPrefix, "png");
    save_canvas_images(canvases, 4, outputPath, outputPrefix, "pdf");

    std::cout << "Saved result plots to: " << outputPath << std::endl;

    AnalysisLog::instance().save("./logs", "resultPlots");
    std::cout << "Saved resultPlots log to: ./logs" << std::endl;
}

void resultPlotsLevyAlpha1Seed(double bins[], size_t numBins){
    ControlVar selectedControlVar = ControlVar::CENTHF;
    qMode modeLCMS = qMode::QLCMS;

    Double_t etaMin = 0.95;
    Double_t etaMax = -0.95;
    Double_t ptMin  = 0.50;

    std::vector<std::pair<FitFunctionType, FitInit>> models = {
        {FitFunctionType::DOUBLE_LEVY, {{0.6, 4.0, 1.0, 0.0, 1.0}}}
    };

    double plotXMin = 0.0;
    double plotXMax = 0.5;

    double plotYMin = 0.9;
    double plotYMax = 2.1;

    double fitMin = 0.04;
    double fitMax = 0.2;
    double fitMinBg = 0.2;

    Double_t q1 = 6.82;
    Double_t q2 = 8.4;

    std::vector<double> sr_R_values, sr_R_errors;
    std::vector<double> sr_lambda_values, sr_lambda_errors;
    std::vector<double> dr_R_values, dr_R_errors;
    std::vector<double> dr_lambda_values, dr_lambda_errors;
    std::vector<std::string> bin_labels;

    for (size_t i = 0; i + 1 < numBins; i++) {
        double bin_low = bins[i];
        double bin_high = bins[i + 1];

        std::string binLabel =
            std::to_string(static_cast<int>(bin_low)) + "-" +
            std::to_string(static_cast<int>(bin_high));

        MultipleFitResults fitResults = doubleRatioMixFitLevyAlphaOneSeed(
            selectedControlVar,
            models,
            bin_low, bin_high,
            etaMin, etaMax, ptMin,
            q1, q2,
            modeLCMS,
            fitMin, fitMax, fitMinBg,
            plotXMin, plotXMax, plotYMin, plotYMax
        );

        for (size_t j = 0; j < fitResults.results.size(); ++j) {
            const auto& res = fitResults.results[j];
            const auto& modelName = fitResults.modelNames[j];

            double lambda = res->Parameter(0);
            double R = res->Parameter(1);
            double lambdaErr = res->ParError(0);
            double RErr = res->ParError(1);

            if (modelName.find("sr_") == 0) {
                sr_R_values.push_back(R);
                sr_R_errors.push_back(RErr);
                sr_lambda_values.push_back(lambda);
                sr_lambda_errors.push_back(lambdaErr);
            } else if (modelName.find("dr_") == 0) {
                dr_R_values.push_back(R);
                dr_R_errors.push_back(RErr);
                dr_lambda_values.push_back(lambda);
                dr_lambda_errors.push_back(lambdaErr);
            }
        }

        bin_labels.push_back(binLabel);
    }

    int n = static_cast<int>(bin_labels.size());
    std::vector<double> x(n), xErr(n, 0.0);
    for (int i = 0; i < n; i++) {
        x[i] = i + 1.0;
    }

    TCanvas* c1_R = new TCanvas("c1_R_alpha1seed", "Single Ratio Levy Seed", 1200, 800);
    TCanvas* c1_lambda = new TCanvas("c1_lambda_alpha1seed", "Single Lambda Levy Seed", 1200, 800);
    TCanvas* c2_R = new TCanvas("c2_R_alpha1seed", "Double Ratio Levy Seed", 1200, 800);
    TCanvas* c2_lambda = new TCanvas("c2_lambda_alpha1seed", "Double Lambda Levy Seed", 1200, 800);
    c1_R->SetTopMargin(0.10);
    c1_lambda->SetTopMargin(0.10);
    c2_R->SetTopMargin(0.10);
    c2_lambda->SetTopMargin(0.10);
    c1_R->SetBottomMargin(0.12);
    c1_lambda->SetBottomMargin(0.12);
    c2_R->SetBottomMargin(0.12);
    c2_lambda->SetBottomMargin(0.12);

    TGraphErrors* gr_sr_R = new TGraphErrors(n, x.data(), sr_R_values.data(), xErr.data(), sr_R_errors.data());
    TGraphErrors* gr_sr_lambda = new TGraphErrors(n, x.data(), sr_lambda_values.data(), xErr.data(), sr_lambda_errors.data());
    TGraphErrors* gr_dr_R = new TGraphErrors(n, x.data(), dr_R_values.data(), xErr.data(), dr_R_errors.data());
    TGraphErrors* gr_dr_lambda = new TGraphErrors(n, x.data(), dr_lambda_values.data(), xErr.data(), dr_lambda_errors.data());

    auto style_graph = [](TGraphErrors* graph, int color, int markerStyle) {
        graph->SetMarkerStyle(markerStyle);
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        graph->SetLineWidth(2);
    };

    style_graph(gr_sr_R, kGreen + 2, 22);
    style_graph(gr_sr_lambda, kGreen + 2, 22);
    style_graph(gr_dr_R, kGreen + 2, 22);
    style_graph(gr_dr_lambda, kGreen + 2, 22);

    auto setup_frame = [&](const char* name, const char* yTitle, double yMin, double yMax) {
        TH1F* frame = new TH1F(name, "", n, 0.5, n + 0.5);
        frame->SetStats(0);
        frame->GetXaxis()->SetTitle("HFsumET bin");
        frame->GetYaxis()->SetTitle(yTitle);
        frame->GetYaxis()->SetRangeUser(yMin, yMax);
        for (int i = 0; i < n; i++) {
            frame->GetXaxis()->SetBinLabel(i + 1, bin_labels[i].c_str());
        }
        return frame;
    };

    auto draw_graph = [&](TCanvas* canvas,
                          TH1F* frame,
                          const char* headerTitle,
                          TGraphErrors* graph) {
        canvas->cd();
        frame->Draw();
        graph->Draw("LP SAME");

        TLegend* legend = new TLegend(0.68, 0.78, 0.90, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(graph, "Levy (#alpha=1 seed)", "lp");
        legend->Draw();
        drawCMSHeaders("#bf{CMS} #it{Work in Progress}", headerTitle, 3.0);
    };

    TH1F* frame_sr_R = setup_frame("frame_sr_R_alpha1seed", "R [fm]", 0.0, 10.0);
    TH1F* frame_sr_lambda = setup_frame("frame_sr_lambda_alpha1seed", "#lambda", 0.0, 1.2);
    TH1F* frame_dr_R = setup_frame("frame_dr_R_alpha1seed", "R [fm]", 0.0, 10.0);
    TH1F* frame_dr_lambda = setup_frame("frame_dr_lambda_alpha1seed", "#lambda", 0.0, 1.2);

    draw_graph(c1_R, frame_sr_R, "PbPb 2.76 TeV | Single Ratio | Levy (#alpha=1 seed)", gr_sr_R);
    draw_graph(c1_lambda, frame_sr_lambda, "PbPb 2.76 TeV | Single Ratio | Levy (#alpha=1 seed)", gr_sr_lambda);
    draw_graph(c2_R, frame_dr_R, "PbPb 2.76 TeV | Double Ratio | Levy (#alpha=1 seed)", gr_dr_R);
    draw_graph(c2_lambda, frame_dr_lambda, "PbPb 2.76 TeV | Double Ratio | Levy (#alpha=1 seed)", gr_dr_lambda);

    TCanvas* canvases[] = {c1_R, c1_lambda, c2_R, c2_lambda};
    const char* outputPath = "./imgs/test/resultPlots/";
    const char* outputPrefix = "resultPlotsLevyAlpha1Seed";

    save_canvas_images(canvases, 4, outputPath, outputPrefix, "png");
    save_canvas_images(canvases, 4, outputPath, outputPrefix, "pdf");

    std::cout << "Saved result plots to: " << outputPath << std::endl;

    AnalysisLog::instance().save("./logs", "resultPlotsLevyAlpha1Seed");
    std::cout << "Saved resultPlotsLevyAlpha1Seed log to: ./logs" << std::endl;
}

int main(){
    // to compile use:
    // g++ -std=c++17 -pthread resultPlots.cpp -o resultPlots `root-config --cflags --libs`
    // to run use:
    // ./resultPlots
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    double bins[] = {3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0};
    //resultPlots(bins, sizeof(bins)/sizeof(bins[0]));
    resultPlotsLevyAlpha1Seed(bins, sizeof(bins)/sizeof(bins[0]));
    return 0;
}
