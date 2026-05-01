#include "../include/ratiosAndFits.h"
#include "../include/statistics.h"

#include "TMatrixDSym.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct StatisticsFitSummary {
    std::string binLabel;
    std::string ratioName;
    std::string fitStage;
    double lambda;
    double lambdaErr;
    double radius;
    double radiusErr;
    double alpha;
    double alphaErr;
    double chi2;
    int ndf;
    double pValue;
};

void writeStatisticsSummary(const std::vector<StatisticsFitSummary>& rows,
                            const char* outputPath)
{
    std::ofstream out(outputPath);
    if (!out.is_open()) {
        std::cerr << "Could not open statistics summary: " << outputPath << std::endl;
        return;
    }

    out << "HFsumET_bin,ratio,fit_stage,lambda,lambda_err,R,R_err,alpha,alpha_err,chi2,ndf,p_value\n";
    out << std::fixed << std::setprecision(8);

    for (const auto& row : rows) {
        out << row.binLabel << ","
            << row.ratioName << ","
            << row.fitStage << ","
            << row.lambda << ","
            << row.lambdaErr << ","
            << row.radius << ","
            << row.radiusErr << ","
            << row.alpha << ","
            << row.alphaErr << ","
            << row.chi2 << ","
            << row.ndf << ","
            << row.pValue << "\n";
    }
}

void writeFTestSummary(const std::vector<FTestSummary>& rows,
                       const char* outputPath)
{
    std::ofstream out(outputPath);
    if (!out.is_open()) {
        std::cerr << "Could not open F-test summary: " << outputPath << std::endl;
        return;
    }

    out << "HFsumET_bin,ratio,chi2_reduced,ndf_reduced,chi2_full,ndf_full,"
        << "numerator_dof,denominator_dof,MSR,MSE,F_statistic,p_value\n";
    out << std::fixed << std::setprecision(10);

    for (const auto& row : rows) {
        out << row.binLabel << ","
            << row.ratioName << ","
            << row.chi2Reduced << ","
            << row.ndfReduced << ","
            << row.chi2Full << ","
            << row.ndfFull << ","
            << row.numeratorDof << ","
            << row.denominatorDof << ","
            << row.msr << ","
            << row.mse << ","
            << row.fStatistic << ","
            << row.pValue << "\n";
    }
}

std::string sanitizeForFileName(const std::string& text)
{
    std::string sanitized;
    for (char c : text) {
        if ((c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9')) {
            sanitized += c;
        } else {
            sanitized += '_';
        }
    }
    return sanitized;
}

void printAndSaveCorrelationMatrix(const std::string& binLabel,
                                   const std::string& ratioName,
                                   const std::string& fitStage,
                                   const FitOutput& fit,
                                   const char* outputDir)
{
    if (!fit.function || !fit.result.Get()) {
        return;
    }

    const int nPar = fit.function->GetNpar();
    TMatrixDSym corr = fit.result->GetCorrelationMatrix();

    std::cout << "\nCorrelation matrix | " << binLabel
              << " | " << ratioName
              << " | " << fitStage << std::endl;

    for (int i = 0; i < nPar; ++i) {
        std::cout << "  [" << i << "] " << fit.function->GetParName(i) << std::endl;
    }
    corr.Print();

    std::ostringstream path;
    path << outputDir << "/corr_"
         << sanitizeForFileName(binLabel) << "_"
         << sanitizeForFileName(ratioName) << "_"
         << sanitizeForFileName(fitStage) << ".csv";

    std::ofstream out(path.str());
    if (!out.is_open()) {
        std::cerr << "Could not save correlation matrix: " << path.str() << std::endl;
        return;
    }

    out << "parameter";
    for (int i = 0; i < nPar; ++i) {
        out << "," << fit.function->GetParName(i);
    }
    out << "\n";

    out << std::fixed << std::setprecision(10);
    for (int i = 0; i < nPar; ++i) {
        out << fit.function->GetParName(i);
        for (int j = 0; j < nPar; ++j) {
            out << "," << corr(i, j);
        }
        out << "\n";
    }
}

void appendFitSummary(std::vector<StatisticsFitSummary>& summaries,
                      const std::string& binLabel,
                      const std::string& ratioName,
                      const std::string& fitStage,
                      const FitOutput& fit,
                      const char* matrixOutputDir)
{
    if (!fit.function || !fit.result.Get()) {
        std::cerr << "Missing " << fitStage << " result for " << ratioName
                  << " in bin " << binLabel << std::endl;
        return;
    }

    StatisticsFitSummary summary;
    summary.binLabel = binLabel;
    summary.ratioName = ratioName;
    summary.fitStage = fitStage;
    summary.lambda = fit.function->GetParameter(0);
    summary.lambdaErr = fit.function->GetParError(0);
    summary.radius = fit.function->GetParameter(1);
    summary.radiusErr = fit.function->GetParError(1);
    summary.alpha = fit.function->GetParameter(2);
    summary.alphaErr = fit.function->GetParError(2);
    summary.chi2 = fit.result->Chi2();
    summary.ndf = fit.result->Ndf();
    summary.pValue = fit.result->Prob();

    std::cout << ratioName << " | " << fitStage
              << " | lambda = " << summary.lambda << " +/- " << summary.lambdaErr
              << ", R = " << summary.radius << " +/- " << summary.radiusErr
              << ", alpha = " << summary.alpha << " +/- " << summary.alphaErr
              << ", chi2/NDF = " << summary.chi2 << "/" << summary.ndf
              << ", p = " << summary.pValue << std::endl;

    printAndSaveCorrelationMatrix(
        binLabel,
        ratioName,
        fitStage,
        fit,
        matrixOutputDir
    );

    summaries.push_back(summary);
}

void runStatisticsForHFsumETBins(
    const std::vector<double>& bins,
    std::vector<StatisticsFitSummary>& summaries,
    std::vector<FTestSummary>& fTestSummaries)
{
    ControlVar selectedControlVar = ControlVar::CENTHF;
    qMode mode = qMode::QLCMS;
    const char* qmodename = (mode == qMode::QINV) ? "qinv" : "qlcms";
    const char* selectionVarName = getSelVarName(selectedControlVar);

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

    const char* matrixOutputDir = "./data/statistics/correlation_matrices";
    gSystem->mkdir(matrixOutputDir, true);

    for (size_t i = 0; i + 1 < bins.size(); ++i) {
        double binLow = bins[i];
        double binHigh = bins[i + 1];
        std::string binLabel =
            std::to_string(static_cast<int>(binLow)) + "-" +
            std::to_string(static_cast<int>(binHigh));

        std::cout << "\n=== HFsumET " << binLabel << " ===" << std::endl;

        TString searchPattern = TString::Format(
            "./data/signal_mix/%s_eta-abs%.2f_CENTHF_%d-%d_cent%dto%d*.root",
            qmodename,
            etaMin,
            (int)binLow,
            (int)binHigh,
            (int)binLow,
            (int)binHigh
        );

        TH1D* q_sigHist_cor = getHistogram(searchPattern.Data(), "hSigSS");
        TH1D* q_mixHist_cor = getHistogram(searchPattern.Data(), "hMixSS");

        TH1D* singleRatio = histhistRatioWithComp(
            q_sigHist_cor, q_mixHist_cor, q1, q2, "sr_cor", "testMix"
        );

        saveRatio(
            singleRatio,
            Form("sr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, binLow, binHigh),
            plotXMin, plotXMax, plotYMin, plotYMax,
            "sr_cor", "PbPb 2.76 TeV | Single Ratio",
            "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix",
            "Single Ratio (Sig/Mix)"
        );

        FitInit levyInit = findModelInitOrDefault(
            models,
            FitFunctionType::DOUBLE_LEVY,
            {{0.6, 4.0, 1.0, 0.0, 1.0}}
        );

        AlphaSeedFitOutput singleRatioLevyFit = fitDoubleLevyAlphaOneThenFreeDetailed(
            singleRatio,
            levyInit,
            fitMin,
            fitMax
        );

        appendFitSummary(
            summaries,
            binLabel,
            "sr_Levy",
            "alpha fixed to 1",
            singleRatioLevyFit.fixedAlphaFit,
            matrixOutputDir
        );
        appendFitSummary(
            summaries,
            binLabel,
            "sr_Levy",
            "alpha free initialized from fixed-alpha fit",
            singleRatioLevyFit.freeAlphaFit,
            matrixOutputDir
        );

        fTestSummaries.push_back(fTestAlphaFixedVsFree(
            binLabel,
            "sr_Levy",
            singleRatioLevyFit.fixedAlphaFit,
            singleRatioLevyFit.freeAlphaFit
        ));

        saveFits(
            singleRatio,
            {singleRatioLevyFit.freeAlphaFit},
            "c_single_ratio_levy_alpha1_seed_fit",
            Form("PbPb 2.76 TeV | Single Ratio | %s: %.0f-%.0f", selectionVarName, binLow, binHigh),
            Form("fit_alpha1seed_sr_%s_%s_%.0f-%.0f", qmodename, selectionVarName, binLow, binHigh),
            "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Sig/Mix",
            mode, plotXMin, plotXMax, plotYMin, plotYMax
        );

        FitInit bgInit{{1.0, 0.1, 2.0, 0.1, 2.0}};
        TFitResultPtr bgRes;
        TF1* bgFit = fitHistogram(
            singleRatio,
            &bgRes,
            FitFunctionType::BACKGROUND,
            bgInit,
            fitMinBg
        );

        TH1D* doubleRatio = histfuncRatio(singleRatio, bgFit, "dr_cor");

        saveRatio(
            doubleRatio,
            Form("dr_cor_%s_%s_%.0f-%.0f", qmodename, selectionVarName, binLow, binHigh),
            plotXMin, plotXMax, plotYMin, plotYMax,
            "dr_cor", "PbPb 2.76 TeV | Double Ratio",
            "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit",
            "Double Ratio (Data/Fit)"
        );

        AlphaSeedFitOutput doubleRatioLevyFit = fitDoubleLevyAlphaOneThenFreeDetailed(
            doubleRatio,
            levyInit,
            fitMin,
            fitMax
        );

        appendFitSummary(
            summaries,
            binLabel,
            "dr_Levy",
            "alpha fixed to 1",
            doubleRatioLevyFit.fixedAlphaFit,
            matrixOutputDir
        );
        appendFitSummary(
            summaries,
            binLabel,
            "dr_Levy",
            "alpha free initialized from fixed-alpha fit",
            doubleRatioLevyFit.freeAlphaFit,
            matrixOutputDir
        );

        fTestSummaries.push_back(fTestAlphaFixedVsFree(
            binLabel,
            "dr_Levy",
            doubleRatioLevyFit.fixedAlphaFit,
            doubleRatioLevyFit.freeAlphaFit
        ));

        saveFits(
            doubleRatio,
            {doubleRatioLevyFit.freeAlphaFit},
            "c_double_ratio_levy_alpha1_seed_fit",
            Form("PbPb 2.76 TeV | Double Ratio | %s: %.0f-%.0f", selectionVarName, binLow, binHigh),
            Form("fit_alpha1seed_dr_%s_%s_%.0f-%.0f", qmodename, selectionVarName, binLow, binHigh),
            "; q_{LCMS} [GeV]; C_{2}(q_{LCMS}) = Data/Fit",
            mode, plotXMin, plotXMax, plotYMin, plotYMax
        );

        delete singleRatioLevyFit.fixedAlphaFit.function;
        delete singleRatioLevyFit.freeAlphaFit.function;
        delete doubleRatioLevyFit.fixedAlphaFit.function;
        delete doubleRatioLevyFit.freeAlphaFit.function;
        delete bgFit;
        delete q_sigHist_cor;
        delete q_mixHist_cor;
        delete singleRatio;
        delete doubleRatio;
    }
}


int main()
{
    // Compile from build/ with:
    // g++ -std=c++17 -pthread statistics.cpp -o statistics `root-config --cflags --libs`
    //
    // Run from build/ with:
    // ./statistics
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    const std::vector<double> hfSumEtBins = {
        3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0
    };

    std::vector<StatisticsFitSummary> summaries;
    std::vector<FTestSummary> fTestSummaries;
    runStatisticsForHFsumETBins(hfSumEtBins, summaries, fTestSummaries);

    gSystem->mkdir("./data/statistics", true);
    const char* fitSummaryPath = "./data/statistics/hfsumet_alpha1_seed_summary.csv";
    const char* fTestSummaryPath = "./data/statistics/hfsumet_f_test_summary.csv";
    writeStatisticsSummary(summaries, fitSummaryPath);
    writeFTestSummary(fTestSummaries, fTestSummaryPath);

    AnalysisLog::instance().save("./logs", "statistics_alpha1_seed");
    plot_ftest();

    std::cout << "\nSaved statistics summary to: " << fitSummaryPath << std::endl;
    std::cout << "Saved F-test summary to: " << fTestSummaryPath << std::endl;
    std::cout << "Saved correlation matrices to: ./data/statistics/correlation_matrices" << std::endl;
    std::cout << "Saved log to: ./logs/statistics_alpha1_seed_log-*.txt" << std::endl;

    return 0;
}
