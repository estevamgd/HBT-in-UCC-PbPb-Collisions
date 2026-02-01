#ifndef NORMALIZER_H
#define NORMALIZER_H

#include "TH1D.h"
#include "TH2D.h"
#include "normalizer.h"

std::vector<TH1D> normalizer(TH1D* histograms[],
                             int numHistograms,
                             Double_t q1,
                             Double_t q2,
                             Double_t scale = 1.0)
{
    std::vector<TH1D> result;
    result.reserve(numHistograms);

    for (int i = 0; i < numHistograms; i++) {

        TH1D* hist = histograms[i];
        if (!hist) continue;

        // Convert q-range to bin indices
        int bin1 = hist->GetXaxis()->FindBin(q1);
        int bin2 = hist->GetXaxis()->FindBin(q2);

        // Protect against out-of-range values
        bin1 = std::max(bin1, 1);
        bin2 = std::min(bin2, hist->GetNbinsX());

        // Integral in [q1, q2], including bin widths
        double integral = hist->Integral(bin1, bin2, "width");
        std::cout << "integral: " << integral << std::endl;
        
        if (integral > 0.0) {
            hist->Scale(scale / integral);
        }

        result.push_back(*hist);
    }

    return result;
}


void normalizer(TH2D* histograms[],
                int numHistograms,
                Double_t q1x, Double_t q2x,
                Double_t q1y, Double_t q2y,
                Double_t scale = 1.0)
{
    for (int i = 0; i < numHistograms; i++) {

        TH2D* hist = histograms[i];
        if (!hist) continue;

        // Convert q-range to bin indices
        int bin1x = hist->GetXaxis()->FindBin(q1x);
        int bin2x = hist->GetXaxis()->FindBin(q2x);
        int bin1y = hist->GetYaxis()->FindBin(q1y);
        int bin2y = hist->GetYaxis()->FindBin(q2y);

        // Protect against out-of-range values
        bin1x = std::max(bin1x, 1);
        bin2x = std::min(bin2x, hist->GetNbinsX());
        
        bin1y = std::max(bin1y, 1);
        bin2y = std::min(bin2y, hist->GetNbinsY());

        // Integral in [q1, q2], including bin widths
        double integral = hist->Integral(bin1x, bin2x, 
            bin1y, bin2y, "width");
        
        if (integral != 0.0) {
            hist->Scale(scale / integral);
        }
    }
}

void normalizerExcludeDEta(TH2D* histograms[],
                           int numHistograms,
                           Double_t detaCut = 0.04,
                           Double_t scale = 1.0)
{
    for (int i = 0; i < numHistograms; i++) {

        TH2D* hist = histograms[i];
        if (!hist) continue;

        TAxis* ax = hist->GetXaxis();
        TAxis* ay = hist->GetYaxis();

        int binEtaMin = 1;
        int binEtaMax = ax->GetNbins();

        int binPhiMin = 1;
        int binPhiMax = ay->GetNbins();

        // bins corresponding to ±detaCut
        int binCutHigh = ax->FindBin(+detaCut);

        // region Δη > +cut
        double integralRight = hist->Integral(
            binCutHigh + 1, binEtaMax,
            binPhiMin, binPhiMax
        );

        double integral = integralRight;
        std::cout << "integral: " << integral << std::endl;

        if (integral != 0.0) {
            hist->Scale(scale / integral);
        }
    }
}

void normalizer(TH1D* histograms[], int numHistograms, Double_t scale = 1.0) {
    for (int i = 0; i < numHistograms; i++) {
        TH1D* hist = histograms[i];
        if (hist->Integral() != 0) {  // Avoid division by zero
            hist->Scale(scale / hist->Integral(), "width");
        }
    }
}

void normalizer2d(TH2D* histograms[], int numHistograms, Double_t scale = 1.0) {
    for (int i = 0; i < numHistograms; i++) {
        TH2D* hist = histograms[i];
        if (hist->Integral() != 0) {  // Avoid division by zero
            hist->Scale(scale / hist->Integral(), "width");
        }
    }
}

#endif