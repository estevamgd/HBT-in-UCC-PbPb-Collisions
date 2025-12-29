#ifndef NORMALIZER_H
#define NORMALIZER_H

#include "TH1D.h"
#include "TH2D.h"
#include "normalizer.h"

void normalizer(TH1D* histograms[],
                int numHistograms,
                Double_t q1,
                Double_t q2,
                Double_t scale = 1.0)
{
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

        if (integral > 0.0) {
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