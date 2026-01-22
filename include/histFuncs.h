#ifndef HISTFUNCS_H
#define HISTFUNCS_H

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <TStyle.h>
#include <fstream>
#include <ctime>
#include <sstream>
#include <filesystem>
#include "my_func.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "TMath.h"
#include <vector>
#include <thread>
#include <mutex>
#include <cmath>
#include <chrono>
#include <string>
#include "data_func.h"
#include <TSystem.h>
#include <TLatex.h>
#include "getFuncs.h"
#include "my_func.h"


/* Histogram layout

Colors: 
kWhite  = 0, kBlack     = 1, kGray      = 920, kRed     = 632, kGreen   = 416, 
kBlue   = 600, kYellow  = 400, kMagenta = 616, kCyan    = 432, kOrange  = 800, 
kSpring = 820, kTeal    = 840, kAzure   = 860, kViolet  = 880, kPink    = 900

Fill Style:
kFDotted1  = 3001, kFDotted2    = 3002, kFDotted3  = 3003,
kFHatched1 = 3004, kHatched2    = 3005, kFHatched3 = 3006,
kFHatched4 = 3007, kFWicker     = 3008, kFScales   = 3009,
kFBricks   = 3010, kFSnowflakes = 3011, kFCircles  = 3012,
kFTiles    = 3013, kFMondrian   = 3014, kFDiamonds = 3015,
kFWaves1   = 3016, kFDashed1    = 3017, kFDashed2  = 3018,
kFAlhambra = 3019, kFWaves2     = 3020, kFStars1   = 3021,
kFStars2   = 3022, kFPyramids   = 3023, kFFrieze   = 3024,
kFMetopes  = 3025, kFEmpty      = 0   , kFSolid    = 1

Line Style:
kSolid = 1, kDashed = 2, kDotted = 3, kDashDotted = 4
*/
TH1D* cHist(const char* name, const char* xAxisTitle, const char* yAxisTitle, 
                        double xMin, double xMax, 
                        double ninterval = 275., double nlength = 0.02, double nscale = 1., 
                        int fill_color = 920, double fill_alpha = 1, int fill_style = 1001, 
                        int line_color = 1, double line_alpha = 1, int line_style = 1001, 
                        int line_width = 1, double label_size = 0.04) {
    // Create the histogram
    double scale = numBins(ninterval, nlength, nscale);
    TH1D* h = new TH1D(name, "", scale, xMin, xMax);

    // Set visual properties
    h->SetFillColorAlpha(fill_color, fill_alpha); 
    h->SetFillStyle(fill_style); 
    h->SetLineColorAlpha(line_color, line_alpha);
    h->SetLineStyle(line_style);
    h->SetLineWidth(line_width); 

    // Set titles and axis labels
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xAxisTitle);
    h->GetYaxis()->SetTitle(yAxisTitle);
    
    // Set label sizes
    h->GetXaxis()->SetLabelSize(label_size);
    h->GetYaxis()->SetLabelSize(label_size);
    
    return h;
}

TH1F* cHistFloat(const char* name, const char* xAxisTitle, const char* yAxisTitle, 
    double xMin, double xMax, 
    double ninterval = 275., double nlength = 0.02, double nscale = 1., 
    int fill_color = 920, double fill_alpha = 1, int fill_style = 1001, 
    int line_color = 1, double line_alpha = 1, int line_style = 1001, 
    int line_width = 1, double label_size = 0.04) {
// Create the histogram
double scale = numBins(ninterval, nlength, nscale);
TH1F* h = new TH1F(name, "", scale, xMin, xMax);

// Set visual properties
h->SetFillColorAlpha(fill_color, fill_alpha); 
h->SetFillStyle(fill_style); 
h->SetLineColorAlpha(line_color, line_alpha);
h->SetLineStyle(line_style);
h->SetLineWidth(line_width); 

// Set titles and axis labels
h->SetTitle("");
h->GetXaxis()->SetTitle(xAxisTitle);
h->GetYaxis()->SetTitle(yAxisTitle);

// Set label sizes
h->GetXaxis()->SetLabelSize(label_size);
h->GetYaxis()->SetLabelSize(label_size);

return h;
}

void drawCMSHeaders(
    const char* leftText  = "#bf{CMS} #it{Work in Progress}",
    const char* rightText = "",
    double rightTextSizeFactor = 1.0
) {
    if (!gPad) return;

    // ---- HARD RESET OF PAD STATE ----
    gPad->SetTitle("");
    gPad->SetName("");
    gPad->SetGrid(0, 0);
    gPad->SetTicks(1, 1);
    gPad->Update();

    const double yPos = 1.0 - gPad->GetTopMargin() + 0.01;

    // ===== Left CMS label =====
    TLatex latexLeft;
    latexLeft.SetNDC();
    latexLeft.SetTextFont(42);
    latexLeft.SetTextSize(0.035);
    latexLeft.SetTextAlign(11);

    latexLeft.DrawLatex(
        gPad->GetLeftMargin(),
        yPos,
        leftText
    );

    // ===== Right header =====
    if (rightText && rightText[0] != '\0') {

        TLatex latexRight;
        latexRight.SetNDC();
        latexRight.SetTextFont(42);
        latexRight.SetTextAlign(31);

        const double xRight = 1.0 - gPad->GetRightMargin();
        const double xLeft  = gPad->GetLeftMargin() + 0.25;
        const double maxWidth = xRight - xLeft;

        double textSize = 0.04 * rightTextSizeFactor;
        latexRight.SetTextSize(textSize);

        latexRight.SetTitle(rightText);
        double textWidth = latexRight.GetXsize();

        while (textWidth > maxWidth && textSize > 0.015) {
            textSize *= 0.95;
            latexRight.SetTextSize(textSize);
            latexRight.SetTitle(rightText);
            textWidth = latexRight.GetXsize();
        }

        latexRight.DrawLatex(xRight, yPos, rightText);
    }

    // ---- FINAL FLUSH ----
    gPad->Modified();
    gPad->Update();
}



void no_statbox(TH1D *histograms[], int numHistograms) {
    for (int i = 0; i < numHistograms; i++) {
        histograms[i]->SetStats(0);
    }
}


#endif