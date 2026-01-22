#ifndef MY_FUNC_H
#define MY_FUNC_H

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

using FourVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;

enum class ControlVar { CENT = 0, MULT = 1, CENTHF = 2 };

enum class qMode { QINV = 0, QLCMS = 1 };

enum LoopMode { SINGLE = 0, BOTH = 1, DOUBLE = 2 }; // Define an enum for modes

static const Int_t colors[9] = {
    kRed-7,
    kYellow-7,
    kGreen-7,
    kBlue-7,
    kPink+10,
    kOrange+10,
    kSpring+10,
    kCyan-7,
    kMagenta-4
};

static const Int_t lineStyles[4] = {
    1, // solid
    2, // dashed
    3, // dotted
    4  // dash-dotted
};

int setLoopMode(LoopMode mode) {
    int xmode = 2;

    switch (mode) {
        case SINGLE:
            xmode = 0;
            break;
        case BOTH:
            xmode = 1;
            break;
        case DOUBLE:
            xmode = 2;
            break;
        default:
            std::cout << "Choose between: SINGLE, BOTH, or DOUBLE." << std::endl;
    }
    return xmode;
}

#endif 