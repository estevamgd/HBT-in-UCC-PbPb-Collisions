#ifndef DATA_FUNC_H
#define DATA_FUNC_H

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "TMath.h"
#include <TF1.h>

using namespace std;

#define HBARC 0.1973269788       // hbar*c in GeV*fm 
#define PI_MASS 0.13957018       // Pion mass in GeV/c^2
#define ALPHAEM 0.0072973525664  // 1.0/137.035999139, fine structure constant
#define PREFACTOR 0.016215165    // ALPHAEM * PI * PI_MASS / HBARC

enum class FitFunctionType { EXPONENTIAL = 0, GAUSSIAN = 1, LEVY = 2, LEVY2 = 3, BACKGROUND = 4, DOUBLE_LEVY = 5 };

// to reject a range in the fit -- in principle did not reject any range
Double_t reject_range_min = 0.01;
Double_t reject_range_max = 0.00001;

// Foward declaration to avoid bugs
Double_t KGamow  (Double_t *x, Double_t *par);
Double_t KMod    (Double_t *x, Double_t *par);
Double_t KCoulomb(Double_t *x, Double_t *par);
Double_t FitExp    (Double_t*, Double_t*);
Double_t FitGauss  (Double_t*, Double_t*);
Double_t FitLevy   (Double_t*, Double_t*);
Double_t FitLevy2  (Double_t*, Double_t*);
Double_t FitBG     (Double_t*, Double_t*);
Double_t FitLevyDR (Double_t*, Double_t*);

// ===== Fit Factory ===== //
struct FitInit {
    std::vector<double> values;
};

struct FitParamInfo {
    int index;            // TF1 parameter index
    std::string label;    // "#lambda", "R", "#alpha"
    std::string unit;     // " fm", "" 
};

struct FitModelConfig {
    const char* name;
    const char* displayName;
    Double_t (*func)(Double_t*, Double_t*);
    int nPar;
    std::vector<std::string> parNames;
    std::vector<std::pair<double,double>> parLimits;
    std::vector<FitParamInfo> legendParams;
};

FitModelConfig getFitModelConfig(FitFunctionType type);

TF1* fitHistogram(
    TH1D* hist,
    FitFunctionType type,
    const FitInit& init,
    double fitMin,
    double fitMax,
    TFitResultPtr* fitResult = nullptr
);

// ===== Functions Definitions ===== //

// Exponential function + long range from: https://github.com/i5albg/hbt_analysis/blob/main/final_HBT.C#L44-L67
Double_t FitExp(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-par[2]*x[0]/HBARC))*(1+par[3]*x[0]);}
    return v;
}

// Gaussian function + long range from: https://github.com/i5albg/hbt_analysis/blob/main/final_HBT.C#L44-L67
Double_t FitGauss(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/HBARC,2.)))*(1+par[3]*x[0]);}
    return v;
} 

// Levy function  from: https://github.com/i5albg/hbt_analysis/blob/main/final_HBT.C#L44-L67
Double_t FitLevy(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/HBARC,par[4])))*(1+par[3]*x[0]);}
    return v;
}

// Levy function
Double_t FitLevy2(Double_t* x, Double_t* par){
    Double_t v = 0;
    Double_t par2[2];
    par2[0] = par[1]; 
    par2[1] = par[2]; 

    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= 1 - par[0] + par[0]*(1 + exp(-pow(par[1]*x[0]/HBARC,par[2])))*KCoulomb(x, par2);}
    return v;
}

// Background function 
Double_t FitBG(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/HBARC,2)))*(1 - par[3]*exp(-pow(par[4]*x[0]/HBARC,2)));}
    return v;
}

// Double Ratio function 
Double_t FitLevyDR(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[4]*(1 - par[0] + par[0]*(1 + exp(-pow(par[1]*x[0]/HBARC,par[2])))*KCoulomb(x, par))*(1+par[3]*x[0]);}
    return v;
}


// Coulomb Correction
Double_t KCoulomb(Double_t *x, Double_t *par){
    Double_t aA = 0.26984, aB = -0.49123, aC = 0.03523, aD = -1.31628, 
    aE = 0.00359, bA = 2.37267, bB = 0.58631, bC = 2.24867, bD = -1.43278, 
    bE = -0.05216, bF = 0.72943, cA = -4.30347, cB = 0.00001, cC = 3.30346, 
    cD = 0.000001, cE = 0.000003, cF = 1.68883, dA = 0.00057, dB = -0.80527, 
    dC = -0.19261, dD = 2.77504, dE = 2.02951, dF = 1.07906;

    Double_t A = pow((aA * par[1] + aB),2) + pow((aC * par[1] + aD),2) + aE * pow((par[1] * par[0] + 1),2);
    Double_t B = (1. + bA*pow(par[0],bB) -pow(par[1],bC)) / (par[1]*par[1]*par[0]*(pow(par[1],bD)+bE*pow(par[0],bF)));
    Double_t C = ((cA  + pow(par[1],cB) + cC*pow(par[0],cD))/cE)*(pow(par[1]/par[0],cF));
    Double_t D = (dA + (pow(par[0],dB) + dC*pow(par[1],dF))/(pow(par[0],dD)*pow(par[1],dE)));
    Double_t par2[6] = {par[0], par[1], A, B, C, D};

    Double_t Aa = 0.126253, Ab = 0.05385, Ac = -0.00913, Ad = -0.01846,
    Ae = 0.00085, Af = 0.00042, Ba = 19.31620, Bb = 5.58961, Bc = 2.26264, 
    Bd = -1.28486, Be = -0.08216, Bf = 0.02384;

    Double_t A2 = Aa + Ab*par[1] + Ac*par[0] + Ad*par[1]*par[0] + Ae*par[0]*par[0] + Af*par[1]*par[1]*par[0]*par[0];
    Double_t B2 = Ba + Bb*par[1] + Bc*par[0] + Bd*par[1]*par[0] + Be*par[0]*par[0] + Bf*par[1]*par[1]*par[0]*par[0];

    Double_t q0 = 0.07;
    Double_t n = 20.0;

    Double_t E = 1. + A2 * exp( - B2 * x[0] ); // exponential smoothing
    Double_t F = 1. / (1. + pow(x[0]/q0,n)); // cutoff factor

    Double_t gamow = KGamow(x, par);
    Double_t mod = KMod(x, par2);
    return 1./(F*(1./gamow)*(1./mod) + (1.-F)*E); 
}

// Modified Coulomb Correction
Double_t KMod(Double_t *x, Double_t *par2){  // Direct implementation from paper
    return 1. + (par2[2]*((PREFACTOR*par2[0])/(par2[1])))/(1. + par2[3]*x[0]*par2[0]/(par2[1]*HBARC) + par2[4]*pow(x[0]*par2[0]/(par2[1]*HBARC),2) + par2[5]*pow(x[0]*par2[0]/(par2[1]*HBARC),4));
}

//  Gamow Factor from: https://github.com/csanadm/coulcorrlevyparam/blob/0c3d63e765762ec06d0eb66e42b95b8561e6a669/coulcorr_param.cc#L84-L88
Double_t KGamow(Double_t *x, Double_t *par){
    double eta = ALPHAEM*PI_MASS/(x[0]); //x[0] is Q = 2k = |p1-p2|
    return ( M_PI*eta*2.0/(exp(2.0*M_PI*eta)-1) );
}

// Function to calculate qinv from two 4-momenta
Double_t GetQ(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, 
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2) {
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Sum4V = p1+p2;
    Double_t q = Sum4V.M2() - 4*PI_MASS*PI_MASS;
    return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

// QLong^2 calculation
Double_t GetQLongSquared(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, 
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2){
    Double_t px1 = p1.Px(), py1 = p1.Py(), pz1 = p1.Pz(), e1  = p1.E(); 
    Double_t px2 = p2.Px(), py2 = p2.Py(), pz2 = p2.Pz(), e2  = p2.E();

    return 4*pow(pz1*e2 - pz2*e1, 2)/(pow((e1 + e2), 2) - pow((pz1 + pz2), 2));
}

// QLCMS calculation
Double_t GetQLCMS(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, 
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2){
    Double_t px1 = p1.Px(), py1 = p1.Py(), pz1 = p1.Pz(), e1  = p1.E(); 
    Double_t px2 = p2.Px(), py2 = p2.Py(), pz2 = p2.Pz(), e2  = p2.E();
    Double_t qlongsqrd = GetQLongSquared(p1, p2);    

    return sqrt(pow(px1 - px2, 2) + pow(py1 - py2, 2) + qlongsqrd);
}

#endif