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

// to reject a range in the fit -- in principle did not reject any range
Double_t reject_range_min = 0.01;
Double_t reject_range_max = 0.00001;

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
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= 1 - par[0] + par[0]*(1 + exp(-pow(par[1]*x[0]/HBARC,par[2])))*KCoulomb(x, par);}
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

//  Gamow Factor from: https://github.com/csanadm/coulcorrlevyparam/blob/0c3d63e765762ec06d0eb66e42b95b8561e6a669/coulcorr_param.cc#L84-L88
Double_t KGamow(Double_t *x, Double_t *par){
    double eta = ALPHAEM*PI_MASS/(x[0]); //x[0] is Q = 2k = |p1-p2|
    return ( M_PI*eta*2.0/(exp(2.0*M_PI*eta)-1) );
}

// Coulomb Correction
Double_t KCoulomb(Double_t *x, Double_t *par){
    Double_t aA = 0.26984, aB = -0.49123, aC = 0.03523, aD = -1.31628, 
    aE = 0.00359, bA = 2.37267, bB = 0.58631, bC = 2.24867, bD = -1.43278, 
    bE = -0.05216, bF = 0.72943, cA = -4.30347, cB = 0.00001, cC = 3.30346, 
    cD = 0.000001, cE = 0.000003, cF = 1.68883, dA = 0.00057, dB = -0.80527, 
    dC = -0.19261, dD = 2.77504, dE = 2.02951, dF = 1.07906;

    Double_t A = pow((aA * par[2] + aB),2) + pow((aC * par[1] + aD),2) + aE * pow((par[2] * par[1] + 1),2);
    Double_t B = (1. + bA*pow(par[1],bB) -pow(par[2],bC)) / (par[2]*par[2]*par[1]*(pow(par[2],bD)+bE*pow(par[1],bF)));
    Double_t C = ((cA  + pow(par[2],cB) + cC*pow(par[1],cD))/cE)*(pow(par[2]/par[1],cF));
    Double_t D = (dA + (pow(par[1],dB) + dC*pow(par[2],dF))/(pow(par[1],dD)*pow(par[2],dE)));
    Double_t par2[4] = {A, B, C, D};

    Double_t Aa = 0.126253, Ab = 0.05385, Ac = -0.00913, Ad = -0.01846,
    Ae = 0.00085, Af = 0.00042, Ba = 19.31620, Bb = 5.58961, Bc = 2.26264, 
    Bd = -1.28486, Be = -0.08216, Bf = 0.02384;

    Double_t A2 = Aa + Ab*par[2] + Ac*par[1] + Ad*par[2]*par[1] + Ae*par[1]*par[1] + Af*par[2]*par[2]*par[1]*par[1];
    Double_t B2 = Ba + Bb*par[2] + Bc*par[1] + Bd*par[2]*par[1] + Be*par[1]*par[1] + Bf*par[2]*par[2]*par[1]*par[1];

    Double_t q0 = 0.07;
    Double_t n = 20.0;

    Double_t E = 1. + A2 * exp( - B2 * x[0] ); // exponential smoothing
    Double_t F = 1. / (1. + pow(x[0]/q0,n)); // cutoff factor

    Double_t kgamow = KGamow(x, par);
    Double_t kmod = KMod(x, par, par2);
    return 1./(F*(1./kgamow)*(1./kmod) + (1.-F)*E); 
}

// Modified Coulomb Correction
Double_t KMod(Double_t *x, Double_t *par, Double_t *par2){  // Direct implementation from paper
    return 1. + (par2[0]*((PREFACTOR*par[1])/(par[2])))/(1. + par2[1]*x[0]*par[1]/(par[2]*HBARC) + par2[2]*pow(x[0]*par[1]/(par[2]*HBARC),2) + par2[3]*pow(x[0]*par[1]/(par[2]*HBARC),4));
}

// Function to calculate qinv from two 4-momenta
Double_t GetQ(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2) {
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Sum4V = p1+p2;
    Double_t q = Sum4V.M2() - 4*PI_MASS*PI_MASS;
    return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

#endif