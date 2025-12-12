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

// Levy function (direct implementation)
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

// Coulomb Correction from: https://github.com/csanadm/coulcorrlevyparam/blob/0c3d63e765762ec06d0eb66e42b95b8561e6a669/coulcorr_param.cc#L16-L81
Double_t KCoulomb(Double_t *x, Double_t *par){
    double q = x[0];
    double R = par[1];
    double alpha = par[2];
   
    double aA =  0.26984;
    double aB = -0.491231;
    double aC =  0.0352313;
    double aD = -1.31628;
    double aE =  0.00359148;

    double bA =  2.37267;
    double bB =  0.586309;
    double bC =  2.24867;
    double bD = -1.43278;
    double bE = -0.0521642;
    double bF =  0.729434;

    double cA = -4.30347 ;
    double cB =  1.17849e-05;
    double cC =  3.30346;
    double cD =  1.27273e-06 ;
    double cE =  3.03399e-06;
    double cF =  1.68883;

    double dA =  0.000568486;
    double dB = -0.805271;
    double dC = -0.192606 ;
    double dD =  2.77504;
    double dE =  2.02951;
    double dF =  1.07906; 
   
    double a_AR = pow((aA * alpha + aB),2) + pow((aC * R + aD),2) + aE * pow((alpha * R + 1),2);
    double b_AR = (1. + bA*pow(R,bB) -pow(alpha,bC)) / (alpha*alpha*R*(pow(alpha,bD)+bE*pow(R,bF)));
    double c_AR = (cA  + pow(alpha,cB) + cC*pow(R,cD))/cE*pow(alpha/R,cF);
    double d_AR = (dA + (pow(R,dB) + dC*pow(alpha,dF))/(pow(R,dD)*pow(alpha,dE)));
    double t = R*q/HBARC/alpha;
    double parametrization = (1. + PREFACTOR*a_AR*R/alpha/(1.+b_AR*t+c_AR*t*t+d_AR*t*t*t*t));
       
    double Aa =  0.126253;
    double Ab =  0.053848;
    double Ac = -0.00912627;
    double Ad = -0.018459;
    double Ae =  0.000851755;
    double Af =  0.000417179;

    double Ba = 19.3162;
    double Bb =  5.58961;
    double Bc =  2.26264;
    double Bd = -1.28486;
    double Be = -0.0821616;
    double Bf =  0.0238446; 
   
    double A = Aa + Ab*alpha + Ac*R + Ad*alpha*R + Ae*R*R + Af*alpha*alpha*R*R;
    double B = Ba + Bb*alpha + Bc*R + Bd*alpha*R + Be*R*R + Bf*alpha*alpha*R*R;
    double exponential_smoothing = 1. + A * exp( - B * q );
   
    double q0 = 0.07;
    double n = 20.0;
    double cutoff = 1. / ( 1.+pow(q/q0,n) );
   
    return 1./(  (1. - cutoff) * exponential_smoothing + cutoff * (1 / KGamow(x,par)) * (1 / parametrization));
}

// Function to calculate qinv from two 4-momenta
Double_t GetQ(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2) {
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Sum4V = p1+p2;
    Double_t q = Sum4V.M2() - 4*PI_MASS*PI_MASS;
    return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

#endif