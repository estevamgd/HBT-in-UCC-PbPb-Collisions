#ifndef DATA_FUNC_H
#define DATA_FUNC_H

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "TMath.h"
#include <TF1.h>

// Pion mass
Float_t pionMass = 0.13957039;

// to reject a range in the fit -- in principle did not reject any range
Double_t reject_range_min = 0.01;
Double_t reject_range_max = 0.00001;

// Exponential function + long range
Double_t func1_exp(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-par[2]*x[0]/0.1973))*(1+par[3]*x[0]);}
    return v;
}

// Gaussian function + long range
Double_t func2_gauss(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,2.)))*(1+par[3]*x[0]);}
    return v;
} 

// Levy function 
Double_t func3_levy(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,par[4])))*(1+par[3]*x[0]);}
    return v;
}

// Function to calculate qinv from two 4-momenta
Double_t GetQ(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2) {
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Sum4V = p1+p2;
    Double_t q = Sum4V.M2() - 4*pionMass*pionMass;
    return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

#endif