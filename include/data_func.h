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
Double_t FitExp(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-par[2]*x[0]/0.1973))*(1+par[3]*x[0]);}
    return v;
}

// Gaussian function + long range
Double_t FitGauss(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,2.)))*(1+par[3]*x[0]);}
    return v;
} 

// Levy function 
Double_t FitLevy(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,par[4])))*(1+par[3]*x[0]);}
    return v;
}

// Levy function (direct implementation)
Double_t FitLevy2(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= 1 - par[0] + par[0]*(1 + exp(-pow(par[1]*x[0],par[2])));} // We still need to add the coulomb correction outside
    return v;
}

// Background function 
Double_t FitBG(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0],2)))*(1 - par[3]*exp(-pow(par[4]*x[0],2)));}
    return v;
}

// Double Ratio function 
Double_t FitLevyDR(Double_t* x, Double_t* par){
    Double_t v = 0;
    if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
    else{v= par[0]*(1 - par[1] + par[1]*(1 + exp(-pow(par[2]*x[0],par[3]))))*(1+par[4]*x[0]);} // We still need to add the coulomb correction outside
    return v;
}

// Function to calculate qinv from two 4-momenta
Double_t GetQ(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p1, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> &p2) {
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Sum4V = p1+p2;
    Double_t q = Sum4V.M2() - 4*pionMass*pionMass;
    return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

#endif