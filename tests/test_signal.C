#include <iostream>
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"



void test_signal() {
    Float_t pionMass;
    TLorentzVector TrackFourVectorA0;
    TLorentzVector TrackFourVectorA1;
    TLorentzVector TrackFourVectorB0;
    TLorentzVector TrackFourVectorB1;
    TLorentzVector TrackFourVectorB2;
    TLorentzVector TrackFourVectorB3;
    TLorentzVector TrackFourVectorC0;
    TLorentzVector TrackFourVectorC1;
    TLorentzVector TrackFourVectorC2;
    TLorentzVector TrackFourVectorD0;
    TLorentzVector TrackFourVectorD1;
    TLorentzVector TrackFourVectorD2;
    TLorentzVector TrackFourVectorD3;

    std::vector<TLorentzVector> TrackFourVectorA;
    std::vector<TLorentzVector> TrackFourVectorB;
    std::vector<TLorentzVector> TrackFourVectorC;
    std::vector<TLorentzVector> TrackFourVectorD;
    
    std::vector<std::vector<TLorentzVector>> AllTrackFourVector;
    std::vector<std::vector<TLorentzVector>> testAllTrackFourVector;

    std::vector<double> valuesA = {0, 1};
    std::vector<double> valuesB = {2, 3, 4, 5};
    std::vector<double> valuesC = {6, 7, 8};
    std::vector<double> valuesD = {9, 10, 11, 12};
    std::vector<std::vector<double>> values = {valuesA, valuesB, valuesC, valuesD};


    //double values[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    TrackFourVectorA0.SetPtEtaPhiM(values[0][0], values[0][0], values[0][0], values[0][0]);
    TrackFourVectorA1.SetPtEtaPhiM(values[0][1], values[0][1], values[0][1], values[0][1]);

    TrackFourVectorB0.SetPtEtaPhiM(values[1][0], values[1][0], values[1][0], values[1][0]);
    TrackFourVectorB1.SetPtEtaPhiM(values[1][1], values[1][1], values[1][1], values[1][1]);
    TrackFourVectorB2.SetPtEtaPhiM(values[1][2], values[1][2], values[1][2], values[1][2]);
    TrackFourVectorB3.SetPtEtaPhiM(values[1][3], values[1][3], values[1][3], values[1][3]);

    TrackFourVectorC0.SetPtEtaPhiM(values[2][0], values[2][0], values[2][0], values[2][0]);
    TrackFourVectorC1.SetPtEtaPhiM(values[2][1], values[2][1], values[2][1], values[2][1]);
    TrackFourVectorC2.SetPtEtaPhiM(values[2][2], values[2][2], values[2][2], values[2][2]);

    TrackFourVectorD0.SetPtEtaPhiM(values[3][0], values[3][0], values[3][0], values[3][0]);
    TrackFourVectorD1.SetPtEtaPhiM(values[3][1], values[3][1], values[3][1], values[3][1]);
    TrackFourVectorD2.SetPtEtaPhiM(values[3][2], values[3][2], values[3][2], values[3][2]);
    TrackFourVectorD3.SetPtEtaPhiM(values[3][3], values[3][3], values[3][3], values[3][3]);

    TrackFourVectorA.push_back(TrackFourVectorA0);
    TrackFourVectorA.push_back(TrackFourVectorA1);
    TrackFourVectorB.push_back(TrackFourVectorB0);
    TrackFourVectorB.push_back(TrackFourVectorB1);
    TrackFourVectorB.push_back(TrackFourVectorB2);
    TrackFourVectorB.push_back(TrackFourVectorB3);
    TrackFourVectorC.push_back(TrackFourVectorC0);
    TrackFourVectorC.push_back(TrackFourVectorC1);
    TrackFourVectorC.push_back(TrackFourVectorC2);
    TrackFourVectorD.push_back(TrackFourVectorD0);
    TrackFourVectorD.push_back(TrackFourVectorD1);
    TrackFourVectorD.push_back(TrackFourVectorD2);
    TrackFourVectorD.push_back(TrackFourVectorD3);

    testAllTrackFourVector.push_back(TrackFourVectorA);
    testAllTrackFourVector.push_back(TrackFourVectorB);
    testAllTrackFourVector.push_back(TrackFourVectorC);
    testAllTrackFourVector.push_back(TrackFourVectorD);
    
    TStopwatch stopwatch, stopwatch1, stopwatch2;
    
    std::cout << "antes: " << AllTrackFourVector.size() << std::endl;

    int totalElements = 0;

    std::cout << "0" << std::endl;
    // --- POOL SETTINGS ---
    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> poolEvents;
    std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>>> phiEvents;
    int poolSize = 3;
    int poolMaxTracks = 0;
    int phiMaxTracks = 0;
    int phiSize = poolSize - (testAllTrackFourVector.size() % poolSize);
    bool phiSwitch = true;
    if (phiSize < 1) phiSwitch = false;
    
    int phiSteps;

    if (phiSwitch) {
        phiSteps = (testAllTrackFourVector.size() - (poolSize - phiSize)) / phiSize; 
    } else {
        phiSteps = 1;
    }

    int x = 0;
    // ------

    TBenchmark benchmark;

    for (Long64_t i = 0; i < testAllTrackFourVector.size(); i++) {

        std::vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>> TrackFourVector;
        for (int j = 0; j < testAllTrackFourVector[i].size(); j++) {
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> vec(values[i][j], values[i][j], values[i][j], pionMass);
            TrackFourVector.push_back(vec);
            // --- MIXING ---
            // Double-loop method
            // começar se o tamanho do vetor de eventos é igual ao poolSize -1
            if (poolEvents.size() == poolSize - 1) {
                // loop entre os eventos dentro do vetor de eventos
                for (int k = 0; k < poolEvents.size(); k++) {
                    // loop entre as tracks do evento
                    for (int z = 0; z < poolEvents[k].size(); z++) {
                    // mixar com a track que acamos de montar (vec)
                        std::cout << "mix2> i: " << i << " j: " << j << " || k: " << k << " z: " << z << std::endl;
                        std::cout << "Track A: " << vec << std::endl;
                        std::cout << "Track B: " << poolEvents[k][z] << std::endl; 
                    }
                }

                // Single-loop method
                for (int w = 0; w < poolMaxTracks*poolEvents.size(); w++) {
                    int k = w / poolMaxTracks;
                    int z = w % poolMaxTracks;

                    if (z > poolEvents[k].size() - 1) continue;
                    std::cout << "mix1> i: " << i << " j: " << j << " || k: " << k << " z: " << z << std::endl; 
                    std::cout << "Track A: " << vec << std::endl;
                    std::cout << "Track B: " << poolEvents[k][z] << std::endl;
                }
            }
        }
        // --- Phi Handler ---
        if (x == phiSteps && phiSwitch && phiEvents.size() < phiSize) {
            phiEvents.push_back(TrackFourVector);
            if (TrackFourVector.size() > phiMaxTracks) phiMaxTracks = TrackFourVector.size();
            x = 0;
        }

        // --- Pool Handler
        poolEvents.push_back(TrackFourVector);
        if (TrackFourVector.size() > poolMaxTracks) poolMaxTracks = TrackFourVector.size();

        // depois de mixar resetamos o vetor de eventos para proxima poool
        if (poolEvents.size() == poolSize) {
            poolEvents.clear();
            poolMaxTracks = 0;
        }

        // --- SIGNAL ---
        if (TrackFourVector.size() > 1) {
            // Double-loop method
            for (size_t p1 = 0; p1 < TrackFourVector.size(); p1++) {
                for (size_t p2 = p1 + 1; p2 < TrackFourVector.size(); p2++) {
                    std::cout << "sig2> i: " << i << " p1: " << p1 << " p2: " << p2 << std::endl;
                    std::cout << "Track A " << TrackFourVector[p1] << std::endl;
                    std::cout << "Track B " << TrackFourVector[p2] << std::endl;
                }
            }

            // Single-loop method
            stopwatch2.Start(kFALSE);
            for (int k = 0; k < TrackFourVector.size()*(TrackFourVector.size() - 1); k++) {
                int p1 = k / (TrackFourVector.size() - 1);
                int p2 = (k % (TrackFourVector.size() - 1))+1;

                if (p1 >= p2) continue;
                std::cout << "sig1> i: " << i << " p1: " << p1 << " p2: " << p2 << std::endl;
                std::cout << "Track A " << TrackFourVector[p1] << std::endl;
                std::cout << "Track B " << TrackFourVector[p2] << std::endl;
            }
        }
        x++;
    }

    // --- PHI MIXING ---
    if (phiSwitch) {
        // add pool left overs to phi
        for (int event = 0; event < poolEvents.size(); event++) {
            phiEvents.push_back(poolEvents[event]);
        }

        if (poolMaxTracks > phiMaxTracks) phiMaxTracks = poolMaxTracks;

        int lastIndex = phiEvents.size() - 1;
        int lastEventSize = phiEvents[lastIndex].size();

        // loop over last event tracks to mix
        for (int lastEventTrackID = 0; lastEventTrackID < lastEventSize; lastEventTrackID++) {
            // -- MIXING ---
            for (int k = 0; k < phiEvents.size() - 1; k++) {
                for (int z = 0; z < phiEvents[k].size(); z++) {
                    std::cout << "phimix2> TrackID: " << lastEventTrackID << " || k: " << k << " z: " << z << std::endl; 
                    std::cout << "track A: " << phiEvents[lastIndex][lastEventTrackID] << std::endl;
                    std::cout << "track B: " << phiEvents[k][z] << std::endl;
                }
            }

            // Single-loop method
            for (int w = 0; w < phiMaxTracks*(phiEvents.size()-1); w++) {
                int k = w / phiMaxTracks;
                int z = w % phiMaxTracks;

                if (z > phiEvents[k].size() - 1) continue;
                std::cout << "phimix1> TrackID: " << lastEventTrackID << " || k: " << k << " z: " << z << std::endl; 
                std::cout << "track A: " << phiEvents[lastIndex][lastEventTrackID] << std::endl;
                std::cout << "track B: " << phiEvents[k][z] << std::endl;
            }
        }
    }
    
    /*benchmark.Stop("4-Vector");

    std::cout << "depois: " << AllTrackFourVector.size() << std::endl;

    std::cout << std::endl;
    std::cout << "1" << std::endl;
    benchmark.Start("TripleLoop");
    for (int i = 0; i < AllTrackFourVector.size(); i++) {
        for (size_t p1 = 0; p1 < AllTrackFourVector[i].size(); p1++) {
            for (size_t p2 = p1 + 1; p2 < AllTrackFourVector[i].size(); p2++) {
                //AllTrackFourVector[i][p1].Print();
                //AllTrackFourVector[i][p2].Print();
            }
        }
    }
    benchmark.Stop("TripleLoop");

    std::cout << std::endl;
    std::cout << "2" << std::endl;
    benchmark.Start("DoubleLoop");
    
    for (int i = 0; i < AllTrackFourVector.size(); i++) {
        
        for (int k = 0; k < AllTrackFourVector[i].size()*(AllTrackFourVector[i].size() - 1); k++) {
            int p1 = k / (AllTrackFourVector[i].size() - 1);
            int p2 = (k % (AllTrackFourVector[i].size() - 1)) + 1;

            if (p2 > p1) {
                std::cout << "i: " << i << " p1: " << p1 << " p2: " << p2 << std::endl;
                //AllTrackFourVector[i][p1].Print();
                //AllTrackFourVector[i][p2].Print();
            }
        }
    }
    benchmark.Stop("DoubleLoop");

    std::cout << std::endl;
    std::cout << "3" << std::endl;
    std::cout << "total elements: " << totalElements << std::endl;
    benchmark.Start("SingleLoop");

    for (int j = 0, a = 0, b = 0, i = 0; j < totalElements; j++) {
        int trigger = AllTrackFourVector[i].size() * (AllTrackFourVector[i].size() - 1);

        int k = j - a;
        
        int p1 = k / (AllTrackFourVector[i].size() - 1);
        int p2 = (k % (AllTrackFourVector[i].size() - 1)) + 1;
        
        if (p2 > p1) {
            std::cout << "i: " << i << " p1: " << p1 << " p2: " << p2 << std::endl;
            //AllTrackFourVector[i][p1].Print();
            //AllTrackFourVector[i][p2].Print();
        }

        b++;
        if (b == trigger) {
            a = j + 1;
            b = 0;
            i++;
        }
    }
    benchmark.Stop("SingleLoop");

     std::cout << "Benchmark Results:" << std::endl;
    benchmark.Show("4-Vector");
    benchmark.Show("TripleLoop");
    benchmark.Show("DoubleLoop");
    benchmark.Show("SingleLoop"); */
}


/*
if (eventTracks.size() >= poolSize) { 
            for (int i = 0; i < eventTracks.size() - 1; i++) {
                for (int j = 0; j < eventTracks[i].size(); j++) {
                    for (int k = i+1; k < eventTracks.size(); k++) {
                        for (int z = 0; z < eventTracks[k].size(); z++) {
                            std::cout << "i: " << i << " j: " << j << " || k: " << k << " z: " << z << std::endl; 
                        }
                    }
                }
            }
            eventTracks.clear();
        }
*/
