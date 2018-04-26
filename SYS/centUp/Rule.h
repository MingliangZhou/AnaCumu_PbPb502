#ifndef Rule_H_
#define Rule_H_

#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TComplex.h"
#include "TMath.h"
#include "TRandom3.h"

char name[300];

const double cutCent[80] = { // centrality cut from 5.02 TeV Pb+Pb
0.069136, 0.07531,  0.081906, 0.088954, 0.09648,  0.104536, 0.113101, 0.122216, 0.131908, 0.142242,
0.15318,  0.164806, 0.177083, 0.190096, 0.203824, 0.218281, 0.233553, 0.249597, 0.266458, 0.284192,
0.302826, 0.322382, 0.342787, 0.364095, 0.386397, 0.409774, 0.434162, 0.459497, 0.485844, 0.513454,
0.542107, 0.571956, 0.602869, 0.634998, 0.668303, 0.702811, 0.738626, 0.775749, 0.814217, 0.85383,
0.894902, 0.937428, 0.981313, 1.02683,  1.07376,  1.12226,  1.17244,  1.22398,  1.27719,  1.33218,
1.38892,  1.44736,  1.50764,  1.57,     1.63416,  1.70029,  1.76876,  1.83945,  1.91222,  1.98743,
2.06495,  2.14525,  2.22772,  2.31305,  2.40148,  2.49292,  2.58742,  2.68513,  2.78688,  2.89266,
3.00247,  3.11632,  3.23537,  3.35982,  3.49022,  3.62674,  3.77047,  3.92323,  4.08541,  4.2649
};

const unsigned int maxNtrk = 10000; // maximum number of tracks in input
const double cutEta = 2.5; // track eta cut
const double cutEtaGap = 0.; // eta gap cut for subevents
const double cutZvtx = 100; // zVtx position cut
const double cutQual = 1; // track quality cut, even tight

const unsigned int nBin = 80; // number of bins for cumulant calculation
const double minBin = 0; // minimum centrality
const double maxBin = 80; // maximum centrality
const unsigned int nRebin[2] = {8, 16}; // number of bins for final results

const unsigned int NV = 5; // harmonics from v0 to v4
const unsigned int NVm = 13; // higher harmonics for Q-cumulant with weights, from v0 to v12
const unsigned int NP = 2; // number of pT cuts for differential particles
const double minPt[NP] = {0.5, 1.0}; // minimum pT for reference particles
const double maxPt[NP] = {5.0, 5.0}; // maximum pT for reference particles
const unsigned int NW = 7; // power of weights for Q-cumulant calcualtion: 0 to 6
const unsigned int NA = 3; // number of permutation for 3-subevent

const unsigned int nMixZvtx = 5; // number of zvtx bins for mixed events
const unsigned int nMixCent = 10; // number of centralities for mixed events
const unsigned int nDepth = 4; // depth of mix event pool

const unsigned int nSample = 20; // number of samples for statistical errors
const unsigned int nFile = 33; // number of files per bundle to calculate statistical errors

#endif


