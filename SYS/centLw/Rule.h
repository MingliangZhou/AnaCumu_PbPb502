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
0.058527, 0.064057, 0.069999, 0.076379, 0.083214, 0.09052,  0.098342, 0.106732, 0.115651, 0.12515,
0.135286, 0.146078, 0.157514, 0.169673, 0.182579, 0.196214, 0.210613, 0.225837, 0.241851, 0.258713,
0.276461, 0.295168, 0.314811, 0.335361, 0.356885, 0.379373, 0.402984, 0.427675, 0.453331, 0.480118,
0.508105, 0.537254, 0.567621, 0.599129, 0.631882, 0.665879, 0.701144, 0.737746, 0.775749, 0.815154,
0.855791, 0.897878, 0.941579, 0.986655, 1.0334,   1.08171,  1.13168,  1.18333,  1.23644,  1.29143,
1.34812,  1.40677,  1.46718,  1.52974,  1.59422,  1.66063,  1.7295,   1.80042,  1.87379,  1.94949,
2.0277,   2.10858,  2.19205,  2.27819,  2.36745,  2.45993,  2.55558,  2.65451,  2.75761,  2.86439,
2.9759,   3.09149,  3.21218,  3.33876,  3.47119,  3.60999,  3.75643,  3.91194,  4.07724,  4.26017
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


