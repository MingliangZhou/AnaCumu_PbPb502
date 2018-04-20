#ifndef Plot_H_
#define Plot_H_

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
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"

char name[200];

const int mS[8] = {20,21,33,34,29,22,23,27};
const int mC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};
const int lS[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
const int lC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};

const unsigned int nCent = 7;
const double binCent[nCent+1] = {0, 5, 10, 20, 30, 40, 60, 80};
const unsigned int nPt = 6;
const double binPt[nPt+1] = {0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0};

class Plot
{
	private:
		TH1D* hEtaEff[nCent][nPt];
		TH1D* hEtaFak[nCent][nPt];
		void readHist_CP(TFile*, TH1D*&, const char*, int iC, int iP);
		void draw_graph(vector<TH1D*>, int iP=-1, int iOpt=-1);
		void styleGraph(TGraph*, int);
		void styleGraph(TH1D*, int);
		void getYrange(TGraph*, double&, double&);

	public:
		Plot();
		~Plot();
		void initialize();
		void execute();
		void finalize();
};

#endif
