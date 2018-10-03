#ifndef Plot_H_
#define Plot_H_

#include "Rule.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"

const int mS[8] = {20,21,33,34,29,22,23,27};
const int mC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};
const int lS[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
const int lC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};

class Plot
{
	private:
		TGraphErrors* c2_1sub[NV][NP];
		TGraphErrors* c4_1sub[NV][NP];
		TGraphErrors* c6_1sub[NV][NP];
		TGraphErrors* nc4_1sub[NV][NP];
		TGraphErrors* nc6_1sub[NV][NP];
		TGraphErrors* sc_1sub[NV][NP];
		TGraphErrors* nsc_1sub[NV][NP];
		TGraphErrors* ac_1sub[NV][NP];
		TGraphErrors* nac_1sub[NV][NP];
		TGraphErrors* isGauss_1sub[NV][NP];
		TGraphErrors* isPower_1sub[NV][NP];
		TGraphErrors* c2_3sub[NV][NP];
		TGraphErrors* c4_3sub[NV][NP];
		TGraphErrors* nc4_3sub[NV][NP];
		TGraphErrors* sc_3sub[NV][NP];
		TGraphErrors* nsc_3sub[NV][NP];
		TGraphErrors* ac_3sub[NV][NP];
		TGraphErrors* nac_3sub[NV][NP];
		void readHist_VP(TFile*, TGraphErrors*&, const char*, int iV, int iP);
		void draw_graph(vector<TGraphErrors*>, int iV=-1, int iP=-1, int iOpt=-1);
		void styleGraph(TGraph*, int);
		void styleGraph(TH1D*, int);
		void getYrange(TGraph*, double&, double&, int);

	public:
		Plot(unsigned int iBin);
		~Plot();
		void initialize(unsigned int iBin);
		void execute(unsigned int iBin);
		void finalize();
};

#endif
