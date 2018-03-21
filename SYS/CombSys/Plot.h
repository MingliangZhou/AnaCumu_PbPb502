#ifndef Plot_H_
#define Plot_H_

#include "../default/Rule.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"

const unsigned int NS = 7; // number of systematic sources

const int mS[8] = {20,21,33,34,29,22,23,27};
const int mC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};
const int lS[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
const int lC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};

class Plot
{
	private:
		TGraphErrors* ratio_c2_1sub[NS][NV][NP];
		TGraphErrors* ratio_c4_1sub[NS][NV][NP];
		TGraphErrors* ratio_c6_1sub[NS][NV][NP];
		TGraphErrors* ratio_nc4_1sub[NS][NV][NP];
		TGraphErrors* ratio_nc6_1sub[NS][NV][NP];
		TGraphErrors* ratio_sc_1sub[NS][NV][NP];
		TGraphErrors* ratio_nsc_1sub[NS][NV][NP];
		TGraphErrors* ratio_ac_1sub[NS][NV][NP];
		TGraphErrors* ratio_nac_1sub[NS][NV][NP];
		TGraphErrors* ratio_isGauss_1sub[NS][NV][NP];
		TGraphErrors* ratio_isPower_1sub[NS][NV][NP];
		TGraphErrors* ratio_c2_3sub[NS][NV][NP];
		TGraphErrors* ratio_c4_3sub[NS][NV][NP];
		TGraphErrors* ratio_nc4_3sub[NS][NV][NP];
		TGraphErrors* ratio_sc_3sub[NS][NV][NP];
		TGraphErrors* ratio_nsc_3sub[NS][NV][NP];
		TGraphErrors* ratio_ac_3sub[NS][NV][NP];
		TGraphErrors* ratio_nac_3sub[NS][NV][NP];
		
		TGraphErrors* sysUp_c2_1sub[NV][NP];
		TGraphErrors* sysUp_c4_1sub[NV][NP];
		TGraphErrors* sysUp_c6_1sub[NV][NP];
		TGraphErrors* sysUp_nc4_1sub[NV][NP];
		TGraphErrors* sysUp_nc6_1sub[NV][NP];
		TGraphErrors* sysUp_sc_1sub[NV][NP];
		TGraphErrors* sysUp_nsc_1sub[NV][NP];
		TGraphErrors* sysUp_ac_1sub[NV][NP];
		TGraphErrors* sysUp_nac_1sub[NV][NP];
		TGraphErrors* sysUp_isGauss_1sub[NV][NP];
		TGraphErrors* sysUp_isPower_1sub[NV][NP];
		TGraphErrors* sysUp_c2_3sub[NV][NP];
		TGraphErrors* sysUp_c4_3sub[NV][NP];
		TGraphErrors* sysUp_nc4_3sub[NV][NP];
		TGraphErrors* sysUp_sc_3sub[NV][NP];
		TGraphErrors* sysUp_nsc_3sub[NV][NP];
		TGraphErrors* sysUp_ac_3sub[NV][NP];
		TGraphErrors* sysUp_nac_3sub[NV][NP];

		TGraphErrors* sysLw_c2_1sub[NV][NP];
		TGraphErrors* sysLw_c4_1sub[NV][NP];
		TGraphErrors* sysLw_c6_1sub[NV][NP];
		TGraphErrors* sysLw_nc4_1sub[NV][NP];
		TGraphErrors* sysLw_nc6_1sub[NV][NP];
		TGraphErrors* sysLw_sc_1sub[NV][NP];
		TGraphErrors* sysLw_nsc_1sub[NV][NP];
		TGraphErrors* sysLw_ac_1sub[NV][NP];
		TGraphErrors* sysLw_nac_1sub[NV][NP];
		TGraphErrors* sysLw_isGauss_1sub[NV][NP];
		TGraphErrors* sysLw_isPower_1sub[NV][NP];
		TGraphErrors* sysLw_c2_3sub[NV][NP];
		TGraphErrors* sysLw_c4_3sub[NV][NP];
		TGraphErrors* sysLw_nc4_3sub[NV][NP];
		TGraphErrors* sysLw_sc_3sub[NV][NP];
		TGraphErrors* sysLw_nsc_3sub[NV][NP];
		TGraphErrors* sysLw_ac_3sub[NV][NP];
		TGraphErrors* sysLw_nac_3sub[NV][NP];

		void readHist_SVP(TFile*, TGraphErrors*&, const char*, int iS, int iV, int iP);
		void readHist_VP(TFile*, TGraphErrors*&, const char*, int iV, int iP);
		void draw_graph(vector<TGraphErrors*>, int iV=-1, int iP=-1, int iOpt=-1, unsigned int iBin=2);
		void styleGraph(TGraph*, int);
		void styleGraph(TH1D*, int);
		void getYrange(TGraph*, double&, double&);

	public:
		Plot(unsigned int iBin);
		~Plot();
		void initialize(unsigned int iBin);
		void execute(unsigned int iBin);
		void finalize();
};

#endif
