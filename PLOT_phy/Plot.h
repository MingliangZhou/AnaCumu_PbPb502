#ifndef Plot_H_
#define Plot_H_

#include "../MAIN_binCent/Rule.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"

const unsigned int NF = 3; // 3 binning cases

const int mS[8] = {20,21,33,34,29,22,23,27};
const int mC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};
const int lS[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
const int lC[8] = {kBlue+2, kRed+2, kGreen+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2, 1};

class Plot
{
	private:
		TGraphErrors* sts_c2_1sub[NF][NV][NP];
		TGraphErrors* sts_c4_1sub[NF][NV][NP];
		TGraphErrors* sts_c6_1sub[NF][NV][NP];
		TGraphErrors* sts_nc4_1sub[NF][NV][NP];
		TGraphErrors* sts_nc6_1sub[NF][NV][NP];
		TGraphErrors* sts_sc_1sub[NF][NV][NP];
		TGraphErrors* sts_nsc_1sub[NF][NV][NP];
		TGraphErrors* sts_ac_1sub[NF][NV][NP];
		TGraphErrors* sts_nac_1sub[NF][NV][NP];
		TGraphErrors* sts_isGauss_1sub[NF][NV][NP];
		TGraphErrors* sts_isPower_1sub[NF][NV][NP];
		TGraphErrors* sts_c2_3sub[NF][NV][NP];
		TGraphErrors* sts_c4_3sub[NF][NV][NP];
		TGraphErrors* sts_nc4_3sub[NF][NV][NP];
		TGraphErrors* sts_sc_3sub[NF][NV][NP];
		TGraphErrors* sts_nsc_3sub[NF][NV][NP];
		TGraphErrors* sts_ac_3sub[NF][NV][NP];
		TGraphErrors* sts_nac_3sub[NF][NV][NP];

		TGraphAsymmErrors* sys_c2_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_c4_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_c6_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nc4_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nc6_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_sc_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nsc_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_ac_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nac_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_isGauss_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_isPower_1sub[NF][NV][NP];
		TGraphAsymmErrors* sys_c2_3sub[NF][NV][NP];
		TGraphAsymmErrors* sys_c4_3sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nc4_3sub[NF][NV][NP];
		TGraphAsymmErrors* sys_sc_3sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nsc_3sub[NF][NV][NP];
		TGraphAsymmErrors* sys_ac_3sub[NF][NV][NP];
		TGraphAsymmErrors* sys_nac_3sub[NF][NV][NP];

		void readHist_VP(TFile*, TGraphErrors*&, const char*, int iV, int iP);
		void readHist_VP(TFile*, TGraphAsymmErrors*&, const char*, int iV, int iP);
		void draw_graph_4by2(vector<TGraphErrors*>, vector<TGraphAsymmErrors*>, int iOpt=-1);
		void draw_graph_4by3(vector<TGraphErrors*>, vector<TGraphAsymmErrors*>, int iV=-1, int iOpt=-1);
		void styleGraph(TGraph*, int);
		void styleGraph(TH1D*, int);
		void getYrange(TGraph*, double&, double&, int);
		void reverseXaxis(TH1D*, double, double, double, bool);

	public:
		Plot();
		~Plot();
		void initialize();
		void execute();
		void finalize();
};

#endif
