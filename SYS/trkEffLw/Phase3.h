#ifndef Phase3_H_
#define Phase3_H_

#include "Rule.h"

class Phase3
{
	private:
		// input
		TH1D* cnt_1sub[NP];
		TH1D* rbn_c2_1sub[nSample+1][NV][NP];
		TH1D* rbn_c4_1sub[nSample+1][NV][NP];
		TH1D* rbn_c6_1sub[nSample+1][NV][NP];
		TH1D* rbn_nc4_1sub[nSample+1][NV][NP];
		TH1D* rbn_nc6_1sub[nSample+1][NV][NP];
		TH1D* rbn_sc_1sub[nSample+1][NV][NP];
		TH1D* rbn_nsc_1sub[nSample+1][NV][NP];
		TH1D* rbn_ac_1sub[nSample+1][NV][NP];
		TH1D* rbn_nac_1sub[nSample+1][NV][NP];
		TH1D* rbn_isGauss_1sub[nSample+1][NV][NP];
		TH1D* rbn_isPower_1sub[nSample+1][NV][NP];
		TH1D* cnt_3sub[NP];
		TH1D* rbn_cnt_3sub[nSample+1][NP];
		TH1D* rbn_c2_3sub[nSample+1][NV][NP];
		TH1D* rbn_c4_3sub[nSample+1][NV][NP];
		TH1D* rbn_nc4_3sub[nSample+1][NV][NP];
		TH1D* rbn_sc_3sub[nSample+1][NV][NP];
		TH1D* rbn_nsc_3sub[nSample+1][NV][NP];
		TH1D* rbn_ac_3sub[nSample+1][NV][NP];
		TH1D* rbn_nac_3sub[nSample+1][NV][NP];
		
		// output
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

		void cal_sts(vector<TH1D*>, TGraphErrors*&, const char*, int iV, int iP);
		void cleanPoint(TGraphErrors*);
		void readHist_VP(TFile*, TH1D*&, const char*, int iV, int iP);

	public:
		Phase3(unsigned int iBin);
		~Phase3();
		void initialize(unsigned int iBin);
		void execute();
		void finalize(unsigned int iBin);
};

#endif
