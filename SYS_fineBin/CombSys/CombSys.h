#ifndef CombSys_H_
#define CombSys_H_

#include "../default/Rule.h"

const unsigned int NS = 9; // number of systematic sources
// 0: default for statistical errors
// 1: low tracking efficiency
// 2: high tracking efficiency
// 3: tight track selection
// 4: loose pileup rejection
// 5: MC closure
// 6: no flattening
// 7: 84% centrality
// 8: 86% centrality

class CombSys
{
	private:
		TGraphErrors* c2_1sub[2][NS][NV][NP];
		TGraphErrors* c4_1sub[2][NS][NV][NP];
		TGraphErrors* c6_1sub[2][NS][NV][NP];
		TGraphErrors* nc4_1sub[2][NS][NV][NP];
		TGraphErrors* nc6_1sub[2][NS][NV][NP];
		TGraphErrors* sc_1sub[2][NS][NV][NP];
		TGraphErrors* nsc_1sub[2][NS][NV][NP];
		TGraphErrors* ac_1sub[2][NS][NV][NP];
		TGraphErrors* nac_1sub[2][NS][NV][NP];
		TGraphErrors* isGauss_1sub[2][NS][NV][NP];
		TGraphErrors* isPower_1sub[2][NS][NV][NP];
		TGraphErrors* c2_3sub[2][NS][NV][NP];
		TGraphErrors* c4_3sub[2][NS][NV][NP];
		TGraphErrors* nc4_3sub[2][NS][NV][NP];
		TGraphErrors* sc_3sub[2][NS][NV][NP];
		TGraphErrors* nsc_3sub[2][NS][NV][NP];
		TGraphErrors* ac_3sub[2][NS][NV][NP];
		TGraphErrors* nac_3sub[2][NS][NV][NP];

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

		TGraphErrors* ratioSm_c2_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_c4_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_c6_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_nc4_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_nc6_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_sc_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_nsc_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_ac_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_nac_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_isGauss_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_isPower_1sub[NS][NV][NP];
		TGraphErrors* ratioSm_c2_3sub[NS][NV][NP];
		TGraphErrors* ratioSm_c4_3sub[NS][NV][NP];
		TGraphErrors* ratioSm_nc4_3sub[NS][NV][NP];
		TGraphErrors* ratioSm_sc_3sub[NS][NV][NP];
		TGraphErrors* ratioSm_nsc_3sub[NS][NV][NP];
		TGraphErrors* ratioSm_ac_3sub[NS][NV][NP];
		TGraphErrors* ratioSm_nac_3sub[NS][NV][NP];

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

		void readHist_VP(TFile*, TGraphErrors*&, const char*, const char*, int iF, int iS, int iV, int iP);
		void cal_ratio(TGraphErrors* gChk, TGraphErrors* gDef, TGraphErrors*& gRatio, const char*, int iS, int iV, int iP);
		void smooth_ratio(TGraphErrors* gRatio, int iB, double percent);
		void reduce_ratio(TGraphErrors* gRatio, int iB, double percent);
		void cal_comb(vector<TGraphErrors*> gVec, TGraphErrors*& gSysUp, TGraphErrors*& gSysLw, const char*, int iV, int iP);
		void smooth();

	public:
		CombSys(unsigned int iBin);
		~CombSys();
		void initialize(unsigned int iBin);
		void execute();
		void finalize(unsigned int iBin);
};


#endif
