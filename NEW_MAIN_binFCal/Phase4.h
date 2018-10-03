#ifndef Phase4_H_
#define Phase4_H_

#include "Rule.h"

class Phase4
{
	private:
		TGraphErrors* gCvt;

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

		TGraphErrors* sts_c2_1sub[NV][NP];
		TGraphErrors* sts_c4_1sub[NV][NP];
		TGraphErrors* sts_c6_1sub[NV][NP];
		TGraphErrors* sts_nc4_1sub[NV][NP];
		TGraphErrors* sts_nc6_1sub[NV][NP];
		TGraphErrors* sts_sc_1sub[NV][NP];
		TGraphErrors* sts_nsc_1sub[NV][NP];
		TGraphErrors* sts_ac_1sub[NV][NP];
		TGraphErrors* sts_nac_1sub[NV][NP];
		TGraphErrors* sts_isGauss_1sub[NV][NP];
		TGraphErrors* sts_isPower_1sub[NV][NP];
		TGraphErrors* sts_c2_3sub[NV][NP];
		TGraphErrors* sts_c4_3sub[NV][NP];
		TGraphErrors* sts_nc4_3sub[NV][NP];
		TGraphErrors* sts_sc_3sub[NV][NP];
		TGraphErrors* sts_nsc_3sub[NV][NP];
		TGraphErrors* sts_ac_3sub[NV][NP];
		TGraphErrors* sts_nac_3sub[NV][NP];

		TGraphAsymmErrors* sys_c2_1sub[NV][NP];
		TGraphAsymmErrors* sys_c4_1sub[NV][NP];
		TGraphAsymmErrors* sys_c6_1sub[NV][NP];
		TGraphAsymmErrors* sys_nc4_1sub[NV][NP];
		TGraphAsymmErrors* sys_nc6_1sub[NV][NP];
		TGraphAsymmErrors* sys_sc_1sub[NV][NP];
		TGraphAsymmErrors* sys_nsc_1sub[NV][NP];
		TGraphAsymmErrors* sys_ac_1sub[NV][NP];
		TGraphAsymmErrors* sys_nac_1sub[NV][NP];
		TGraphAsymmErrors* sys_isGauss_1sub[NV][NP];
		TGraphAsymmErrors* sys_isPower_1sub[NV][NP];
		TGraphAsymmErrors* sys_c2_3sub[NV][NP];
		TGraphAsymmErrors* sys_c4_3sub[NV][NP];
		TGraphAsymmErrors* sys_nc4_3sub[NV][NP];
		TGraphAsymmErrors* sys_sc_3sub[NV][NP];
		TGraphAsymmErrors* sys_nsc_3sub[NV][NP];
		TGraphAsymmErrors* sys_ac_3sub[NV][NP];
		TGraphAsymmErrors* sys_nac_3sub[NV][NP];

		void readHist_VP(TFile*, TGraphErrors*&, const char*, const char*, int iV, int iP);
		void combine(TGraphErrors* gSts, TGraphErrors* gSysUp, TGraphErrors* gSysLw, TGraphAsymmErrors*& gSys, const char*, int iV, int iP);
		void trim(TGraphErrors* gSts, TGraphAsymmErrors* gSys);

	public:
		Phase4(unsigned int iBin);
		~Phase4();
		void initialize(unsigned int iBin);
		void execute();
		void finalize(unsigned int iBin);
};

#endif
