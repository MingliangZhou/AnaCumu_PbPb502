#ifndef Cumu_H_
#define Cumu_H_

#include "Rule.h"

class Cumu
{
	private:
		void run_c_1sub(); // run standard cumulant runculation
		void run_sc_1sub(); // run standard symmetric cumulant runculation
		void run_ac_1sub(); // run standard asymmetric cumulant runculation
		void run_c_3sub(); // run 3-subevent cumulant runculation
		void run_sc_3sub(); // run 3-subevent symmetric cumulant runculation
		void run_ac_3sub(); // run 3-subevent asymmetric cumulant runculation
		void run_nc_1sub(); // run standard normalized cumulant runculation
		void run_nsc_1sub(); // run standard normalized symmetric cumulant runculation
		void run_nac_1sub(); // run standard normalized asymmetric cumulant runculation
		void run_fluc_1sub(); // run standard cumulant flow fluctuation runculation
		void run_nc_3sub(); // run 3-subevent normalized cumulant runculation
		void run_nsc_3sub(); // run 3-subevent normalized symmetric cumulant runculation
		void run_nac_3sub(); // run 3-subevent normalized asymmetric cumulant runculation
		void cal_c_1sub(TH1D* pc2, TH1D* pc4, TH1D* pc6, TH1D* c2, TH1D* c4, TH1D* c6);
		void cal_sc(TH1D* pc2_1, TH1D* pc2_2, TH1D* psc4, TH1D* sc4);
		void cal_c_3sub(TH1D* pc2_1, TH1D* pc2_2, TH1D* pc4, TH1D* c2, TH1D* c4);
		void cal_nc(TH1D* c2, TH1D* c4, TH1D* c6, TH1D* nc4, TH1D* nc6);
		void cal_nsc(TH1D* c2_1, TH1D* c2_2, TH1D* sc4, TH1D* nsc4);
		void cal_nac(TH1D* c4, TH1D* c2, TH1D* ac4, TH1D* nac4);
		void cal_fluc(TH1D* c2, TH1D* c4, TH1D* c6, TH1D* isGauss, TH1D* isPower);
		void rebin_cumu(TH1D* hIn, TH1D* hWght, TH1D* hOut);
		void merge3sub(TH1D* h1, TH1D* h1w, TH1D* h2, TH1D* h2w, TH1D* h3, TH1D* h3w);
		void readHist_VP(TFile*, TH1D*&, const char*, int iV, int iP);
		void readHist_AVP(TFile*, TH1D*&, const char*, int iA, int iV, int iP);
		void initHist_VP(TH1D*&, const char*, int iV, int iP, int nBins);
		void initHist_AVP(TH1D*&, const char*, int iA, int iV, int iP, int nBins);

	public:
		// input
		TH1D* cnt_1sub[NP]; // number of events passing each pT cuts in standard
		TH1D* pc2_1sub_mean[NV][NP]; // mean of <corr2>_{n|n}
		TH1D* pc2_1sub_wght[NV][NP]; // weight of <corr2>_{n|n}
		TH1D* pc4_1sub_mean[NV][NP]; // mean of <corr4>_{n|n}
		TH1D* pc4_1sub_wght[NV][NP]; // weight of <corr4>_{n|n}
		TH1D* pc6_1sub_mean[NV][NP]; // mean of <corr6>_{n|n}
		TH1D* pc6_1sub_wght[NV][NP]; // weight of <corr6>_{n|n}
		TH1D* psc4_1sub_mean[NV][NP];
		TH1D* psc4_1sub_wght[NV][NP];
		TH1D* pac3_1sub_mean[NV][NP];
		TH1D* pac3_1sub_wght[NV][NP];

		TH1D* pc2_1_1sub_BG_mean[NA][NV][NP];
		TH1D* pc2_1_1sub_BG_wght[NA][NV][NP];
		TH1D* pc2_2_1sub_BG_mean[NA][NV][NP];
		TH1D* pc2_2_1sub_BG_wght[NA][NV][NP];
		TH1D* pc2_3_1sub_BG_mean[NA][NV][NP];
		TH1D* pc2_3_1sub_BG_wght[NA][NV][NP];
		TH1D* pc2_4_1sub_BG_mean[NA][NV][NP];
		TH1D* pc2_4_1sub_BG_wght[NA][NV][NP];
		TH1D* pc4_1sub_BG_mean[NA][NV][NP];
		TH1D* pc4_1sub_BG_wght[NA][NV][NP];

		TH1D* cnt_3sub[NA][NP]; // number of events passing each pT cuts in 3-subevent
		TH1D* pc2_1_3sub_mean[NA][NV][NP]; // mean of <2>_{a|b}
		TH1D* pc2_1_3sub_wght[NA][NV][NP]; // weight for <2>_{a|b}
		TH1D* pc2_2_3sub_mean[NA][NV][NP]; // mean of <2>_{a|c}
		TH1D* pc2_2_3sub_wght[NA][NV][NP]; // weight for <2>_{a|c}
		TH1D* pc4_3sub_mean[NA][NV][NP]; // mean of <4>_{a,a|b,c}
		TH1D* pc4_3sub_wght[NA][NV][NP]; // weight of <4>_{a,a|b,c}
		TH1D* psc4_3sub_mean[NA][NV][NP];
		TH1D* psc4_3sub_wght[NA][NV][NP];
		TH1D* pac3_3sub_mean[NA][NV][NP];
		TH1D* pac3_3sub_wght[NA][NV][NP];

		// output
		TH1D* raw_cnt_1sub[NP];
		TH1D* rbn_cnt_1sub[NP];
		TH1D* ori_c2_1sub[NV][NP]; // before rebin
		TH1D* raw_c2_1sub[NV][NP]; // rebin first time
		TH1D* rbn_c2_1sub[NV][NP]; // rebin second time
		TH1D* ori_c4_1sub[NV][NP];
		TH1D* raw_c4_1sub[NV][NP];
		TH1D* rbn_c4_1sub[NV][NP];
		TH1D* ori_c6_1sub[NV][NP];
		TH1D* raw_c6_1sub[NV][NP];
		TH1D* rbn_c6_1sub[NV][NP];
		TH1D* raw_nc4_1sub[NV][NP];
		TH1D* rbn_nc4_1sub[NV][NP];
		TH1D* raw_nc6_1sub[NV][NP];
		TH1D* rbn_nc6_1sub[NV][NP];
		TH1D* ori_sc_1sub[NV][NP];
		TH1D* raw_sc_1sub[NV][NP];
		TH1D* rbn_sc_1sub[NV][NP];
		TH1D* raw_nsc_1sub[NV][NP];
		TH1D* rbn_nsc_1sub[NV][NP];
		TH1D* ori_ac_1sub[NV][NP];
		TH1D* raw_ac_1sub[NV][NP];
		TH1D* rbn_ac_1sub[NV][NP];
		TH1D* raw_nac_1sub[NV][NP];
		TH1D* rbn_nac_1sub[NV][NP];
		TH1D* raw_isGauss_1sub[NV][NP];
		TH1D* rbn_isGauss_1sub[NV][NP];
		TH1D* raw_isPower_1sub[NV][NP];
		TH1D* rbn_isPower_1sub[NV][NP];

		TH1D* raw_cnt_3sub[NA][NP];
		TH1D* rbn_cnt_3sub[NA][NP];
		TH1D* ori_c2_3sub[NA][NV][NP];
		TH1D* raw_c2_3sub[NA][NV][NP];
		TH1D* rbn_c2_3sub[NA][NV][NP];
		TH1D* ori_c4_3sub[NA][NV][NP];
		TH1D* raw_c4_3sub[NA][NV][NP];
		TH1D* raw_p4_3sub[NA][NV][NP]; // special for normalized asymmetric cumulant
		TH1D* rbn_c4_3sub[NA][NV][NP];
		TH1D* raw_nc4_3sub[NA][NV][NP];
		TH1D* rbn_nc4_3sub[NA][NV][NP];
		TH1D* raw_nc6_3sub[NA][NV][NP];
		TH1D* rbn_nc6_3sub[NA][NV][NP];
		TH1D* ori_sc_3sub[NA][NV][NP];
		TH1D* raw_sc_3sub[NA][NV][NP];
		TH1D* rbn_sc_3sub[NA][NV][NP];
		TH1D* raw_nsc_3sub[NA][NV][NP];
		TH1D* rbn_nsc_3sub[NA][NV][NP];
		TH1D* ori_ac_3sub[NA][NV][NP];
		TH1D* raw_ac_3sub[NA][NV][NP];
		TH1D* rbn_ac_3sub[NA][NV][NP];
		TH1D* raw_nac_3sub[NA][NV][NP];
		TH1D* rbn_nac_3sub[NA][NV][NP];
		
		Cumu(TFile*, unsigned int iRebin);
		~Cumu();
		void cal_all();
		void writeHist(TFile*&);
};

#endif
