#ifndef MultiCorr_H_
#define MultiCorr_H_

#include "Event.h"

class MultiCorr
{
	private:
		void initHist_VP(TH1D*&, const char*, int iV, int iP); // initialize multi-particle correlation for standard
		void initHist_AVP(TH1D*&, const char*, int iA, int iV, int iP); // initialize multi-particle correlation for 1-subevent

		double cal_2pc_1sub(Event* evt, int iV, int iP); // 2-particle correlation in standard
		double cal_w2pc_1sub(Event* evt, int iP); // weight for 2-particle correlation in standard
		double cal_4pc_1sub(Event* evt, int iV, int iP); // 4-particle correlation in standard
		double cal_w4pc_1sub(Event* evt, int iP); // weight for 4-particle correlation in standard
		double cal_6pc_1sub(Event* evt, int iV, int iP); // 6-particle correlation in standard
		double cal_w6pc_1sub(Event* evt, int iP); // weight for 6-particle correlation in standard
		double cal_4psc_1sub(Event* evt, int iV, int iP); // 4-particle symmetric correlation 1+2-1-2, 1+3-1-3, 2+3-2-3, 2+4-2,4
		double cal_3pac_1sub(Event* evt, int iV, int iP); // 3-particle asymmetric correlation 1+1-2, 1+2-3, 2+2-4
		double cal_w3pc_1sub(Event* evt, int iP); // weight for 3-particle asymmetric correlation
		
		double cal_2pc_1sub_BG(Event* evt1, Event* evt2, int iV, int iP); // background
		double cal_w2pc_1sub_BG(Event* evt1, Event* evt2, int iP); // weight for background
		double cal_4pc_1sub_BG(Event* evt1, Event* evt2, Event* evt3, Event* evt4, int iV, int iP); // background
		double cal_w4pc_1sub_BG(Event* evt1, Event* evt2, Event* evt3, Event* evt4, int iP); // weight for background

		double cal_2pc_3sub(Event* evt, int iV, int iP, int iO); // 3-subevent: 2 options
		double cal_w2pc_3sub(Event* evt, int iP, int iO); // weight for 3-subevent: 2 options
		double cal_4pc_3sub(Event* evt, int iV, int iP); // 3-subevent
		double cal_w4pc_3sub(Event* evt, int iP); // weight for 3-subevent
		double cal_4psc_3sub(Event* evt, int iV, int iP); // 4-particle symmetric correlation
		double cal_3pac_3sub(Event* evt, int iV, int iP); // 3-subevent
		double cal_w3pc_3sub(Event* evt, int iV); // weight for 3-subevent

		void swap(double&, double&, double&); // swap double for 3-subevent
		void swap(TComplex&, TComplex&, TComplex&); // swap complex for 3-subevent
		void swap(Event*& evt1, Event*& evt2, Event*& evt3); // swap event for background

	public:
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

		MultiCorr();
		~MultiCorr();
		void fill_1sub(Event*);
		void fill_1sub_BG(vector<Event*>);
		void fill_3sub(Event*);
		void writeHist(TFile*&);
};

#endif
