#ifndef Event_H_
#define Event_H_

#include "Tool.h"
#include "INPUT/isGoodEvent_PbPb2015.C"

class Event
{
	private:
		unsigned int m_runNo; // run number
		unsigned int m_lbn; // lumi-block number
		unsigned int m_bcid; // bunch crossing number
		bool m_pileUp; // is pileup collision?
		char m_TrgTp; // trigger info
		int m_centrality; // centrality
		float m_fcalEt; // FCal E_T from both sides [TeV]
		float m_zPriVtx; // z position of vertex
		int m_nTrk; // number of tracks
		unsigned int m_trkInfo[maxNtrk]; // all tracking info

		double m_evtWght;
		double m_cent;
		double m_trkWei[maxNtrk]; // track weight
		double m_trkPt[maxNtrk]; // track pT
		double m_trkEta[maxNtrk]; // track eta
		double m_trkPhi[maxNtrk]; // track phi
		unsigned int m_trkQual[maxNtrk]; // track quality: &2 loose, &4 tight
		unsigned int m_trkChg[maxNtrk]; // track charge

		double convertPt(int ptbin) // convert pT to double
		{
			float ptgev=110;
			if(ptbin<20) ptgev = 0.1*ptbin+0.05;
			else if(ptbin<30) ptgev = 0.2*(ptbin-20)+2.1;
			else if(ptbin<38) ptgev = 0.5*(ptbin-30)+4.2;
			else if(ptbin<54) ptgev = 2*(ptbin-38)+8.8;
			else if(ptbin<60) ptgev = 5*(ptbin-54)+42;
			else if(ptbin<63) ptgev = 10*(ptbin-60)+74;
			return ptgev;
		}
		double convertEta(int charEta) { return double((charEta+0.5)*cutEta*2/500-cutEta); } // convert eta to double
		double convertPhi(int charPhi) { return double((charPhi+0.5)*2*TMath::Pi()/(0x1ff+1)); } // convert phi to double
		int detSubEvt(int, double); // determine the subevent tag: B, A, C

	public:
		double S[NP][NW];     TComplex Q[NVm][NP][NW]; // RFP
		double nTrk3sub[NP][3][NW]; TComplex Q3sub[NVm][NP][3][NW]; // harmonic | subevt | power of weight

		Event(TChain*);
		~Event();
		bool isGoodEvent();
		void setEvent(Tool*);
		double evtZvtx() const { return m_zPriVtx; }
		double evtCls() const { return m_cent; }
		double evtWght() const { return m_evtWght; }
		bool isMB() const { return m_TrgTp&1; }
		bool isHMT() const { return m_TrgTp&2; }
		double evtFCal() const {return m_fcalEt; }
		double evtCent() const {return m_cent; }
		double evtNchRec() const {return S[0][0]; }
		double evtNch() const {return S[0][1]; }
};

#endif
