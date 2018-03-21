#ifndef Event_H_
#define Event_H_

#include "Tool.h"
#include "../../MAIN_binCent/INPUT/isGoodEvent_PbPb2015.C"

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
		float m_fcalEtA; // FCal E_T from both sides [TeV]
		float m_fcalEtC; // FCal E_T from both sides [TeV]
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

		double convertPt(int charPt) { return double((charPt+0.5)/10.); }
		double convertEta(int charEta) { return double((charEta+0.5)/250*2*cutEta-cutEta); }
		double convertPhi(int charPhi) { return double((charPhi+0.5)/pow(2,14)*2*TMath::Pi()-TMath::Pi()); }
		int detSubEvt(int, double); // determine the subevent tag: B, A, C

	public:
		double S[NP][NW];     TComplex Q[NVm][NP][NW]; // RFP
		double nTrk3sub[NP][3][NW]; TComplex Q3sub[NVm][NP][3][NW]; // harmonic | subevt | power of weight

		Event(TChain*);
		~Event();
		bool isGoodEvent();
		void setEvent(Tool*);
		double evtZvtx() const { return m_zPriVtx; }
		double evtCls() const { return m_centrality; }
		double evtWght() const { return m_evtWght; }
};

#endif
