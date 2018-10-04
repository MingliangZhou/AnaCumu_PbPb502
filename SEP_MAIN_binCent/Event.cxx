#include "Event.h"

Event::Event(TChain* treeIn)
{
	treeIn->SetBranchAddress("RunNumber",  &m_runNo);
	treeIn->SetBranchAddress("lbn",        &m_lbn);
	treeIn->SetBranchAddress("bcid",       &m_bcid);
	treeIn->SetBranchAddress("m_PileUp",   &m_pileUp);
	treeIn->SetBranchAddress("m_TrgTp",    &m_TrgTp);
	treeIn->SetBranchAddress("Centrality", &m_centrality);
	treeIn->SetBranchAddress("Fcal_Et",    &m_fcalEt);
	treeIn->SetBranchAddress("vx_z",       &m_zPriVtx);
	treeIn->SetBranchAddress("trk_n",      &m_nTrk);
	treeIn->SetBranchAddress("trk_data",   m_trkInfo);
}

Event::~Event()
{
}

bool Event::isGoodEvent()
{
	if(m_pileUp) return false; // pileup rejection
	if(!(m_TrgTp&1) && !(m_TrgTp&2)) return false; // minbias trigger and ultra-central trigger
	if(!isGoodEvent_PbPb2015(m_runNo,m_lbn,m_bcid)) return false; // GRL
	if(fabs(m_zPriVtx)>=cutZvtx) return false; // zVtx cut

	return true;
}

void Event::setEvent(Tool* tool)
{
	m_evtWght = tool->detTrigWght(m_fcalEt);
	m_cent = tool->detCent(m_fcalEt);

	unsigned int infoTmp = 0;
	double Qx[NVm][NP][NW] = {{{0}}};       double Qy[NVm][NP][NW] = {{{0}}};
	int tag3sub = 0; double Qx3sub[NVm][NP][3][NW] = {{{{0}}}}; double Qy3sub[NVm][NP][3][NW] = {{{{0}}}};
	memset(S,0,sizeof S); memset(nTrk3sub,0,sizeof nTrk3sub);
	for(int iT=0; iT<m_nTrk; iT++)
	{
		infoTmp = m_trkInfo[iT];
		m_trkPhi[iT] = convertPhi(infoTmp & 0x1ff);   infoTmp >>= 9;
		m_trkEta[iT] = convertEta(infoTmp & 0x1ff);   infoTmp >>= 9;
		m_trkChg[iT] = infoTmp & 0x1;                        infoTmp >>= 1;
		m_trkPt[iT]  = convertPt(infoTmp & 0x3f);     infoTmp >>= 6;
		m_trkQual[iT]= infoTmp & 0x7;

		double trkEff = tool->detTrkEff(m_fcalEt,m_trkPt[iT],m_trkEta[iT]);
		double trkFak = tool->detTrkFak(m_fcalEt,m_trkPt[iT],m_trkEta[iT]);
		double trkFla = tool->detFlat(m_zPriVtx,m_centrality,m_trkChg[iT],m_trkPt[iT],m_trkEta[iT],m_trkPhi[iT]);
		m_trkWei[iT] = trkFla*(1-trkFak)/trkEff;

		if(fabs(m_trkEta[iT])>=cutEta) continue;
		if((m_trkQual[iT] & int(pow(2,cutQual)))==0) continue;

		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iW=0; iW<NW; iW++)
			{
				if(m_trkPt[iT]<minPt[iP] || m_trkPt[iT]>=maxPt[iP]) continue;
				S[iP][iW] += pow(m_trkWei[iT],iW);
				tag3sub = detSubEvt(3,m_trkEta[iT]);
				if(tag3sub>=0) nTrk3sub[iP][tag3sub][iW] += pow(m_trkWei[iT],iW);
				for(unsigned int iV=0; iV<NVm; iV++)
				{
					Qx[iV][iP][iW] += cos(iV*m_trkPhi[iT])*pow(m_trkWei[iT],iW);
					Qy[iV][iP][iW] += sin(iV*m_trkPhi[iT])*pow(m_trkWei[iT],iW);
					if(tag3sub>=0)
					{
						Qx3sub[iV][iP][tag3sub][iW] += cos(iV*m_trkPhi[iT])*pow(m_trkWei[iT],iW);
						Qy3sub[iV][iP][tag3sub][iW] += sin(iV*m_trkPhi[iT])*pow(m_trkWei[iT],iW);
					}
				}
			}
		}
	}

	for(unsigned int iV=0; iV<NVm; iV++)
	{
		for(unsigned int iW=0; iW<NW; iW++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				Q[iV][iP][iW] = TComplex(Qx[iV][iP][iW],Qy[iV][iP][iW]);
				for(int i=0; i<3; i++) Q3sub[iV][iP][i][iW] = TComplex(Qx3sub[iV][iP][i][iW],Qy3sub[iV][iP][i][iW]);
			}
		}
	}
}

int Event::detSubEvt(int numSubEvent, double eta)
{
	int tag = -1;
	switch(numSubEvent)
	{
		case 1:
			{
				tag = 0;
				return tag;
			}
		case 2:
			{
				if(eta<0-cutEtaGap/2) tag = 0;
				else if(eta>=0+cutEtaGap/2) tag = 1;
				else tag = -1;
				return tag;
			}
		case 3:
			{
				if(fabs(eta)<cutEta/3-cutEtaGap/2) tag = 0;
				else if(eta< 0-cutEta/3.-cutEtaGap/2.) tag = 1;
				else if(eta>=0+cutEta/3.+cutEtaGap/2.) tag = 2;
				else tag = -1;
				return tag;
			}
		case 4:
			{
				if(eta<-cutEta/2.-cutEtaGap/2.) tag = 0;
				else if(eta>=-cutEta/2.+cutEtaGap/2. && eta<0-cutEtaGap/2.) tag = 1;
				else if(eta>=0+cutEtaGap/2. && eta<+cutEta/2.-cutEtaGap/2.) tag = 2;
				else if(eta>=+cutEta/2.+cutEtaGap/2.) tag = 3;
				else tag = -1;
				return tag;
			}
		default:
			{
				std::cout<<"Wrong number of subevents!"<<std::endl;
				return -1;
			}
	}
}

