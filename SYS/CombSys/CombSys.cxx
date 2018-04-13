#include "CombSys.h"

CombSys::CombSys(unsigned int iBin)
{
	initialize(iBin);
	execute();
	finalize(iBin);
}

CombSys::~CombSys()
{
}

void CombSys::execute()
{
	cout<<"execute..."<<endl;

	for(unsigned int iS=0; iS<NS; iS++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				cal_ratio(c2_1sub[0][iS][iV][iP],c2_1sub[1][iS][iV][iP],ratio_c2_1sub[iS][iV][iP],"ratio_c2_1sub",iS,iV,iP);
				cal_ratio(c4_1sub[0][iS][iV][iP],c4_1sub[1][iS][iV][iP],ratio_c4_1sub[iS][iV][iP],"ratio_c4_1sub",iS,iV,iP);
				cal_ratio(c6_1sub[0][iS][iV][iP],c6_1sub[1][iS][iV][iP],ratio_c6_1sub[iS][iV][iP],"ratio_c6_1sub",iS,iV,iP);
				cal_ratio(nc4_1sub[0][iS][iV][iP],nc4_1sub[1][iS][iV][iP],ratio_nc4_1sub[iS][iV][iP],"ratio_nc4_1sub",iS,iV,iP);
				cal_ratio(nc6_1sub[0][iS][iV][iP],nc6_1sub[1][iS][iV][iP],ratio_nc6_1sub[iS][iV][iP],"ratio_nc6_1sub",iS,iV,iP);
				cal_ratio(sc_1sub[0][iS][iV][iP],sc_1sub[1][iS][iV][iP],ratio_sc_1sub[iS][iV][iP],"ratio_sc_1sub",iS,iV,iP);
				cal_ratio(nsc_1sub[0][iS][iV][iP],nsc_1sub[1][iS][iV][iP],ratio_nsc_1sub[iS][iV][iP],"ratio_nsc_1sub",iS,iV,iP);
				cal_ratio(ac_1sub[0][iS][iV][iP],ac_1sub[1][iS][iV][iP],ratio_ac_1sub[iS][iV][iP],"ratio_ac_1sub",iS,iV,iP);
				cal_ratio(nac_1sub[0][iS][iV][iP],nac_1sub[1][iS][iV][iP],ratio_nac_1sub[iS][iV][iP],"ratio_nac_1sub",iS,iV,iP);
				cal_ratio(isGauss_1sub[0][iS][iV][iP],isGauss_1sub[1][iS][iV][iP],ratio_isGauss_1sub[iS][iV][iP],"ratio_isGauss_1sub",iS,iV,iP);
				cal_ratio(isPower_1sub[0][iS][iV][iP],isPower_1sub[1][iS][iV][iP],ratio_isPower_1sub[iS][iV][iP],"ratio_isPower_1sub",iS,iV,iP);
				cal_ratio(c2_3sub[0][iS][iV][iP],c2_3sub[1][iS][iV][iP],ratio_c2_3sub[iS][iV][iP],"ratio_c2_3sub",iS,iV,iP);
				cal_ratio(c4_3sub[0][iS][iV][iP],c4_3sub[1][iS][iV][iP],ratio_c4_3sub[iS][iV][iP],"ratio_c4_3sub",iS,iV,iP);
				cal_ratio(nc4_3sub[0][iS][iV][iP],nc4_3sub[1][iS][iV][iP],ratio_nc4_3sub[iS][iV][iP],"ratio_nc4_3sub",iS,iV,iP);
				cal_ratio(sc_3sub[0][iS][iV][iP],sc_3sub[1][iS][iV][iP],ratio_sc_3sub[iS][iV][iP],"ratio_sc_3sub",iS,iV,iP);
				cal_ratio(nsc_3sub[0][iS][iV][iP],nsc_3sub[1][iS][iV][iP],ratio_nsc_3sub[iS][iV][iP],"ratio_nsc_3sub",iS,iV,iP);
				cal_ratio(ac_3sub[0][iS][iV][iP],ac_3sub[1][iS][iV][iP],ratio_ac_3sub[iS][iV][iP],"ratio_ac_3sub",iS,iV,iP);
				cal_ratio(nac_3sub[0][iS][iV][iP],nac_3sub[1][iS][iV][iP],ratio_nac_3sub[iS][iV][iP],"ratio_nac_3sub",iS,iV,iP);
			}
		}
	}

	smooth();

	vector<TGraphErrors*> gVec;
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c2_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_c2_1sub[iV][iP],sysLw_c2_1sub[iV][iP],"c2_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c4_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_c4_1sub[iV][iP],sysLw_c4_1sub[iV][iP],"c4_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c6_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_c6_1sub[iV][iP],sysLw_c6_1sub[iV][iP],"c6_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nc4_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nc4_1sub[iV][iP],sysLw_nc4_1sub[iV][iP],"nc4_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nc6_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nc6_1sub[iV][iP],sysLw_nc6_1sub[iV][iP],"nc6_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_sc_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_sc_1sub[iV][iP],sysLw_sc_1sub[iV][iP],"sc_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nsc_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nsc_1sub[iV][iP],sysLw_nsc_1sub[iV][iP],"nsc_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_ac_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_ac_1sub[iV][iP],sysLw_ac_1sub[iV][iP],"ac_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nac_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nac_1sub[iV][iP],sysLw_nac_1sub[iV][iP],"nac_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_isGauss_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_isGauss_1sub[iV][iP],sysLw_isGauss_1sub[iV][iP],"isGauss_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_isPower_1sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_isPower_1sub[iV][iP],sysLw_isPower_1sub[iV][iP],"isPower_1sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c2_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_c2_3sub[iV][iP],sysLw_c2_3sub[iV][iP],"c2_3sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c4_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_c4_3sub[iV][iP],sysLw_c4_3sub[iV][iP],"c4_3sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nc4_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nc4_3sub[iV][iP],sysLw_nc4_3sub[iV][iP],"nc4_3sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_sc_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_sc_3sub[iV][iP],sysLw_sc_3sub[iV][iP],"sc_3sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nsc_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nsc_3sub[iV][iP],sysLw_nsc_3sub[iV][iP],"nsc_3sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_ac_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_ac_3sub[iV][iP],sysLw_ac_3sub[iV][iP],"ac_3sub",iV,iP); gVec.clear();
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nac_3sub[iS][iV][iP]);
			cal_comb(gVec,sysUp_nac_3sub[iV][iP],sysLw_nac_3sub[iV][iP],"nac_3sub",iV,iP); gVec.clear();
		}
	}
}

void CombSys::initialize(unsigned int iBin)
{
	cout<<"initialize..."<<endl;

	TFile* fIn[2][NS];
	for(unsigned int iF=0; iF<2; iF++)
	{
		for(unsigned int iS=0; iS<NS; iS++)
		{
			if(iF==0 && iS==0) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==0 && iS==1) sprintf(name,"../../trkEffLw/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==0 && iS==2) sprintf(name,"../../trkEffUp/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==0 && iS==3) sprintf(name,"../../trkSel/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==0 && iS==4) sprintf(name,"../../pileup/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==0 && iS==5) sprintf(name,"../../mcTruth/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==0 && iS==6) sprintf(name,"../../flat/OUTPUT/Phase3/Phase3_bin%d.root",iBin);

			if(iF==1 && iS==0) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==1 && iS==1) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==1 && iS==2) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==1 && iS==3) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==1 && iS==4) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==1 && iS==5) sprintf(name,"../../mcRecon/OUTPUT/Phase3/Phase3_bin%d.root",iBin);
			if(iF==1 && iS==6) sprintf(name,"../../default/OUTPUT/Phase3/Phase3_bin%d.root",iBin);

			fIn[iF][iS] = new TFile(name,"READ");
			for(unsigned int iV=0; iV<NV; iV++)
			{
				for(unsigned int iP=0; iP<NP; iP++)
				{
					readHist_VP(fIn[iF][iS],c2_1sub[iF][iS][iV][iP],"c2_1sub","c2_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],c4_1sub[iF][iS][iV][iP],"c4_1sub","c4_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],c6_1sub[iF][iS][iV][iP],"c6_1sub","c6_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nc4_1sub[iF][iS][iV][iP],"nc4_1sub","nc4_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nc6_1sub[iF][iS][iV][iP],"nc6_1sub","nc6_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],sc_1sub[iF][iS][iV][iP],"sc_1sub","sc_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nsc_1sub[iF][iS][iV][iP],"nsc_1sub","nsc_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],ac_1sub[iF][iS][iV][iP],"ac_1sub","ac_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nac_1sub[iF][iS][iV][iP],"nac_1sub","nac_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],isGauss_1sub[iF][iS][iV][iP],"isGauss_1sub","isGauss_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],isPower_1sub[iF][iS][iV][iP],"isPower_1sub","isPower_1sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],c2_3sub[iF][iS][iV][iP],"c2_3sub","c2_3sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],c4_3sub[iF][iS][iV][iP],"c4_3sub","c4_3sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nc4_3sub[iF][iS][iV][iP],"nc4_3sub","nc4_3sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],sc_3sub[iF][iS][iV][iP],"sc_3sub","sc_3sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nsc_3sub[iF][iS][iV][iP],"nsc_3sub","nsc_3sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],ac_3sub[iF][iS][iV][iP],"ac_3sub","ac_3sub",iF,iS,iV,iP);
					readHist_VP(fIn[iF][iS],nac_3sub[iF][iS][iV][iP],"nac_3sub","nac_3sub",iF,iS,iV,iP);
				}
			}
		}
	}
}

void CombSys::finalize(unsigned int iBin)
{
	cout<<"finalize..."<<endl;
	
	sprintf(name,"../OUTPUT/Sys_bin%d.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();
	for(unsigned int iF=0; iF<2; iF++)
	{
		for(unsigned int iS=0; iS<NS; iS++)
		{
			for(unsigned int iV=0; iV<NV; iV++)
			{
				for(unsigned int iP=0; iP<NP; iP++)
				{
					c2_1sub[iF][iS][iV][iP]->Write();
					c4_1sub[iF][iS][iV][iP]->Write();
					c6_1sub[iF][iS][iV][iP]->Write();
					nc4_1sub[iF][iS][iV][iP]->Write();
					nc6_1sub[iF][iS][iV][iP]->Write();
					sc_1sub[iF][iS][iV][iP]->Write();
					nsc_1sub[iF][iS][iV][iP]->Write();
					ac_1sub[iF][iS][iV][iP]->Write();
					nac_1sub[iF][iS][iV][iP]->Write();
					isGauss_1sub[iF][iS][iV][iP]->Write();
					isPower_1sub[iF][iS][iV][iP]->Write();
					c2_3sub[iF][iS][iV][iP]->Write();
					c4_3sub[iF][iS][iV][iP]->Write();
					nc4_3sub[iF][iS][iV][iP]->Write();
					sc_3sub[iF][iS][iV][iP]->Write();
					nsc_3sub[iF][iS][iV][iP]->Write();
					ac_3sub[iF][iS][iV][iP]->Write();
					nac_3sub[iF][iS][iV][iP]->Write();
				}
			}
		}
	}
	for(unsigned int iS=0; iS<NS; iS++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				ratio_c2_1sub[iS][iV][iP]->Write();
				ratio_c4_1sub[iS][iV][iP]->Write();
				ratio_c6_1sub[iS][iV][iP]->Write();
				ratio_nc4_1sub[iS][iV][iP]->Write();
				ratio_nc6_1sub[iS][iV][iP]->Write();
				ratio_sc_1sub[iS][iV][iP]->Write();
				ratio_nsc_1sub[iS][iV][iP]->Write();
				ratio_ac_1sub[iS][iV][iP]->Write();
				ratio_nac_1sub[iS][iV][iP]->Write();
				ratio_isGauss_1sub[iS][iV][iP]->Write();
				ratio_isPower_1sub[iS][iV][iP]->Write();
				ratio_c2_3sub[iS][iV][iP]->Write();
				ratio_c4_3sub[iS][iV][iP]->Write();
				ratio_nc4_3sub[iS][iV][iP]->Write();
				ratio_sc_3sub[iS][iV][iP]->Write();
				ratio_nsc_3sub[iS][iV][iP]->Write();
				ratio_ac_3sub[iS][iV][iP]->Write();
				ratio_nac_3sub[iS][iV][iP]->Write();
			}
		}
	}
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			sysUp_c2_1sub[iV][iP]->Write();
			sysUp_c4_1sub[iV][iP]->Write();
			sysUp_c6_1sub[iV][iP]->Write();
			sysUp_nc4_1sub[iV][iP]->Write();
			sysUp_nc6_1sub[iV][iP]->Write();
			sysUp_sc_1sub[iV][iP]->Write();
			sysUp_nsc_1sub[iV][iP]->Write();
			sysUp_ac_1sub[iV][iP]->Write();
			sysUp_nac_1sub[iV][iP]->Write();
			sysUp_isGauss_1sub[iV][iP]->Write();
			sysUp_isPower_1sub[iV][iP]->Write();
			sysUp_c2_3sub[iV][iP]->Write();
			sysUp_c4_3sub[iV][iP]->Write();
			sysUp_nc4_3sub[iV][iP]->Write();
			sysUp_sc_3sub[iV][iP]->Write();
			sysUp_nsc_3sub[iV][iP]->Write();
			sysUp_ac_3sub[iV][iP]->Write();
			sysUp_nac_3sub[iV][iP]->Write();
		}
	}
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			sysLw_c2_1sub[iV][iP]->Write();
			sysLw_c4_1sub[iV][iP]->Write();
			sysLw_c6_1sub[iV][iP]->Write();
			sysLw_nc4_1sub[iV][iP]->Write();
			sysLw_nc6_1sub[iV][iP]->Write();
			sysLw_sc_1sub[iV][iP]->Write();
			sysLw_nsc_1sub[iV][iP]->Write();
			sysLw_ac_1sub[iV][iP]->Write();
			sysLw_nac_1sub[iV][iP]->Write();
			sysLw_isGauss_1sub[iV][iP]->Write();
			sysLw_isPower_1sub[iV][iP]->Write();
			sysLw_c2_3sub[iV][iP]->Write();
			sysLw_c4_3sub[iV][iP]->Write();
			sysLw_nc4_3sub[iV][iP]->Write();
			sysLw_sc_3sub[iV][iP]->Write();
			sysLw_nsc_3sub[iV][iP]->Write();
			sysLw_ac_3sub[iV][iP]->Write();
			sysLw_nac_3sub[iV][iP]->Write();
		}
	}
	fOut->Close();
}

void CombSys::cal_ratio(TGraphErrors* gChk, TGraphErrors* gDef, TGraphErrors*& gRatio, const char* gName, int iS, int iV, int iP)
{
	int NB = gDef->GetN();
	double x[NB];
	double y[NB];
	double xErr[NB];
	double yErr[NB];
	for(int iB=0; iB<NB; iB++)
	{
		double yDef = 0;
		double yChk = 0;
		double eDef = 0;
		double eChk = 0;
		gChk->GetPoint(iB,x[iB],yChk);
		gDef->GetPoint(iB,x[iB],yDef);
		eChk = gChk->GetErrorY(iB);
		eDef = gDef->GetErrorY(iB);

		xErr[iB] = 0;
		if(yDef!=0) y[iB] = yChk/yDef - 1.;
		else y[iB] = 0;
		if(yDef!=0 && yChk!=0) yErr[iB] = fabs(yChk/yDef*sqrt(pow(eChk/yChk,2)+pow(eDef/yDef,2)));
		else yErr[iB] = 0;
	}
	gRatio = new TGraphErrors(NB,x,y,xErr,yErr);
	sprintf(name,"%s_Sys%d_Har%d_Pt%d",gName,iS,iV,iP);
	gRatio->SetName(name);
}

void CombSys::smooth_ratio(TGraphErrors* gRatio, int iB, double percent)
{
	double x; double y;
	gRatio->GetPoint(iB,x,y);
	double yErr = gRatio->GetErrorY(iB);
	if(y>0)
	{
		if(fabs(y-percent*yErr)<yErr) y = percent*yErr;
		else y -= yErr;
	}
	else
	{
		if(fabs(y-percent*yErr)<yErr) y = -percent*yErr;
		else y += yErr;
	}
	gRatio->SetPoint(iB,x,y);
}

void CombSys::reduce_ratio(TGraphErrors* gRatio, int iB, double percent)
{
	double x; double y;
	gRatio->GetPoint(iB,x,y);
	y *= percent;
	gRatio->SetPoint(iB,x,y);
}

void CombSys::cal_comb(vector<TGraphErrors*> gVec, TGraphErrors*& gSysUp, TGraphErrors*& gSysLw, const char* gName, int iV, int iP)
{
	int NB = gVec.at(0)->GetN();
	double x[NB];
	double xErr[NB];
	double yUp[NB];
	double yLw[NB];
	double yUpErr[NB];
	double yLwErr[NB];
	for(int iB=0; iB<NB; iB++)
	{
		gVec.at(0)->GetPoint(iB,x[iB],yUp[iB]);
		xErr[iB] = 0;
		yUp[iB] = 0;
		yLw[iB] = 0;
		yUpErr[iB] = 0;
		yLwErr[iB] = 0;
		for(unsigned int iS=1; iS<NS; iS++)
		{
			double yTmp; double eTmp;
			gVec.at(iS)->GetPoint(iB,eTmp,yTmp);
			eTmp = gVec.at(iS)->GetErrorY(iB);
			if(yTmp>=0)
			{
				yUp[iB] += pow(yTmp,2);
				yUpErr[iB] += pow(eTmp,2);
			}
			else
			{
				yLw[iB] += pow(yTmp,2);
				yLwErr[iB] += pow(eTmp,2);
			}
		}
		yUp[iB] = +sqrt(yUp[iB]);
		yLw[iB] = -sqrt(yLw[iB]);
		yUpErr[iB] = sqrt(yUpErr[iB]);
		yLwErr[iB] = sqrt(yLwErr[iB]);
	}

	gSysUp = new TGraphErrors(NB,x,yUp,xErr,yUpErr);
	sprintf(name,"sysUp_%s_Har%d_Pt%d",gName,iV,iP);
	gSysUp->SetName(name);

	gSysLw = new TGraphErrors(NB,x,yLw,xErr,yLwErr);
	sprintf(name,"sysLw_%s_Har%d_Pt%d",gName,iV,iP);
	gSysLw->SetName(name);
}

void CombSys::readHist_VP(TFile* fIn, TGraphErrors*& hIn, const char* hName, const char* hNewName, int iF, int iS, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
	sprintf(name,"%s_File%d_Sys%d_Har%d_Pt%d",hNewName,iF,iS,iV,iP);
	hIn->SetName(name);
}

void CombSys::smooth()
{
	// sys_c2_1sub_Har4_Pt0
	reduce_ratio(ratio_c2_1sub[2][4][0],7,0.5);
	reduce_ratio(ratio_c2_1sub[3][4][0],7,0.5);

	// sys_c2_1sub_Har4_Pt1
	reduce_ratio(ratio_c2_1sub[2][4][1],7,0.5);
	reduce_ratio(ratio_c2_1sub[3][4][1],7,0.5);

	// sys_c4_1sub_Har1_Pt0
	smooth_ratio(ratio_c4_1sub[6][1][0],0,1);
	smooth_ratio(ratio_c4_1sub[6][1][0],1,0.5);
	smooth_ratio(ratio_c4_1sub[6][1][0],2,0.5);
	
	// sys_c4_1sub_Har1_Pt1
	smooth_ratio(ratio_c4_1sub[3][1][1],0,0.1);
	smooth_ratio(ratio_c4_1sub[6][1][1],2,0.2);
	reduce_ratio(ratio_c4_1sub[3][1][1],0,0.3);

	// sys_c4_1sub_Har3_Pt0
	smooth_ratio(ratio_c4_1sub[3][3][0],6,0.2);
	smooth_ratio(ratio_c4_1sub[6][3][0],7,0.1);

	// sys_c4_1sub_Har3_Pt1
	smooth_ratio(ratio_c4_1sub[3][3][1],6,0.1);
	smooth_ratio(ratio_c4_1sub[6][3][1],7,0.1);

	// sys_c4_1sub_Har4_Pt0
	smooth_ratio(ratio_c4_1sub[2][4][0],7,0.1);
	smooth_ratio(ratio_c4_1sub[3][4][0],7,0.1);

	// sys_c4_1sub_Har4_Pt1
	smooth_ratio(ratio_c4_1sub[2][4][1],7,0.2);
	smooth_ratio(ratio_c4_1sub[3][4][1],7,0.2);

	// sys_c4_3sub_Har1_Pt0
	smooth_ratio(ratio_c4_3sub[6][1][0],0,0.001);
	reduce_ratio(ratio_c4_3sub[6][1][0],0,0.5);
	reduce_ratio(ratio_c4_3sub[3][1][0],0,0.2);

	// sys_c4_3sub_Har1_Pt1
	smooth_ratio(ratio_c4_3sub[6][1][1],3,0.5);
	reduce_ratio(ratio_c4_3sub[6][1][1],3,0.3);
	reduce_ratio(ratio_c4_3sub[3][1][1],3,0.3);
	reduce_ratio(ratio_c4_3sub[4][1][1],5,0.3);
	reduce_ratio(ratio_c4_3sub[4][1][1],6,0.3);

	// sys_c4_3sub_Har3_Pt0
	smooth_ratio(ratio_c4_3sub[3][3][0],7,0.2);
	smooth_ratio(ratio_c4_3sub[6][3][0],7,0.2);

	// sys_c4_3sub_Har3_Pt1
	smooth_ratio(ratio_c4_3sub[4][3][1],6,0.1);
	smooth_ratio(ratio_c4_3sub[3][3][1],7,0.05);

	// sys_c4_3sub_Har4_Pt0
	smooth_ratio(ratio_c4_3sub[6][4][0],2,0.1);
	smooth_ratio(ratio_c4_3sub[3][4][0],2,0.2);

	// sys_c4_3sub_Har4_Pt1
	smooth_ratio(ratio_c4_3sub[6][4][1],2,0.2);
	smooth_ratio(ratio_c4_3sub[3][4][1],6,0.2);
	smooth_ratio(ratio_c4_3sub[3][4][1],7,0.2);
	smooth_ratio(ratio_c4_3sub[6][4][1],6,0.2);
	smooth_ratio(ratio_c4_3sub[6][4][1],7,0.2);

	// sys_c6_1sub_Har2_Pt1
	smooth_ratio(ratio_c6_1sub[3][2][1],7,0.05);

	// sys_c6_1sub_Har3_Pt0
	reduce_ratio(ratio_c6_1sub[3][3][0],5,0.5);

	// sys_c6_1sub_Har3_Pt1
	smooth_ratio(ratio_c6_1sub[3][3][1],5,0.1);
	smooth_ratio(ratio_c6_1sub[4][3][1],6,0.1);

	// sys_c6_1sub_Har4_Pt1
	reduce_ratio(ratio_c6_1sub[3][4][1],6,0.5);
	reduce_ratio(ratio_c6_1sub[4][4][1],6,0.5);

	// sys_nc4_1sub_Har3_Pt0
	smooth_ratio(ratio_nc4_1sub[3][3][0],6,0.1);
	smooth_ratio(ratio_nc4_1sub[6][3][0],7,0.1);

	// sys_nc4_1sub_Har3_Pt1
	smooth_ratio(ratio_nc4_1sub[3][3][1],6,0.1);
	smooth_ratio(ratio_nc4_1sub[6][3][1],7,0.1);

	// sys_nc4_1sub_Har4_Pt0
	smooth_ratio(ratio_nc4_1sub[6][4][0],0,0.5);
	smooth_ratio(ratio_nc4_1sub[6][4][0],1,0.5);
	smooth_ratio(ratio_nc4_1sub[2][4][0],7,0.2);
	smooth_ratio(ratio_nc4_1sub[3][4][0],7,0.2);

	// sys_nc4_1sub_Har4_Pt1
	smooth_ratio(ratio_nc4_1sub[2][4][1],7,0.05);
	smooth_ratio(ratio_nc4_1sub[3][4][1],7,0.05);
	
	// sys_nc4_3sub_Har3_Pt0
	smooth_ratio(ratio_nc4_3sub[3][3][0],6,0.1);
	smooth_ratio(ratio_nc4_3sub[4][3][0],6,0.1);
	smooth_ratio(ratio_nc4_3sub[6][3][0],6,0.1);
	smooth_ratio(ratio_nc4_3sub[6][3][0],7,0.1);

	// sys_nc4_3sub_Har3_Pt1
	smooth_ratio(ratio_nc4_3sub[4][3][1],6,0.1);
	smooth_ratio(ratio_nc4_3sub[6][3][1],6,0.1);
	reduce_ratio(ratio_nc4_3sub[4][3][1],6,0.2);
	
	// sys_nc4_3sub_Har4_Pt0
	smooth_ratio(ratio_nc4_3sub[3][4][0],2,0.2);
	smooth_ratio(ratio_nc4_3sub[6][4][0],2,0.2);
	reduce_ratio(ratio_nc4_3sub[6][4][0],2,0.5);

	// sys_nc4_3sub_Har4_Pt1
	smooth_ratio(ratio_nc4_3sub[6][4][1],2,0.2);
	smooth_ratio(ratio_nc4_3sub[3][4][1],6,0.1);
	smooth_ratio(ratio_nc4_3sub[3][4][1],7,0.1);
	smooth_ratio(ratio_nc4_3sub[6][4][1],6,0.1);
	smooth_ratio(ratio_nc4_3sub[6][4][1],7,0.1);
	
	// sys_nc6_1sub_Har3_Pt1
	smooth_ratio(ratio_nc6_1sub[3][3][1],5,0.2);
	reduce_ratio(ratio_nc6_1sub[4][3][1],6,0.1);

	// sys_nc6_1sub_Har4_Pt1
	reduce_ratio(ratio_nc6_1sub[3][4][1],6,0.5);
	reduce_ratio(ratio_nc6_1sub[4][4][1],6,0.5);

	// sys_nsc_1sub_Har2_Pt1
	smooth_ratio(ratio_nsc_1sub[3][2][1],6,0.2);
	smooth_ratio(ratio_nsc_1sub[3][2][1],7,0.2);
	smooth_ratio(ratio_nsc_1sub[6][2][1],7,0.2);

	// sys_nsc_3sub_Har2_Pt1
	smooth_ratio(ratio_nsc_3sub[3][2][1],7,0.05);

	// sys_nsc_3sub_Har3_Pt1
	smooth_ratio(ratio_nsc_3sub[3][3][1],7,0.05);

	// sys_sc_3sub_Har2_Pt1
	smooth_ratio(ratio_sc_3sub[3][2][1],7,0.05);
	smooth_ratio(ratio_sc_3sub[6][2][1],7,0.1);

	// sys_sc_3sub_Har3_Pt1
	smooth_ratio(ratio_sc_3sub[3][3][1],7,0.05);
	
	// special treatment for flattening
	for(unsigned int iS=6; iS<7; iS++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				for(int iB=0; iB<ratio_c2_1sub[0][0][0]->GetN(); iB++)
				{
					reduce_ratio(ratio_c2_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_c4_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_c6_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nc4_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nc6_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_sc_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nsc_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_ac_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nac_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_isGauss_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_isPower_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_c2_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_c4_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nc4_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_sc_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nsc_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_ac_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nac_3sub[iS][iV][iP],iB,0.5);
				}
			}
		}
	}
	
	// speical treatment for MC closure
	for(unsigned int iS=5; iS<6; iS++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				for(int iB=0; iB<ratio_c2_1sub[0][0][0]->GetN(); iB++)
				{
					smooth_ratio(ratio_c2_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_c4_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_c6_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nc4_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nc6_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_sc_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nsc_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_ac_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nac_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_isGauss_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_isPower_1sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_c2_3sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_c4_3sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nc4_3sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_sc_3sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nsc_3sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_ac_3sub[iS][iV][iP],iB,0.05);
					smooth_ratio(ratio_nac_3sub[iS][iV][iP],iB,0.05);

					reduce_ratio(ratio_sc_1sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_nsc_1sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_ac_1sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_nac_1sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_sc_3sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_nsc_3sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_ac_3sub[iS][iV][iP],iB,1E-9);
					reduce_ratio(ratio_nac_3sub[iS][iV][iP],iB,1E-9);
					
					reduce_ratio(ratio_c4_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_c6_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nc4_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nc6_1sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_c4_3sub[iS][iV][iP],iB,0.5);
					reduce_ratio(ratio_nc4_3sub[iS][iV][iP],iB,0.5);
				}
			}
		}
	}
	reduce_ratio(ratio_c2_3sub[5][2][1],1,0.3);
	reduce_ratio(ratio_c2_3sub[5][4][1],3,0.5);
	reduce_ratio(ratio_c2_3sub[5][4][1],4,0.5);
	reduce_ratio(ratio_c4_1sub[5][3][0],7,0.001);
	reduce_ratio(ratio_c4_1sub[5][3][1],7,0.1);
	reduce_ratio(ratio_c4_3sub[5][3][1],7,0.5);
	reduce_ratio(ratio_c4_3sub[5][4][1],6,0.05);
	reduce_ratio(ratio_c6_1sub[5][3][1],6,0.5);
	reduce_ratio(ratio_nc4_1sub[5][3][0],7,0.001);
	reduce_ratio(ratio_nc4_1sub[5][3][1],7,0.01);
	reduce_ratio(ratio_nc4_3sub[5][3][1],7,0.0005);
	reduce_ratio(ratio_nc4_3sub[5][4][1],6,0.05);
	reduce_ratio(ratio_nc6_1sub[5][3][1],6,0.5);
	/*
	reduce_ratio(ratio_nsc_1sub[5][2][0],3,0.5);
	reduce_ratio(ratio_nsc_1sub[5][2][0],5,0.003);
	reduce_ratio(ratio_nsc_3sub[5][2][0],2,0.5);
	reduce_ratio(ratio_nsc_3sub[5][3][0],3,0.2);
	reduce_ratio(ratio_nsc_3sub[5][3][0],4,0.5);
	reduce_ratio(ratio_nsc_3sub[5][3][0],5,0.01);
	reduce_ratio(ratio_nsc_3sub[5][3][1],1,0.02);
	reduce_ratio(ratio_nsc_3sub[5][3][1],2,0.02);
	reduce_ratio(ratio_nsc_3sub[5][3][1],5,0.005);
	reduce_ratio(ratio_sc_1sub[5][2][0],3,0.1);
	reduce_ratio(ratio_sc_1sub[5][2][0],5,0.001);
	reduce_ratio(ratio_sc_1sub[5][2][1],0,0.2);
	reduce_ratio(ratio_sc_1sub[5][2][1],2,0.5);
	reduce_ratio(ratio_sc_1sub[5][2][1],4,0.5);
	reduce_ratio(ratio_sc_1sub[5][2][1],7,0.1);
	reduce_ratio(ratio_sc_1sub[5][3][0],0,0.5);
	reduce_ratio(ratio_sc_1sub[5][3][0],1,0.2);
	reduce_ratio(ratio_sc_1sub[5][3][0],2,0.2);
	reduce_ratio(ratio_sc_1sub[5][3][1],1,0.2);
	reduce_ratio(ratio_sc_1sub[5][3][1],3,0.5);
	reduce_ratio(ratio_sc_1sub[5][3][1],4,0.5);
	reduce_ratio(ratio_sc_3sub[5][2][0],3,0.2);
	reduce_ratio(ratio_sc_3sub[5][2][0],4,0.04);
	reduce_ratio(ratio_sc_3sub[5][2][0],5,0.5);
	reduce_ratio(ratio_sc_3sub[5][2][1],0,0.5);
	reduce_ratio(ratio_sc_3sub[5][2][1],1,0.2);
	reduce_ratio(ratio_sc_3sub[5][2][1],3,0.04);
	reduce_ratio(ratio_sc_3sub[5][2][1],4,0.04);
	reduce_ratio(ratio_sc_3sub[5][2][1],7,0.5);
	reduce_ratio(ratio_sc_3sub[5][3][0],0,0.3);
	reduce_ratio(ratio_sc_3sub[5][3][0],1,0.1);
	reduce_ratio(ratio_sc_3sub[5][3][0],3,0.2);
	reduce_ratio(ratio_sc_3sub[5][3][0],4,0.01);
	reduce_ratio(ratio_sc_3sub[5][3][0],5,0.08);
	reduce_ratio(ratio_sc_3sub[5][3][1],1,0.3);
	reduce_ratio(ratio_sc_3sub[5][3][1],3,0.5);
	reduce_ratio(ratio_sc_3sub[5][3][1],4,0.1);
	reduce_ratio(ratio_sc_3sub[5][3][1],5,0.002);
	reduce_ratio(ratio_sc_3sub[5][3][1],6,0.05);
	*/
}
