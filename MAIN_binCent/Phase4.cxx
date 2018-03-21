#include "Phase4.h"

Phase4::Phase4(unsigned int iBin)
{
	initialize(iBin);
	execute();
	finalize(iBin);
}

Phase4::~Phase4()
{
}

void Phase4::execute()
{
	cout<<"execute..."<<endl;
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			combine(sts_c2_1sub[iV][iP],sysUp_c2_1sub[iV][iP],sysLw_c2_1sub[iV][iP],sys_c2_1sub[iV][iP],"sys_c2_1sub",iV,iP);
			combine(sts_c4_1sub[iV][iP],sysUp_c4_1sub[iV][iP],sysLw_c4_1sub[iV][iP],sys_c4_1sub[iV][iP],"sys_c4_1sub",iV,iP);
			combine(sts_c6_1sub[iV][iP],sysUp_c6_1sub[iV][iP],sysLw_c6_1sub[iV][iP],sys_c6_1sub[iV][iP],"sys_c6_1sub",iV,iP);
			combine(sts_nc4_1sub[iV][iP],sysUp_nc4_1sub[iV][iP],sysLw_nc4_1sub[iV][iP],sys_nc4_1sub[iV][iP],"sys_nc4_1sub",iV,iP);
			combine(sts_nc6_1sub[iV][iP],sysUp_nc6_1sub[iV][iP],sysLw_nc6_1sub[iV][iP],sys_nc6_1sub[iV][iP],"sys_nc6_1sub",iV,iP);
			combine(sts_sc_1sub[iV][iP],sysUp_sc_1sub[iV][iP],sysLw_sc_1sub[iV][iP],sys_sc_1sub[iV][iP],"sys_sc_1sub",iV,iP);
			combine(sts_nsc_1sub[iV][iP],sysUp_nsc_1sub[iV][iP],sysLw_nsc_1sub[iV][iP],sys_nsc_1sub[iV][iP],"sys_nsc_1sub",iV,iP);
			combine(sts_ac_1sub[iV][iP],sysUp_ac_1sub[iV][iP],sysLw_ac_1sub[iV][iP],sys_ac_1sub[iV][iP],"sys_ac_1sub",iV,iP);
			combine(sts_nac_1sub[iV][iP],sysUp_nac_1sub[iV][iP],sysLw_nac_1sub[iV][iP],sys_nac_1sub[iV][iP],"sys_nac_1sub",iV,iP);
			combine(sts_isGauss_1sub[iV][iP],sysUp_isGauss_1sub[iV][iP],sysLw_isGauss_1sub[iV][iP],sys_isGauss_1sub[iV][iP],"sys_isGauss_1sub",iV,iP);
			combine(sts_isPower_1sub[iV][iP],sysUp_isPower_1sub[iV][iP],sysLw_isPower_1sub[iV][iP],sys_isPower_1sub[iV][iP],"sys_isPower_1sub",iV,iP);
			combine(sts_c2_3sub[iV][iP],sysUp_c2_3sub[iV][iP],sysLw_c2_3sub[iV][iP],sys_c2_3sub[iV][iP],"sys_c2_3sub",iV,iP);
			combine(sts_c4_3sub[iV][iP],sysUp_c4_3sub[iV][iP],sysLw_c4_3sub[iV][iP],sys_c4_3sub[iV][iP],"sys_c4_3sub",iV,iP);
			combine(sts_nc4_3sub[iV][iP],sysUp_nc4_3sub[iV][iP],sysLw_nc4_3sub[iV][iP],sys_nc4_3sub[iV][iP],"sys_nc4_3sub",iV,iP);
			combine(sts_sc_3sub[iV][iP],sysUp_sc_3sub[iV][iP],sysLw_sc_3sub[iV][iP],sys_sc_3sub[iV][iP],"sys_sc_3sub",iV,iP);
			combine(sts_nsc_3sub[iV][iP],sysUp_nsc_3sub[iV][iP],sysLw_nsc_3sub[iV][iP],sys_nsc_3sub[iV][iP],"sys_nsc_3sub",iV,iP);
			combine(sts_ac_3sub[iV][iP],sysUp_ac_3sub[iV][iP],sysLw_ac_3sub[iV][iP],sys_ac_3sub[iV][iP],"sys_ac_3sub",iV,iP);
			combine(sts_nac_3sub[iV][iP],sysUp_nac_3sub[iV][iP],sysLw_nac_3sub[iV][iP],sys_nac_3sub[iV][iP],"sys_nac_3sub",iV,iP);

			trim(sts_c2_1sub[iV][iP],sys_c2_1sub[iV][iP]);
			trim(sts_c4_1sub[iV][iP],sys_c4_1sub[iV][iP]);
			trim(sts_c6_1sub[iV][iP],sys_c6_1sub[iV][iP]);
			trim(sts_nc4_1sub[iV][iP],sys_nc4_1sub[iV][iP]);
			trim(sts_nc6_1sub[iV][iP],sys_nc6_1sub[iV][iP]);
			trim(sts_sc_1sub[iV][iP],sys_sc_1sub[iV][iP]);
			trim(sts_nsc_1sub[iV][iP],sys_nsc_1sub[iV][iP]);
			trim(sts_ac_1sub[iV][iP],sys_ac_1sub[iV][iP]);
			trim(sts_nac_1sub[iV][iP],sys_nac_1sub[iV][iP]);
			trim(sts_isGauss_1sub[iV][iP],sys_isGauss_1sub[iV][iP]);
			trim(sts_isPower_1sub[iV][iP],sys_isPower_1sub[iV][iP]);
			trim(sts_c2_3sub[iV][iP],sys_c2_3sub[iV][iP]);
			trim(sts_c4_3sub[iV][iP],sys_c4_3sub[iV][iP]);
			trim(sts_nc4_3sub[iV][iP],sys_nc4_3sub[iV][iP]);
			trim(sts_sc_3sub[iV][iP],sys_sc_3sub[iV][iP]);
			trim(sts_nsc_3sub[iV][iP],sys_nsc_3sub[iV][iP]);
			trim(sts_ac_3sub[iV][iP],sys_ac_3sub[iV][iP]);
			trim(sts_nac_3sub[iV][iP],sys_nac_3sub[iV][iP]);
		}
	}
}

void Phase4::initialize(unsigned int iBin)
{
	cout<<"initialize..."<<endl;

	sprintf(name,"../OUTPUT/Phase3/Phase3_bin%d.root",iBin);
	TFile* fSts = new TFile(name,"READ");
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			readHist_VP(fSts,sts_c2_1sub[iV][iP],"c2_1sub","sts_c2_1sub",iV,iP);
			readHist_VP(fSts,sts_c4_1sub[iV][iP],"c4_1sub","sts_c4_1sub",iV,iP);
			readHist_VP(fSts,sts_c6_1sub[iV][iP],"c6_1sub","sts_c6_1sub",iV,iP);
			readHist_VP(fSts,sts_nc4_1sub[iV][iP],"nc4_1sub","sts_nc4_1sub",iV,iP);
			readHist_VP(fSts,sts_nc6_1sub[iV][iP],"nc6_1sub","sts_nc6_1sub",iV,iP);
			readHist_VP(fSts,sts_sc_1sub[iV][iP],"sc_1sub","sts_sc_1sub",iV,iP);
			readHist_VP(fSts,sts_nsc_1sub[iV][iP],"nsc_1sub","sts_nsc_1sub",iV,iP);
			readHist_VP(fSts,sts_ac_1sub[iV][iP],"ac_1sub","sts_ac_1sub",iV,iP);
			readHist_VP(fSts,sts_nac_1sub[iV][iP],"nac_1sub","sts_nac_1sub",iV,iP);
			readHist_VP(fSts,sts_isGauss_1sub[iV][iP],"isGauss_1sub","sts_isGauss_1sub",iV,iP);
			readHist_VP(fSts,sts_isPower_1sub[iV][iP],"isPower_1sub","sts_isPower_1sub",iV,iP);
			readHist_VP(fSts,sts_c2_3sub[iV][iP],"c2_3sub","sts_c2_3sub",iV,iP);
			readHist_VP(fSts,sts_c4_3sub[iV][iP],"c4_3sub","sts_c4_3sub",iV,iP);
			readHist_VP(fSts,sts_nc4_3sub[iV][iP],"nc4_3sub","sts_nc4_3sub",iV,iP);
			readHist_VP(fSts,sts_sc_3sub[iV][iP],"sc_3sub","sts_sc_3sub",iV,iP);
			readHist_VP(fSts,sts_nsc_3sub[iV][iP],"nsc_3sub","sts_nsc_3sub",iV,iP);
			readHist_VP(fSts,sts_ac_3sub[iV][iP],"ac_3sub","sts_ac_3sub",iV,iP);
			readHist_VP(fSts,sts_nac_3sub[iV][iP],"nac_3sub","sts_nac_3sub",iV,iP);
		}
	}

	sprintf(name,"/phenix/plhf/mzhou/AnaCumu_PbPb502/SYS/CombSys/OUTPUT/Sys_bin0.root");
	TFile* fSys = new TFile(name,"READ");
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			int iPsel = 0;
			if(iP<2) iPsel = 0;
			else iPsel = 1;
			readHist_VP(fSys,sysUp_c2_1sub[iV][iP],"sysUp_c2_1sub","sysUp_c2_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_c4_1sub[iV][iP],"sysUp_c4_1sub","sysUp_c4_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_c6_1sub[iV][iP],"sysUp_c6_1sub","sysUp_c6_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nc4_1sub[iV][iP],"sysUp_nc4_1sub","sysUp_nc4_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nc6_1sub[iV][iP],"sysUp_nc6_1sub","sysUp_nc6_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_sc_1sub[iV][iP],"sysUp_sc_1sub","sysUp_sc_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nsc_1sub[iV][iP],"sysUp_nsc_1sub","sysUp_nsc_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_ac_1sub[iV][iP],"sysUp_ac_1sub","sysUp_ac_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nac_1sub[iV][iP],"sysUp_nac_1sub","sysUp_nac_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_isGauss_1sub[iV][iP],"sysUp_isGauss_1sub","sysUp_isGauss_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_isPower_1sub[iV][iP],"sysUp_isPower_1sub","sysUp_isPower_1sub",iV,iPsel);
			readHist_VP(fSys,sysUp_c2_3sub[iV][iP],"sysUp_c2_3sub","sysUp_c2_3sub",iV,iPsel);
			readHist_VP(fSys,sysUp_c4_3sub[iV][iP],"sysUp_c4_3sub","sysUp_c4_3sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nc4_3sub[iV][iP],"sysUp_nc4_3sub","sysUp_nc4_3sub",iV,iPsel);
			readHist_VP(fSys,sysUp_sc_3sub[iV][iP],"sysUp_sc_3sub","sysUp_sc_3sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nsc_3sub[iV][iP],"sysUp_nsc_3sub","sysUp_nsc_3sub",iV,iPsel);
			readHist_VP(fSys,sysUp_ac_3sub[iV][iP],"sysUp_ac_3sub","sysUp_ac_3sub",iV,iPsel);
			readHist_VP(fSys,sysUp_nac_3sub[iV][iP],"sysUp_nac_3sub","sysUp_nac_3sub",iV,iPsel);
			
			readHist_VP(fSys,sysLw_c2_1sub[iV][iP],"sysLw_c2_1sub","sysLw_c2_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_c4_1sub[iV][iP],"sysLw_c4_1sub","sysLw_c4_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_c6_1sub[iV][iP],"sysLw_c6_1sub","sysLw_c6_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nc4_1sub[iV][iP],"sysLw_nc4_1sub","sysLw_nc4_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nc6_1sub[iV][iP],"sysLw_nc6_1sub","sysLw_nc6_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_sc_1sub[iV][iP],"sysLw_sc_1sub","sysLw_sc_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nsc_1sub[iV][iP],"sysLw_nsc_1sub","sysLw_nsc_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_ac_1sub[iV][iP],"sysLw_ac_1sub","sysLw_ac_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nac_1sub[iV][iP],"sysLw_nac_1sub","sysLw_nac_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_isGauss_1sub[iV][iP],"sysLw_isGauss_1sub","sysLw_isGauss_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_isPower_1sub[iV][iP],"sysLw_isPower_1sub","sysLw_isPower_1sub",iV,iPsel);
			readHist_VP(fSys,sysLw_c2_3sub[iV][iP],"sysLw_c2_3sub","sysLw_c2_3sub",iV,iPsel);
			readHist_VP(fSys,sysLw_c4_3sub[iV][iP],"sysLw_c4_3sub","sysLw_c4_3sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nc4_3sub[iV][iP],"sysLw_nc4_3sub","sysLw_nc4_3sub",iV,iPsel);
			readHist_VP(fSys,sysLw_sc_3sub[iV][iP],"sysLw_sc_3sub","sysLw_sc_3sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nsc_3sub[iV][iP],"sysLw_nsc_3sub","sysLw_nsc_3sub",iV,iPsel);
			readHist_VP(fSys,sysLw_ac_3sub[iV][iP],"sysLw_ac_3sub","sysLw_ac_3sub",iV,iPsel);
			readHist_VP(fSys,sysLw_nac_3sub[iV][iP],"sysLw_nac_3sub","sysLw_nac_3sub",iV,iPsel);
		}
	}

	TFile* fCvt = new TFile("../INPUT/hist_cvt.root","READ");
	gCvt = (TGraphErrors*)fCvt->Get("gCvt_Cent_FCal");
}

void Phase4::finalize(unsigned int iBin)
{
	cout<<"finalize..."<<endl;

	sprintf(name,"../OUTPUT/Phase4/hist_PbPb502_binCent_bin%d.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			sts_c2_1sub[iV][iP]->Write();
			sts_c4_1sub[iV][iP]->Write();
			sts_c6_1sub[iV][iP]->Write();
			sts_nc4_1sub[iV][iP]->Write();
			sts_nc6_1sub[iV][iP]->Write();
			sts_sc_1sub[iV][iP]->Write();
			sts_nsc_1sub[iV][iP]->Write();
			sts_ac_1sub[iV][iP]->Write();
			sts_nac_1sub[iV][iP]->Write();
			sts_isGauss_1sub[iV][iP]->Write();
			sts_isPower_1sub[iV][iP]->Write();
			sts_c2_3sub[iV][iP]->Write();
			sts_c4_3sub[iV][iP]->Write();
			sts_nc4_3sub[iV][iP]->Write();
			sts_sc_3sub[iV][iP]->Write();
			sts_nsc_3sub[iV][iP]->Write();
			sts_ac_3sub[iV][iP]->Write();
			sts_nac_3sub[iV][iP]->Write();

			sys_c2_1sub[iV][iP]->Write();
			sys_c4_1sub[iV][iP]->Write();
			sys_c6_1sub[iV][iP]->Write();
			sys_nc4_1sub[iV][iP]->Write();
			sys_nc6_1sub[iV][iP]->Write();
			sys_sc_1sub[iV][iP]->Write();
			sys_nsc_1sub[iV][iP]->Write();
			sys_ac_1sub[iV][iP]->Write();
			sys_nac_1sub[iV][iP]->Write();
			sys_isGauss_1sub[iV][iP]->Write();
			sys_isPower_1sub[iV][iP]->Write();
			sys_c2_3sub[iV][iP]->Write();
			sys_c4_3sub[iV][iP]->Write();
			sys_nc4_3sub[iV][iP]->Write();
			sys_sc_3sub[iV][iP]->Write();
			sys_nsc_3sub[iV][iP]->Write();
			sys_ac_3sub[iV][iP]->Write();
			sys_nac_3sub[iV][iP]->Write();
		}
	}
	fOut->Close();
}

void Phase4::combine(TGraphErrors* gSts, TGraphErrors* gSysUp, TGraphErrors* gSysLw, TGraphAsymmErrors*& gSys, const char* hName, int iV, int iP)
{
	int NB = gSts->GetN();
	double x[NB];
	double xElw[NB];
	double xEup[NB];
	double y[NB];
	double yElw[NB];
	double yEup[NB];
	double xLw;
	double xUp;
	gSysLw->GetPoint(0,xLw,x[0]);
	gSysLw->GetPoint(gSysLw->GetN()-1,xUp,x[0]);
	for(int iB=0; iB<NB; iB++)
	{
		gSts->GetPoint(iB,x[iB],y[iB]);
		xElw[iB] = 0;
		xEup[iB] = 0;
		double xCvt = x[iB];
		//double xCvt = gCvt->Eval(x[iB]);
		if(xCvt<xLw) xCvt = xLw;
		if(xCvt>xUp) xCvt = xUp;
		yElw[iB] = fabs(y[iB]*gSysLw->Eval(xCvt));
		yEup[iB] = fabs(y[iB]*gSysUp->Eval(xCvt));
	}
	gSys = new TGraphAsymmErrors(NB,x,y,xElw,xEup,yElw,yEup);
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	gSys->SetName(name);
}

void Phase4::readHist_VP(TFile* fIn, TGraphErrors*& hIn, const char* hName, const char* hNewName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
	sprintf(name,"%s_Har%d_Pt%d",hNewName,iV,iP);
	hIn->SetName(name);
}

void Phase4::trim(TGraphErrors* gSts, TGraphAsymmErrors* gSys)
{
	for(int i=0; i<gSts->GetN(); i++)
	{
		if(gSys->GetErrorYhigh(i)<0.2*gSts->GetErrorY(i)) gSys->SetPointEYhigh(i,0.2*gSts->GetErrorY(i));
		if(gSys->GetErrorYlow(i)< 0.2*gSts->GetErrorY(i)) gSys->SetPointEYlow(i,0.2*gSts->GetErrorY(i));
	}
}


