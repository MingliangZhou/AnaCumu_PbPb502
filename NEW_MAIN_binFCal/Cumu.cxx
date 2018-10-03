#include "Cumu.h"

Cumu::Cumu(TFile* fIn, unsigned int iRebin)
{
/*
	for(unsigned int iP=0; iP<NP; iP++)
	{
		sprintf(name,"cnt_1sub_Pt%d",iP);
		cnt_1sub[iP] = (TH1D*)fIn->Get(name);
		sprintf(name,"raw_cnt_1sub_Pt%d",iP);
		raw_cnt_1sub[iP] = (TH1D*)cnt_1sub[iP]->Clone(name);
		raw_cnt_1sub[iP]->Rebin(nBin/nRebin[0]);
		sprintf(name,"rbn_cnt_1sub_Pt%d",iP);
		rbn_cnt_1sub[iP] = (TH1D*)cnt_1sub[iP]->Clone(name);
		rbn_cnt_1sub[iP]->Rebin(nBin/nRebin[iRebin]);
		for(unsigned int iV=0; iV<NV; iV++)
		{
			readHist_VP(fIn,pc2_1sub_mean[iV][iP],"pc2_1sub_mean",iV,iP);
			readHist_VP(fIn,pc2_1sub_wght[iV][iP],"pc2_1sub_wght",iV,iP);
			readHist_VP(fIn,pc4_1sub_mean[iV][iP],"pc4_1sub_mean",iV,iP);
			readHist_VP(fIn,pc4_1sub_wght[iV][iP],"pc4_1sub_wght",iV,iP);
			readHist_VP(fIn,pc6_1sub_mean[iV][iP],"pc6_1sub_mean",iV,iP);
			readHist_VP(fIn,pc6_1sub_wght[iV][iP],"pc6_1sub_wght",iV,iP);
			readHist_VP(fIn,psc4_1sub_mean[iV][iP],"psc4_1sub_mean",iV,iP);
			readHist_VP(fIn,psc4_1sub_wght[iV][iP],"psc4_1sub_wght",iV,iP);
			readHist_VP(fIn,pac3_1sub_mean[iV][iP],"pac3_1sub_mean",iV,iP);
			readHist_VP(fIn,pac3_1sub_wght[iV][iP],"pac3_1sub_wght",iV,iP);
		}
	}
*/
	for(unsigned int iP=0; iP<NP; iP++)
	{
		for(unsigned int iA=0; iA<NA; iA++)
		{
			sprintf(name,"cnt_3sub_Pm%d_Pt%d",iA,iP);
			cnt_3sub[iA][iP] = (TH1D*)fIn->Get(name);
			sprintf(name,"raw_cnt_3sub_Pm%d_Pt%d",iA,iP);
			raw_cnt_3sub[iA][iP] = (TH1D*)cnt_3sub[iA][iP]->Clone(name);
			raw_cnt_3sub[iA][iP]->Rebin(nBin/nRebin[0]);
			sprintf(name,"rbn_cnt_3sub_Pm%d_Pt%d",iA,iP);
			rbn_cnt_3sub[iA][iP] = (TH1D*)cnt_3sub[iA][iP]->Clone(name);
			rbn_cnt_3sub[iA][iP]->Rebin(nBin/nRebin[iRebin]);
			for(unsigned int iV=0; iV<NV; iV++)
			{
				readHist_AVP(fIn,pc2_1_3sub_mean[iA][iV][iP],"pc2_1_3sub_mean",iA,iV,iP);
				readHist_AVP(fIn,pc2_1_3sub_wght[iA][iV][iP],"pc2_1_3sub_wght",iA,iV,iP);
				readHist_AVP(fIn,pc2_2_3sub_mean[iA][iV][iP],"pc2_2_3sub_mean",iA,iV,iP);
				readHist_AVP(fIn,pc2_2_3sub_wght[iA][iV][iP],"pc2_2_3sub_wght",iA,iV,iP);
				readHist_AVP(fIn,pc4_3sub_mean[iA][iV][iP],"pc4_3sub_mean",iA,iV,iP);
				readHist_AVP(fIn,pc4_3sub_wght[iA][iV][iP],"pc4_3sub_wght",iA,iV,iP);
				readHist_AVP(fIn,psc4_3sub_Type1_mean[iA][iV][iP],"psc4_3sub_Type1_mean",iA,iV,iP);
				readHist_AVP(fIn,psc4_3sub_Type2_mean[iA][iV][iP],"psc4_3sub_Type2_mean",iA,iV,iP);
				readHist_AVP(fIn,psc4_3sub_wght[iA][iV][iP],"psc4_3sub_wght",iA,iV,iP);
				readHist_AVP(fIn,pac3_3sub_mean[iA][iV][iP],"pac4_3sub_mean",iA,iV,iP);
				readHist_AVP(fIn,pac3_3sub_wght[iA][iV][iP],"pac4_3sub_wght",iA,iV,iP);
			}
		}
	}
	
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			initHist_VP(ori_c2_1sub[iV][iP],"ori_c2_1sub",iV,iP,nBin);
			initHist_VP(raw_c2_1sub[iV][iP],"raw_c2_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_c2_1sub[iV][iP],"rbn_c2_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(ori_c4_1sub[iV][iP],"ori_c4_1sub",iV,iP,nBin);
			initHist_VP(raw_c4_1sub[iV][iP],"raw_c4_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_c4_1sub[iV][iP],"rbn_c4_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(ori_c6_1sub[iV][iP],"ori_c6_1sub",iV,iP,nBin);
			initHist_VP(raw_c6_1sub[iV][iP],"raw_c6_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_c6_1sub[iV][iP],"rbn_c6_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_nc4_1sub[iV][iP],"raw_nc4_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_nc4_1sub[iV][iP],"rbn_nc4_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_nc6_1sub[iV][iP],"raw_nc6_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_nc6_1sub[iV][iP],"rbn_nc6_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(ori_sc_1sub[iV][iP],"ori_sc_1sub",iV,iP,nBin);
			initHist_VP(raw_sc_1sub[iV][iP],"raw_sc_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_sc_1sub[iV][iP],"rbn_sc_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_nsc_1sub[iV][iP],"raw_nsc_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_nsc_1sub[iV][iP],"rbn_nsc_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_ac_1sub[iV][iP],"raw_ac_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_ac_1sub[iV][iP],"rbn_ac_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_nac_1sub[iV][iP],"raw_nac_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_nac_1sub[iV][iP],"rbn_nac_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_isGauss_1sub[iV][iP],"raw_isGauss_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_isGauss_1sub[iV][iP],"rbn_isGauss_1sub",iV,iP,nRebin[iRebin]);
			initHist_VP(raw_isPower_1sub[iV][iP],"raw_isPower_1sub",iV,iP,nRebin[0]);
			initHist_VP(rbn_isPower_1sub[iV][iP],"rbn_isPower_1sub",iV,iP,nRebin[iRebin]);
		}
	}

	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				initHist_AVP(ori_c2_3sub[iA][iV][iP],"ori_c2_3sub",iA,iV,iP,nBin);
				initHist_AVP(raw_c2_3sub[iA][iV][iP],"raw_c2_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_c2_3sub[iA][iV][iP],"rbn_c2_3sub",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(ori_c4_3sub[iA][iV][iP],"ori_c4_3sub",iA,iV,iP,nBin);
				initHist_AVP(raw_c4_3sub[iA][iV][iP],"raw_c4_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(raw_p4_3sub[iA][iV][iP],"raw_p4_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_c4_3sub[iA][iV][iP],"rbn_c4_3sub",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(raw_nc4_3sub[iA][iV][iP],"raw_nc4_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_nc4_3sub[iA][iV][iP],"rbn_nc4_3sub",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(raw_nc6_3sub[iA][iV][iP],"raw_nc6_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_nc6_3sub[iA][iV][iP],"rbn_nc6_3sub",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(ori_sc_3sub_Type1[iA][iV][iP],"ori_sc_3sub_Type1",iA,iV,iP,nBin);
				initHist_AVP(ori_sc_3sub_Type2[iA][iV][iP],"ori_sc_3sub_Type2",iA,iV,iP,nBin);
				initHist_AVP(raw_sc_3sub_Type1[iA][iV][iP],"raw_sc_3sub_Type1",iA,iV,iP,nRebin[0]);
				initHist_AVP(raw_sc_3sub_Type2[iA][iV][iP],"raw_sc_3sub_Type2",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_sc_3sub_Type1[iA][iV][iP],"rbn_sc_3sub_Type1",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(rbn_sc_3sub_Type2[iA][iV][iP],"rbn_sc_3sub_Type2",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(raw_nsc_3sub_Type1[iA][iV][iP],"raw_nsc_3sub_Type1",iA,iV,iP,nRebin[0]);
				initHist_AVP(raw_nsc_3sub_Type2[iA][iV][iP],"raw_nsc_3sub_Type2",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_nsc_3sub_Type1[iA][iV][iP],"rbn_nsc_3sub_Type1",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(rbn_nsc_3sub_Type2[iA][iV][iP],"rbn_nsc_3sub_Type2",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(raw_ac_3sub[iA][iV][iP],"raw_ac_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_ac_3sub[iA][iV][iP],"rbn_ac_3sub",iA,iV,iP,nRebin[iRebin]);
				initHist_AVP(raw_nac_3sub[iA][iV][iP],"raw_nac_3sub",iA,iV,iP,nRebin[0]);
				initHist_AVP(rbn_nac_3sub[iA][iV][iP],"rbn_nac_3sub",iA,iV,iP,nRebin[iRebin]);
			}
		}
	}
}

Cumu::~Cumu()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			delete ori_c2_1sub[iV][iP];
			delete raw_c2_1sub[iV][iP];
			delete rbn_c2_1sub[iV][iP];
			delete ori_c4_1sub[iV][iP];
			delete raw_c4_1sub[iV][iP];
			delete rbn_c4_1sub[iV][iP];
			delete ori_c6_1sub[iV][iP];
			delete raw_c6_1sub[iV][iP];
			delete rbn_c6_1sub[iV][iP];
			delete raw_nc4_1sub[iV][iP];
			delete rbn_nc4_1sub[iV][iP];
			delete raw_nc6_1sub[iV][iP];
			delete rbn_nc6_1sub[iV][iP];
			delete ori_sc_1sub[iV][iP];
			delete raw_sc_1sub[iV][iP];
			delete rbn_sc_1sub[iV][iP];
			delete raw_nsc_1sub[iV][iP];
			delete rbn_nsc_1sub[iV][iP];
			delete raw_ac_1sub[iV][iP];
			delete rbn_ac_1sub[iV][iP];
			delete raw_nac_1sub[iV][iP];
			delete rbn_nac_1sub[iV][iP];
			delete raw_isGauss_1sub[iV][iP];
			delete rbn_isGauss_1sub[iV][iP];
			delete raw_isPower_1sub[iV][iP];
			delete rbn_isPower_1sub[iV][iP];
		}
	}

	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				delete ori_c2_3sub[iA][iV][iP];
				delete raw_c2_3sub[iA][iV][iP];
				delete rbn_c2_3sub[iA][iV][iP];
				delete ori_c4_3sub[iA][iV][iP];
				delete raw_c4_3sub[iA][iV][iP];
				delete rbn_c4_3sub[iA][iV][iP];
				delete raw_nc4_3sub[iA][iV][iP];
				delete rbn_nc4_3sub[iA][iV][iP];
				delete raw_nc6_3sub[iA][iV][iP];
				delete rbn_nc6_3sub[iA][iV][iP];
				delete ori_sc_3sub_Type1[iA][iV][iP];
				delete ori_sc_3sub_Type2[iA][iV][iP];
				delete raw_sc_3sub_Type1[iA][iV][iP];
				delete raw_sc_3sub_Type2[iA][iV][iP];
				delete rbn_sc_3sub_Type1[iA][iV][iP];
				delete rbn_sc_3sub_Type2[iA][iV][iP];
				delete raw_nsc_3sub_Type1[iA][iV][iP];
				delete raw_nsc_3sub_Type2[iA][iV][iP];
				delete rbn_nsc_3sub_Type1[iA][iV][iP];
				delete rbn_nsc_3sub_Type2[iA][iV][iP];
				delete raw_ac_3sub[iA][iV][iP];
				delete rbn_ac_3sub[iA][iV][iP];
				delete raw_nac_3sub[iA][iV][iP];
				delete rbn_nac_3sub[iA][iV][iP];
			}
		}
	}
}

void Cumu::cal_all()
{
/*
	run_c_1sub();
	run_sc_1sub();
	run_ac_1sub();

	run_nc_1sub();
	run_nsc_1sub();
	run_nac_1sub();
	run_fluc_1sub();
*/

	run_c_3sub();
	run_sc_3sub();
	run_ac_3sub();

	run_nc_3sub();
	run_nsc_3sub();
	run_nac_3sub();
}

void Cumu::run_c_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			pc2_1sub_mean[iV][iP]->Divide(pc2_1sub_wght[iV][iP]);
			pc4_1sub_mean[iV][iP]->Divide(pc4_1sub_wght[iV][iP]);
			pc6_1sub_mean[iV][iP]->Divide(pc6_1sub_wght[iV][iP]);

			cal_c_1sub(pc2_1sub_mean[iV][iP], pc4_1sub_mean[iV][iP], pc6_1sub_mean[iV][iP], ori_c2_1sub[iV][iP], ori_c4_1sub[iV][iP], ori_c6_1sub[iV][iP]);

			rebin_cumu(ori_c2_1sub[iV][iP],cnt_1sub[iP],raw_c2_1sub[iV][iP]);
			rebin_cumu(ori_c4_1sub[iV][iP],cnt_1sub[iP],raw_c4_1sub[iV][iP]);
			rebin_cumu(ori_c6_1sub[iV][iP],cnt_1sub[iP],raw_c6_1sub[iV][iP]);

			rebin_cumu(raw_c2_1sub[iV][iP],raw_cnt_1sub[iP],rbn_c2_1sub[iV][iP]);
			rebin_cumu(raw_c4_1sub[iV][iP],raw_cnt_1sub[iP],rbn_c4_1sub[iV][iP]);
			rebin_cumu(raw_c6_1sub[iV][iP],raw_cnt_1sub[iP],rbn_c6_1sub[iV][iP]);
		}
	}
}

void Cumu::run_sc_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			psc4_1sub_mean[iV][iP]->Divide(psc4_1sub_wght[iV][iP]);

			if(iV==0) cal_sc(pc2_1sub_mean[1][iP], pc2_1sub_mean[2][iP], psc4_1sub_mean[iV][iP], ori_sc_1sub[iV][iP]);
			if(iV==1) cal_sc(pc2_1sub_mean[1][iP], pc2_1sub_mean[3][iP], psc4_1sub_mean[iV][iP], ori_sc_1sub[iV][iP]);
			if(iV==2) cal_sc(pc2_1sub_mean[2][iP], pc2_1sub_mean[3][iP], psc4_1sub_mean[iV][iP], ori_sc_1sub[iV][iP]);
			if(iV==3) cal_sc(pc2_1sub_mean[2][iP], pc2_1sub_mean[4][iP], psc4_1sub_mean[iV][iP], ori_sc_1sub[iV][iP]);

			rebin_cumu(ori_sc_1sub[iV][iP],cnt_1sub[iP],raw_sc_1sub[iV][iP]);
			
			rebin_cumu(raw_sc_1sub[iV][iP],raw_cnt_1sub[iP],rbn_sc_1sub[iV][iP]);
		}
	}
}

void Cumu::run_ac_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			pac3_1sub_mean[iV][iP]->Divide(pac3_1sub_wght[iV][iP]);
			
			rebin_cumu(pac3_1sub_mean[iV][iP],cnt_1sub[iP],raw_ac_1sub[iV][iP]);
			
			rebin_cumu(raw_ac_1sub[iV][iP],raw_cnt_1sub[iP],rbn_ac_1sub[iV][iP]);
		}
	}
}

void Cumu::run_c_3sub()
{
	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				pc2_1_3sub_mean[iA][iV][iP]->Divide(pc2_1_3sub_wght[iA][iV][iP]);
				pc2_2_3sub_mean[iA][iV][iP]->Divide(pc2_2_3sub_wght[iA][iV][iP]);
				pc4_3sub_mean[iA][iV][iP]->Divide(pc4_3sub_wght[iA][iV][iP]);

				cal_c_3sub(pc2_1_3sub_mean[iA][iV][iP], pc2_2_3sub_mean[iA][iV][iP], pc4_3sub_mean[iA][iV][iP], ori_c2_3sub[iA][iV][iP], ori_c4_3sub[iA][iV][iP]);

				rebin_cumu(ori_c2_3sub[iA][iV][iP],cnt_3sub[iA][iP],raw_c2_3sub[iA][iV][iP]);
				rebin_cumu(ori_c4_3sub[iA][iV][iP],cnt_3sub[iA][iP],raw_c4_3sub[iA][iV][iP]);

				rebin_cumu(raw_c2_3sub[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_c2_3sub[iA][iV][iP]);
				rebin_cumu(raw_c4_3sub[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_c4_3sub[iA][iV][iP]);
				
				rebin_cumu(pc4_3sub_mean[iA][iV][iP],cnt_3sub[iA][iP],raw_p4_3sub[iA][iV][iP]);
			}
		}
	}
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			merge3sub(rbn_c4_3sub[0][iV][iP],rbn_cnt_3sub[0][iP],rbn_c4_3sub[1][iV][iP],rbn_cnt_3sub[1][iP],rbn_c4_3sub[2][iV][iP],rbn_cnt_3sub[2][iP]);
			merge3sub(raw_p4_3sub[0][iV][iP],raw_cnt_3sub[0][iP],raw_p4_3sub[1][iV][iP],raw_cnt_3sub[1][iP],raw_p4_3sub[2][iV][iP],raw_cnt_3sub[2][iP]);
		}
	}
}

void Cumu::run_sc_3sub()
{
	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				psc4_3sub_Type1_mean[iA][iV][iP]->Divide(psc4_3sub_wght[iA][iV][iP]);

				if(iV==0) cal_sc(pc2_1_3sub_mean[iA][1][iP],pc2_2_3sub_mean[iA][2][iP],psc4_3sub_Type1_mean[iA][iV][iP],ori_sc_3sub_Type1[iA][iV][iP]);
				if(iV==1) cal_sc(pc2_1_3sub_mean[iA][1][iP],pc2_2_3sub_mean[iA][3][iP],psc4_3sub_Type1_mean[iA][iV][iP],ori_sc_3sub_Type1[iA][iV][iP]);
				if(iV==2) cal_sc(pc2_1_3sub_mean[iA][2][iP],pc2_2_3sub_mean[iA][3][iP],psc4_3sub_Type1_mean[iA][iV][iP],ori_sc_3sub_Type1[iA][iV][iP]);
				if(iV==3) cal_sc(pc2_1_3sub_mean[iA][2][iP],pc2_2_3sub_mean[iA][4][iP],psc4_3sub_Type1_mean[iA][iV][iP],ori_sc_3sub_Type1[iA][iV][iP]);

				psc4_3sub_Type2_mean[iA][iV][iP]->Divide(psc4_3sub_wght[iA][iV][iP]);

				if(iV==0) cal_sc(pc2_1_3sub_mean[iA][1][iP],pc2_2_3sub_mean[iA][2][iP],psc4_3sub_Type2_mean[iA][iV][iP],ori_sc_3sub_Type2[iA][iV][iP]);
				if(iV==1) cal_sc(pc2_1_3sub_mean[iA][1][iP],pc2_2_3sub_mean[iA][3][iP],psc4_3sub_Type2_mean[iA][iV][iP],ori_sc_3sub_Type2[iA][iV][iP]);
				if(iV==2) cal_sc(pc2_1_3sub_mean[iA][2][iP],pc2_2_3sub_mean[iA][3][iP],psc4_3sub_Type2_mean[iA][iV][iP],ori_sc_3sub_Type2[iA][iV][iP]);
				if(iV==3) cal_sc(pc2_1_3sub_mean[iA][2][iP],pc2_2_3sub_mean[iA][4][iP],psc4_3sub_Type2_mean[iA][iV][iP],ori_sc_3sub_Type2[iA][iV][iP]);

				rebin_cumu(ori_sc_3sub_Type1[iA][iV][iP],cnt_3sub[iA][iP],raw_sc_3sub_Type1[iA][iV][iP]);
				rebin_cumu(raw_sc_3sub_Type1[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_sc_3sub_Type1[iA][iV][iP]);

				rebin_cumu(ori_sc_3sub_Type2[iA][iV][iP],cnt_3sub[iA][iP],raw_sc_3sub_Type2[iA][iV][iP]);
				rebin_cumu(raw_sc_3sub_Type2[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_sc_3sub_Type2[iA][iV][iP]);
			}
		}
	}
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			//merge3sub(rbn_sc_3sub[0][iV][iP],rbn_cnt_3sub[0][iP],rbn_sc_3sub[1][iV][iP],rbn_cnt_3sub[1][iP],rbn_sc_3sub[2][iV][iP],rbn_cnt_3sub[2][iP]);
		}
	}
}

void Cumu::run_ac_3sub()
{
	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				pac3_3sub_mean[iA][iV][iP]->Divide(pac3_3sub_wght[iA][iV][iP]);

				rebin_cumu(pac3_3sub_mean[iA][iV][iP],cnt_3sub[iA][iP],raw_ac_3sub[iA][iV][iP]);
				rebin_cumu(raw_ac_3sub[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_ac_3sub[iA][iV][iP]);
			}
		}
	}
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			merge3sub(rbn_ac_3sub[0][iV][iP],rbn_cnt_3sub[0][iP],rbn_ac_3sub[1][iV][iP],rbn_cnt_3sub[1][iP],rbn_ac_3sub[2][iV][iP],rbn_cnt_3sub[2][iP]);
		}
	}
}

void Cumu::run_nc_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			cal_nc(raw_c2_3sub[1][iV][iP],raw_c4_1sub[iV][iP],raw_c6_1sub[iV][iP],raw_nc4_1sub[iV][iP],raw_nc6_1sub[iV][iP]);

			rebin_cumu(raw_nc4_1sub[iV][iP],raw_cnt_1sub[iP],rbn_nc4_1sub[iV][iP]);
			rebin_cumu(raw_nc6_1sub[iV][iP],raw_cnt_1sub[iP],rbn_nc6_1sub[iV][iP]);
		}
	}
}

void Cumu::run_nsc_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			if(iV==0) cal_nsc(raw_c2_3sub[1][1][iP], raw_c2_3sub[1][2][iP], raw_sc_1sub[iV][iP], raw_nsc_1sub[iV][iP]);
			if(iV==1) cal_nsc(raw_c2_3sub[1][1][iP], raw_c2_3sub[1][3][iP], raw_sc_1sub[iV][iP], raw_nsc_1sub[iV][iP]);
			if(iV==2) cal_nsc(raw_c2_3sub[1][2][iP], raw_c2_3sub[1][3][iP], raw_sc_1sub[iV][iP], raw_nsc_1sub[iV][iP]);
			if(iV==3) cal_nsc(raw_c2_3sub[1][2][iP], raw_c2_3sub[1][4][iP], raw_sc_1sub[iV][iP], raw_nsc_1sub[iV][iP]);

			rebin_cumu(raw_nsc_1sub[iV][iP],raw_cnt_1sub[iP],rbn_nsc_1sub[iV][iP]);
		}
	}
}

void Cumu::run_nac_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			if(iV==0) cal_nac(raw_c4_3sub[0][1][iP],raw_c2_3sub[1][1][iP],raw_c2_3sub[1][2][iP],raw_ac_1sub[iV][iP],raw_nac_1sub[iV][iP]);
			if(iV==1) cal_nac(raw_c4_3sub[0][1][iP],raw_c2_3sub[1][1][iP],raw_c2_3sub[1][3][iP],raw_ac_1sub[iV][iP],raw_nac_1sub[iV][iP]);
			if(iV==2) cal_nac(raw_c4_3sub[0][2][iP],raw_c2_3sub[1][2][iP],raw_c2_3sub[1][4][iP],raw_ac_1sub[iV][iP],raw_nac_1sub[iV][iP]);

			rebin_cumu(raw_nac_1sub[iV][iP],raw_cnt_1sub[iP],rbn_nac_1sub[iV][iP]);
		}
	}
}

void Cumu::run_fluc_1sub()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			cal_fluc(raw_c2_3sub[1][iV][iP],raw_c4_1sub[iV][iP],raw_c6_1sub[iV][iP],raw_isGauss_1sub[iV][iP],raw_isPower_1sub[iV][iP]);
			
			rebin_cumu(raw_isGauss_1sub[iV][iP],raw_cnt_1sub[iP],rbn_isGauss_1sub[iV][iP]);
			rebin_cumu(raw_isPower_1sub[iV][iP],raw_cnt_1sub[iP],rbn_isPower_1sub[iV][iP]);
		}
	}
}

void Cumu::run_nc_3sub()
{
	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				cal_nc(raw_c2_3sub[1][iV][iP],raw_c4_3sub[iA][iV][iP],raw_c6_1sub[iV][iP],raw_nc4_3sub[iA][iV][iP],raw_nc6_3sub[iA][iV][iP]);

				rebin_cumu(raw_nc4_3sub[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_nc4_3sub[iA][iV][iP]);
				rebin_cumu(raw_nc6_3sub[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_nc6_3sub[iA][iV][iP]);
			}
		}
	}
}

void Cumu::run_nsc_3sub()
{
	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				if(iV==0) cal_nsc(raw_c2_3sub[1][1][iP], raw_c2_3sub[1][2][iP], raw_sc_3sub_Type1[iA][iV][iP], raw_nsc_3sub_Type1[iA][iV][iP]);
				if(iV==1) cal_nsc(raw_c2_3sub[1][1][iP], raw_c2_3sub[1][3][iP], raw_sc_3sub_Type1[iA][iV][iP], raw_nsc_3sub_Type1[iA][iV][iP]);
				if(iV==2) cal_nsc(raw_c2_3sub[1][2][iP], raw_c2_3sub[1][3][iP], raw_sc_3sub_Type1[iA][iV][iP], raw_nsc_3sub_Type1[iA][iV][iP]);
				if(iV==3) cal_nsc(raw_c2_3sub[1][2][iP], raw_c2_3sub[1][4][iP], raw_sc_3sub_Type1[iA][iV][iP], raw_nsc_3sub_Type1[iA][iV][iP]);

				rebin_cumu(raw_nsc_3sub_Type1[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_nsc_3sub_Type1[iA][iV][iP]);
				
				if(iV==0) cal_nsc(raw_c2_3sub[1][1][iP], raw_c2_3sub[1][2][iP], raw_sc_3sub_Type2[iA][iV][iP], raw_nsc_3sub_Type2[iA][iV][iP]);
				if(iV==1) cal_nsc(raw_c2_3sub[1][1][iP], raw_c2_3sub[1][3][iP], raw_sc_3sub_Type2[iA][iV][iP], raw_nsc_3sub_Type2[iA][iV][iP]);
				if(iV==2) cal_nsc(raw_c2_3sub[1][2][iP], raw_c2_3sub[1][3][iP], raw_sc_3sub_Type2[iA][iV][iP], raw_nsc_3sub_Type2[iA][iV][iP]);
				if(iV==3) cal_nsc(raw_c2_3sub[1][2][iP], raw_c2_3sub[1][4][iP], raw_sc_3sub_Type2[iA][iV][iP], raw_nsc_3sub_Type2[iA][iV][iP]);

				rebin_cumu(raw_nsc_3sub_Type2[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_nsc_3sub_Type2[iA][iV][iP]);
			}
		}
	}
}

void Cumu::run_nac_3sub()
{
	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				if(iV==0) cal_nac(raw_c4_3sub[0][1][iP],raw_c2_3sub[1][1][iP],raw_c2_3sub[1][2][iP],raw_ac_3sub[iA][iV][iP],raw_nac_3sub[iA][iV][iP]);
				if(iV==1) cal_nac(raw_c4_3sub[0][1][iP],raw_c2_3sub[1][1][iP],raw_c2_3sub[1][3][iP],raw_ac_3sub[iA][iV][iP],raw_nac_3sub[iA][iV][iP]);
				if(iV==2) cal_nac(raw_c4_3sub[0][2][iP],raw_c2_3sub[1][2][iP],raw_c2_3sub[1][4][iP],raw_ac_3sub[iA][iV][iP],raw_nac_3sub[iA][iV][iP]);

				rebin_cumu(raw_nac_3sub[iA][iV][iP],raw_cnt_3sub[iA][iP],rbn_nac_3sub[iA][iV][iP]);
			}
		}
	}
}

void Cumu::cal_c_1sub(TH1D* pc2, TH1D* pc4, TH1D* pc6, TH1D* c2, TH1D* c4, TH1D* c6)
{
	int NX = c2->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p2 = pc2->GetBinContent(iX+1);
		double p4 = pc4->GetBinContent(iX+1);
		double p6 = pc6->GetBinContent(iX+1);
		c2->SetBinContent(iX+1,p2);
		c4->SetBinContent(iX+1,p4-2.*pow(p2,2));
		c6->SetBinContent(iX+1,p6-9.*p4*p2+12.*pow(p2,3));
	}
}

void Cumu::cal_sc(TH1D* pc2_1, TH1D* pc2_2, TH1D* psc4, TH1D* sc4)
{
	int NX = pc2_1->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p2_1 = pc2_1->GetBinContent(iX+1);
		double p2_2 = pc2_2->GetBinContent(iX+1);
		double p4 = psc4->GetBinContent(iX+1);
		sc4->SetBinContent(iX+1,p4-p2_1*p2_2);
	}
}

void Cumu::cal_c_3sub(TH1D* pc2_1, TH1D* pc2_2, TH1D* pc4, TH1D* c2, TH1D* c4)
{
	int NX = pc2_1->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p2_1 = pc2_1->GetBinContent(iX+1);
		double p2_2 = pc2_2->GetBinContent(iX+1);
		double p4 = pc4->GetBinContent(iX+1);
		c2->SetBinContent(iX+1,p2_1);
		c4->SetBinContent(iX+1,p4-2.*p2_1*p2_2);
	}
}

void Cumu::cal_nc(TH1D* c2, TH1D* c4, TH1D* c6, TH1D* nc4, TH1D* nc6)
{
	int NX = c2->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p2 = c2->GetBinContent(iX+1);
		double p4 = c4->GetBinContent(iX+1);
		double p6 = c6->GetBinContent(iX+1);
		if(p2==0) continue;
		double pc4 = p4/pow(p2,2);
		double pc6 = p6/pow(p2,3)/4;
		nc4->SetBinContent(iX+1,pc4);
		nc6->SetBinContent(iX+1,pc6);
	}
}

void Cumu::cal_nsc(TH1D* c2_1, TH1D* c2_2, TH1D* sc4, TH1D* nsc4)
{
	int NX = c2_1->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p2_1 = c2_1->GetBinContent(iX+1);
		double p2_2 = c2_2->GetBinContent(iX+1);
		double p4 = sc4->GetBinContent(iX+1);
		if(p2_1!=0 && p2_2!=0) nsc4->SetBinContent(iX+1,p4/p2_1/p2_2);
	}
}

void Cumu::cal_nac(TH1D* c4, TH1D* c2_1, TH1D* c2_2, TH1D* ac3, TH1D* nac3)
{
	int NX = c4->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p4 = c4->GetBinContent(iX+1);
		double p2_1 = c2_1->GetBinContent(iX+1);
		double p2_2 = c2_2->GetBinContent(iX+1);
		double pa3 = ac3->GetBinContent(iX+1);
		p4 += 2.*pow(p2_1,2);
		if(p4*p2_2>0) nac3->SetBinContent(iX+1,pa3/sqrt(p4*p2_2));
		else if(p4*p2_2<0) nac3->SetBinContent(iX+1,pa3/sqrt(-p4*p2_2));
	}
}

void Cumu::cal_fluc(TH1D* c2, TH1D* c4, TH1D* c6, TH1D* isGauss, TH1D* isPower)
{
	int NX = c2->GetNbinsX();
	for(int iX=0; iX<NX; iX++)
	{
		double p2 = c2->GetBinContent(iX+1);
		double p4 = c4->GetBinContent(iX+1);
		double p6 = c6->GetBinContent(iX+1);
		double gauss = 0;
		double power = 0;
		if(p4!=0)
		{
			if(p4<0) gauss = p6/(4*pow(-p4,3./2));
			else gauss = -p6/(4*pow(p4,3./2));
		}
		if(p2!=0 && p4!=0)
		{
			power = p6*(1-2*pow(p2,2)/p4)/(12*-p4*p2);
		}
		isGauss->SetBinContent(iX+1,gauss);
		isPower->SetBinContent(iX+1,power);
	}
}

void Cumu::rebin_cumu(TH1D* hIn, TH1D* hWght, TH1D* hOut)
{
	unsigned int nBins = hIn->GetNbinsX();
	unsigned int nRebins = hOut->GetNbinsX();
	double sum;
	double cnt;
	for(unsigned int iB=0; iB<nRebins; iB++)
	{
		sum = 0;
		cnt = 0;
		for(unsigned int i=0; i<nBins; i++)
		{
			double xLw = hOut->GetBinCenter(iB+1)-0.5*hOut->GetBinWidth(iB+1);
			double xUp = hOut->GetBinCenter(iB+1)+0.5*hOut->GetBinWidth(iB+1);
			if(hIn->GetBinCenter(i+1)<xLw || hIn->GetBinCenter(i+1)>=xUp) continue;
			sum += hWght->GetBinContent(i+1)*hIn->GetBinContent(i+1);
			cnt += hWght->GetBinContent(i+1);
		}
		if(cnt>0) hOut->SetBinContent(iB+1,sum/cnt);
		else hOut->SetBinContent(iB+1,0);
	}
}

void Cumu::merge3sub(TH1D* h1, TH1D* h1w, TH1D* h2, TH1D* h2w, TH1D* h3, TH1D* h3w)
{
	unsigned int nBins = h1->GetNbinsX();
	double sum;
	double cnt;
	for(unsigned int iB=0; iB<nBins; iB++)
	{
		sum = 0;
		cnt = 0;

		sum += h1w->GetBinContent(iB+1)*h1->GetBinContent(iB+1);
		sum += h2w->GetBinContent(iB+1)*h2->GetBinContent(iB+1);
		sum += h3w->GetBinContent(iB+1)*h3->GetBinContent(iB+1);
		cnt += h1w->GetBinContent(iB+1);
		cnt += h2w->GetBinContent(iB+1);
		cnt += h3w->GetBinContent(iB+1);

		if(cnt>0) h1->SetBinContent(iB+1,sum/cnt);
		else h1->SetBinContent(iB+1,0);
	}
}

void Cumu::writeHist(TFile*& fIn)
{
  fIn->cd();
	/*
	for(unsigned int iP=0; iP<NP; iP++)
	{
		cnt_1sub[iP]->Write();
		raw_cnt_1sub[iP]->Write();
		rbn_cnt_1sub[iP]->Write();
	}
	*/
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			ori_c2_1sub[iV][iP]->Write();
			raw_c2_1sub[iV][iP]->Write();
			rbn_c2_1sub[iV][iP]->Write();
			ori_c4_1sub[iV][iP]->Write();
			raw_c4_1sub[iV][iP]->Write();
			rbn_c4_1sub[iV][iP]->Write();
			ori_c6_1sub[iV][iP]->Write();
			raw_c6_1sub[iV][iP]->Write();
			rbn_c6_1sub[iV][iP]->Write();
			raw_nc4_1sub[iV][iP]->Write();
			rbn_nc4_1sub[iV][iP]->Write();
			raw_nc6_1sub[iV][iP]->Write();
			rbn_nc6_1sub[iV][iP]->Write();
			ori_sc_1sub[iV][iP]->Write();
			raw_sc_1sub[iV][iP]->Write();
			rbn_sc_1sub[iV][iP]->Write();
			raw_nsc_1sub[iV][iP]->Write();
			rbn_nsc_1sub[iV][iP]->Write();
			raw_ac_1sub[iV][iP]->Write();
			rbn_ac_1sub[iV][iP]->Write();
			raw_nac_1sub[iV][iP]->Write();
			rbn_nac_1sub[iV][iP]->Write();
			raw_isGauss_1sub[iV][iP]->Write();
			rbn_isGauss_1sub[iV][iP]->Write();
			raw_isPower_1sub[iV][iP]->Write();
			rbn_isPower_1sub[iV][iP]->Write();
		}
	}

	for(unsigned int iA=0; iA<NA; iA++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				ori_c2_3sub[iA][iV][iP]->Write();
				raw_c2_3sub[iA][iV][iP]->Write();
				rbn_c2_3sub[iA][iV][iP]->Write();
				ori_c4_3sub[iA][iV][iP]->Write();
				raw_c4_3sub[iA][iV][iP]->Write();
				rbn_c4_3sub[iA][iV][iP]->Write();
				raw_nc4_3sub[iA][iV][iP]->Write();
				rbn_nc4_3sub[iA][iV][iP]->Write();
				//raw_nc6_3sub[iA][iV][iP]->Write();
				//rbn_nc6_3sub[iA][iV][iP]->Write();
				ori_sc_3sub_Type1[iA][iV][iP]->Write();
				ori_sc_3sub_Type2[iA][iV][iP]->Write();
				raw_sc_3sub_Type1[iA][iV][iP]->Write();
				raw_sc_3sub_Type2[iA][iV][iP]->Write();
				rbn_sc_3sub_Type1[iA][iV][iP]->Write();
				rbn_sc_3sub_Type2[iA][iV][iP]->Write();
				raw_nsc_3sub_Type1[iA][iV][iP]->Write();
				raw_nsc_3sub_Type2[iA][iV][iP]->Write();
				rbn_nsc_3sub_Type1[iA][iV][iP]->Write();
				rbn_nsc_3sub_Type2[iA][iV][iP]->Write();
				raw_ac_3sub[iA][iV][iP]->Write();
				rbn_ac_3sub[iA][iV][iP]->Write();
				raw_nac_3sub[iA][iV][iP]->Write();
				rbn_nac_3sub[iA][iV][iP]->Write();
			}
		}
	}
  fIn->Close();
}

void Cumu::readHist_VP(TFile* fIn, TH1D*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TH1D*)fIn->Get(name);
}

void Cumu::readHist_AVP(TFile* fIn, TH1D*& hIn, const char* hName, int iA, int iV, int iP)
{
	sprintf(name,"%s_Pm%d_Har%d_Pt%d",hName,iA,iV,iP);
	hIn = (TH1D*)fIn->Get(name);
}

void Cumu::initHist_VP(TH1D*& hIn, const char* hName, int iV, int iP, int nBins)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = new TH1D(name,"",nBins,minBin,maxBin);
}

void Cumu::initHist_AVP(TH1D*& hIn, const char* hName, int iA, int iV, int iP, int nBins)
{
	sprintf(name,"%s_Pm%d_Har%d_Pt%d",hName,iA,iV,iP);
	hIn = new TH1D(name,"",nBins,minBin,maxBin);
}

