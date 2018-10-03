#include "Phase3.h"

Phase3::Phase3(unsigned int iBin)
{
	initialize(iBin);
	execute();
	finalize(iBin);
}

Phase3::~Phase3()
{
}

void Phase3::execute()
{
	cout<<"execute..."<<endl;

	vector<TH1D*> hSamples;
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_c2_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],c2_1sub[iV][iP],"c2_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_c4_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],c4_1sub[iV][iP],"c4_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_c6_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],c6_1sub[iV][iP],"c6_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nc4_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],nc4_1sub[iV][iP],"nc4_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nc6_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],nc6_1sub[iV][iP],"nc6_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_sc_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],sc_1sub[iV][iP],"sc_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nsc_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],nsc_1sub[iV][iP],"nsc_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_ac_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],ac_1sub[iV][iP],"ac_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nac_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],nac_1sub[iV][iP],"nac_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_isGauss_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],isGauss_1sub[iV][iP],"isGauss_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_isPower_1sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_1sub[iP],isPower_1sub[iV][iP],"isPower_1sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_c2_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],c2_3sub[iV][iP],"c2_3sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_c4_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],c4_3sub[iV][iP],"c4_3sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nc4_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],nc4_3sub[iV][iP],"nc4_3sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_sc_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],sc_3sub[iV][iP],"sc_3sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nsc_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],nsc_3sub[iV][iP],"nsc_3sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_ac_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],ac_3sub[iV][iP],"ac_3sub",iV,iP); hSamples.clear();

			for(unsigned int iSample=0; iSample<=nSample; iSample++) hSamples.push_back(rbn_nac_3sub[iSample][iV][iP]);
			cal_sts(hSamples,xMean_3sub[iP],nac_3sub[iV][iP],"nac_3sub",iV,iP); hSamples.clear();
		}
	}
}

void Phase3::initialize(unsigned int iBin)
{
	cout<<"initialize..."<<endl;

	TH1::AddDirectory(kFALSE);
	TFile* fIn[nSample+1];
	for(unsigned int iSample=0; iSample<=nSample; iSample++)
	{
		if(iSample<nSample) sprintf(name,"../OUTPUT/Phase2/bin%d/Phase2_Sample%d.root",iBin,iSample);
		else sprintf(name,"../OUTPUT/Phase2/bin%d/Phase2_All.root",iBin);
		fIn[iSample] = new TFile(name,"READ");
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				readHist_VP(fIn[iSample],rbn_c2_1sub[iSample][iV][iP],"rbn_c2_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_c4_1sub[iSample][iV][iP],"rbn_c4_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_c6_1sub[iSample][iV][iP],"rbn_c6_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_nc4_1sub[iSample][iV][iP],"rbn_nc4_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_nc6_1sub[iSample][iV][iP],"rbn_nc6_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_sc_1sub[iSample][iV][iP],"rbn_sc_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_nsc_1sub[iSample][iV][iP],"rbn_nsc_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_ac_1sub[iSample][iV][iP],"rbn_ac_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_nac_1sub[iSample][iV][iP],"rbn_nac_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_isGauss_1sub[iSample][iV][iP],"rbn_isGauss_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_isPower_1sub[iSample][iV][iP],"rbn_isPower_1sub",iV,iP);
				readHist_VP(fIn[iSample],rbn_c2_3sub[iSample][iV][iP],"rbn_c2_3sub_Pm1",iV,iP);
				readHist_VP(fIn[iSample],rbn_c4_3sub[iSample][iV][iP],"rbn_c4_3sub_Pm2",iV,iP);
				readHist_VP(fIn[iSample],rbn_nc4_3sub[iSample][iV][iP],"rbn_nc4_3sub_Pm2",iV,iP);
				readHist_VP(fIn[iSample],rbn_sc_3sub[iSample][iV][iP],"rbn_sc_3sub_Pm2",iV,iP);
				readHist_VP(fIn[iSample],rbn_nsc_3sub[iSample][iV][iP],"rbn_nsc_3sub_Pm2",iV,iP);
				readHist_VP(fIn[iSample],rbn_ac_3sub[iSample][iV][iP],"rbn_ac_3sub_Pm2",iV,iP);
				readHist_VP(fIn[iSample],rbn_nac_3sub[iSample][iV][iP],"rbn_nac_3sub_Pm2",iV,iP);
			}
		}
	}

	TFile* fX = new TFile("../../MAIN_binFCal/OUTPUT/Phase1/Phase1_All.root","READ");
	for(unsigned int iP=0; iP<NP; iP++)
	{
		sprintf(name,"cnt_1sub_Pt%d",iP); cnt_1sub[iP] = (TH1D*)fX->Get(name);
		sprintf(name,"xMean_1sub_Pt%d",iP); xMean_1sub[iP] = new TH1D(name,"",nRebin[iBin],minBin,maxBin);
		for(unsigned int iB=0; iB<nRebin[iBin]; iB++)
		{
			double sum = 0;
			double cnt = 0;
			for(unsigned int i=0; i<nBin; i++)
			{
				double xLw = xMean_1sub[iP]->GetBinCenter(iB+1)-0.5*xMean_1sub[iP]->GetBinWidth(iB+1);
				double xUp = xMean_1sub[iP]->GetBinCenter(iB+1)+0.5*xMean_1sub[iP]->GetBinWidth(iB+1);
				if(cnt_1sub[iP]->GetBinCenter(i+1)<xLw || cnt_1sub[iP]->GetBinCenter(i+1)>=xUp) continue;
				sum += cnt_1sub[iP]->GetBinCenter(i+1)*cnt_1sub[iP]->GetBinContent(i+1);
				cnt += cnt_1sub[iP]->GetBinContent(i+1);
			}
			if(cnt>0) xMean_1sub[iP]->SetBinContent(iB+1,sum/cnt);
			else xMean_1sub[iP]->SetBinContent(iB+1,0);
		}
		sprintf(name,"cnt_3sub_Pm0_Pt%d",iP); cnt_3sub[iP] = (TH1D*)fX->Get(name);
		sprintf(name,"xMean_3sub_Pt%d",iP); xMean_3sub[iP] = new TH1D(name,"",nRebin[iBin],minBin,maxBin);
		for(unsigned int iB=0; iB<nRebin[iBin]; iB++)
		{
			double sum = 0;
			double cnt = 0;
			for(unsigned int i=0; i<nBin; i++)
			{
				double xLw = xMean_3sub[iP]->GetBinCenter(iB+1)-0.5*xMean_3sub[iP]->GetBinWidth(iB+1);
				double xUp = xMean_3sub[iP]->GetBinCenter(iB+1)+0.5*xMean_3sub[iP]->GetBinWidth(iB+1);
				if(cnt_3sub[iP]->GetBinCenter(i+1)<xLw || cnt_3sub[iP]->GetBinCenter(i+1)>=xUp) continue;
				sum += cnt_3sub[iP]->GetBinCenter(i+1)*cnt_3sub[iP]->GetBinContent(i+1);
				cnt += cnt_3sub[iP]->GetBinContent(i+1);
			}
			if(cnt>0) xMean_3sub[iP]->SetBinContent(iB+1,sum/cnt);
			else xMean_3sub[iP]->SetBinContent(iB+1,0);
		}
	}

}

void Phase3::finalize(unsigned int iBin)
{
	cout<<"finalize..."<<endl;

	sprintf(name,"../OUTPUT/Phase3/Phase3_bin%d.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			cleanPoint(c2_1sub[iV][iP]); c2_1sub[iV][iP]->Write();
			cleanPoint(c4_1sub[iV][iP]); c4_1sub[iV][iP]->Write();
			cleanPoint(c6_1sub[iV][iP]); c6_1sub[iV][iP]->Write();
			cleanPoint(nc4_1sub[iV][iP]); nc4_1sub[iV][iP]->Write();
			cleanPoint(nc6_1sub[iV][iP]); nc6_1sub[iV][iP]->Write();
			cleanPoint(sc_1sub[iV][iP]); sc_1sub[iV][iP]->Write();
			cleanPoint(nsc_1sub[iV][iP]); nsc_1sub[iV][iP]->Write();
			cleanPoint(ac_1sub[iV][iP]); ac_1sub[iV][iP]->Write();
			cleanPoint(nac_1sub[iV][iP]); nac_1sub[iV][iP]->Write();
			cleanPoint(isGauss_1sub[iV][iP]); isGauss_1sub[iV][iP]->Write();
			cleanPoint(isPower_1sub[iV][iP]); isPower_1sub[iV][iP]->Write();
			cleanPoint(c2_3sub[iV][iP]); c2_3sub[iV][iP]->Write();
			cleanPoint(c4_3sub[iV][iP]); c4_3sub[iV][iP]->Write();
			cleanPoint(nc4_3sub[iV][iP]); nc4_3sub[iV][iP]->Write();
			cleanPoint(sc_3sub[iV][iP]); sc_3sub[iV][iP]->Write();
			cleanPoint(nsc_3sub[iV][iP]); nsc_3sub[iV][iP]->Write();
			cleanPoint(ac_3sub[iV][iP]); ac_3sub[iV][iP]->Write();
			cleanPoint(nac_3sub[iV][iP]); nac_3sub[iV][iP]->Write();
		}
	}
	fOut->Close();
}

void Phase3::cal_sts(vector<TH1D*> vIn, TH1D* xMean, TGraphErrors*& gOut, const char* gName, int iV, int iP)
{
	unsigned int NB = (vIn.at(0))->GetNbinsX();
	double x[NB];
	double y[NB];
	double xErr[NB];
	double yErr[NB];
	for(unsigned int iB=0; iB<NB; iB++)
	{
		x[iB] = xMean->GetBinContent(iB+1);
		y[iB] = (vIn.at(nSample))->GetBinContent(iB+1);
		xErr[iB] = 0;
		double m1 = 0;
		double m2 = 0;
		for(unsigned int iS=0; iS<nSample; iS++)
		{
			m1 += pow((vIn.at(iS))->GetBinContent(iB+1),1);
			m2 += pow((vIn.at(iS))->GetBinContent(iB+1),2);
		}
		double mean = m1/nSample;
		double sigm = sqrt((m2-pow(m1,2)/nSample)/(nSample-1));
		double nS = 0;
		m1 = 0;
		m2 = 0;
		for(unsigned int iS=0; iS<nSample; iS++)
		{
			if(fabs((vIn.at(iS))->GetBinContent(iB+1)-mean)>2.5*sigm && fabs(mean)>1E-7) continue; // remove outliers
			m1 += pow((vIn.at(iS))->GetBinContent(iB+1),1);
			m2 += pow((vIn.at(iS))->GetBinContent(iB+1),2);
			nS ++;
		}
		if(nS<nSample-2 && iV!=0) cout<<iV<<" | "<<iP<<" || "<<iB<<" | "<<mean<<" | "<<sigm<<" | "<<(vIn.at(nSample))->GetBinContent(iB+1)<<": "<<nS<<endl;
		yErr[iB] = sqrt((m2-pow(m1,2)/nS)/nS/(nS-1));
	}
	gOut = new TGraphErrors(NB,x,y,xErr,yErr);
	sprintf(name,"%s_Har%d_Pt%d",gName,iV,iP);
	gOut->SetName(name);
}

void Phase3::cleanPoint(TGraphErrors* gIn)
{
	for(int i=gIn->GetN(); i>=0; i--)
	{
		double x;
		double y;
		gIn->GetPoint(i,x,y);
		if(x<minBin || x>maxBin) gIn->RemovePoint(i);
	}
}

void Phase3::readHist_VP(TFile* fIn, TH1D*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TH1D*)fIn->Get(name);
}


