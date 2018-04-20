#include "Tool.h"



Tool::Tool()
{
	TFile* fInEff = new TFile("INPUT/trkEff_LOOSE.root","READ");
	allTrue = (TH3D*)fInEff->Get("allTrue");
	allReco = (TH3D*)fInEff->Get("allReco");
	matchedReco = (TH3D*)fInEff->Get("matchedReco");
	unMatchedReco = (TH3D*)fInEff->Get("unMatchedReco");
	//matchedReco->Divide(allTrue);
	//unMatchedReco->Divide(allReco);
	

	for(unsigned int iC=0; iC<nCent; iC++)
	{
		for(unsigned int iP=0; iP<nPt; iP++)
		{
			sprintf(name,"hEtaEff_Cent%d_Pt%d",iC,iP); hEtaEff[iC][iP] = new TH1D(name,"",50,-2.5,2.5);
			sprintf(name,"hEtaFak_Cent%d_Pt%d",iC,iP); hEtaFak[iC][iP] = new TH1D(name,"",50,-2.5,2.5);
		}
	}

	for(unsigned int iC=0; iC<nCent; iC++)
	{
		for(unsigned int iP=0; iP<nPt; iP++)
		{
			cout<<iC+1<<"/"<<nCent<<" | "<<iP+1<<"/"<<nPt<<endl;
			rebin_ratio(allTrue, matchedReco, hEtaEff[iC][iP], iC, iP);
			rebin_ratio(allReco, unMatchedReco, hEtaFak[iC][iP], iC, iP);
		}
	}

	TFile* fOut = new TFile("trkEff_mon_LOOSE.root","RECREATE");
	fOut->cd();
	for(unsigned int iC=0; iC<nCent; iC++)
	{
		for(unsigned int iP=0; iP<nPt; iP++)
		{
			hEtaEff[iC][iP]->Write();
			hEtaFak[iC][iP]->Write();
		}
	}
	fOut->Close();
}

Tool::~Tool()
{
}

void Tool::rebin_ratio(TH3D* hRef, TH3D* hChk, TH1D* hRat, int iC, int iP)
{
	unsigned int NX = hRef->GetNbinsX();
	unsigned int NY = hRef->GetNbinsY();
	unsigned int NZ = hRef->GetNbinsZ();
	unsigned int NB = hRat->GetNbinsX();
	for(unsigned int iB=0; iB<NB; iB++)
	{
		double cnt = 0;
		double sum = 0;
		double t_xLw = cutCent[80-int(binCent[iC+1])];
		double t_xUp = cutCent[80-int(binCent[iC])];
		double t_yLw = binPt[iP];
		double t_yUp = binPt[iP+1];
		double t_zLw = hRat->GetBinLowEdge(iB+1);
		double t_zUp = hRat->GetBinLowEdge(iB+1) + hRat->GetBinWidth(iB+1);
		for(unsigned int iX=0; iX<NX; iX++)
		{
			for(unsigned int iY=0; iY<NY; iY++)
			{
				//for(unsigned int iZ=0; iZ<NZ; iZ++)
				{
					double xLw = hRef->GetXaxis()->GetBinLowEdge(iX+1);
					double xUp = hRef->GetXaxis()->GetBinLowEdge(iX+1) + hRef->GetXaxis()->GetBinWidth(iX+1);
					double yLw = hRef->GetYaxis()->GetBinLowEdge(iY+1);
					double yUp = hRef->GetYaxis()->GetBinLowEdge(iY+1) + hRef->GetYaxis()->GetBinWidth(iY+1);
					//double zLw = hRef->GetZaxis()->GetBinLowEdge(iZ+1);
					//double zUp = hRef->GetZaxis()->GetBinLowEdge(iZ+1) + hRef->GetZaxis()->GetBinWidth(iZ+1);

					if(xLw<t_xLw || xUp>t_xUp) continue;
					if(yLw<t_yLw-1e-3 || yUp>t_yUp+1e-3) continue;
					//if(zLw<t_zLw || zUp>t_zUp) continue;

					//if(iB==0) cout<<t_xLw<<" | "<<t_xUp<<" ||| "<<xLw<<" | "<<xUp<<endl;
					//if(iB==0) cout<<t_yLw<<" | "<<t_yUp<<" ||| "<<yLw<<" | "<<yUp<<endl;

					cnt += hChk->GetBinContent(iX+1,iY+1,iB+1);
					sum += hRef->GetBinContent(iX+1,iY+1,iB+1);
				}
			}
		}
		//cout<<(t_zLw+t_zUp)/2<<": "<<cnt<<"/"<<sum<<endl;
		if(sum>0) hRat->SetBinContent(iB+1,cnt/sum);
	}
}





