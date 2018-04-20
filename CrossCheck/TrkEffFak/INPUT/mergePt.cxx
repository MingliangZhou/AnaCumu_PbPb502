#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

void mergePt()
{
	TFile* fIn = new TFile("looseEff4M.root","READ");

	TH3D* allTrue = (TH3D*)fIn->Get("allTrue"); allTrue->SetName("old_allTrue");
	TH3D* allReco = (TH3D*)fIn->Get("allReco"); allReco->SetName("old_allReco");
	TH3D* unmatchedRecoToPrimary = (TH3D*)fIn->Get("unmatchedRecoToPrimary");
	TH3D* matchedRecoNorm = (TH3D*)fIn->Get("matchedRecoTruePt");

	unsigned int NX = allTrue->GetNbinsX();
	unsigned int NY = allTrue->GetNbinsY();
	unsigned int NZ = allTrue->GetNbinsZ();

	double binX[NX+1];
	double binY[19] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10., 20.}; // FIX
	double binZ[NZ+1];
	
	for(unsigned int iX=0; iX<NX; iX++)
	{
		binX[iX] = allTrue->GetXaxis()->GetBinLowEdge(iX+1);
		if(iX==NX-1) binX[iX+1] = allTrue->GetXaxis()->GetBinUpEdge(iX+1);
	}
	for(unsigned int iZ=0; iZ<NZ; iZ++)
	{
		binZ[iZ] = allTrue->GetZaxis()->GetBinLowEdge(iZ+1);
		if(iZ==NZ-1) binZ[iZ+1] = allTrue->GetZaxis()->GetBinUpEdge(iZ+1);
	}

	TH3D* new_allTrue = new TH3D("allTrue","",NX,binX,18,binY,NZ,binZ);
	TH3D* new_allReco = new TH3D("allReco","",NX,binX,18,binY,NZ,binZ);
	TH3D* new_matchedRecoNorm = new TH3D("matchedReco","",NX,binX,18,binY,NZ,binZ);
	TH3D* new_unmatchedRecoToPrimary = new TH3D("unMatchedReco","",NX,binX,18,binY,NZ,binZ);

	// allTrue
	for(unsigned iX=0; iX<NX; iX++)
	{
		for(unsigned int iZ=0; iZ<NZ; iZ++)
		{
			for(unsigned int iB=0; iB<18; iB++)
			{
				double sum = 0;
				for(unsigned int iY=0; iY<NY; iY++)
				{
					if(allTrue->GetYaxis()->GetBinCenter(iY+1)>=binY[iB] && allTrue->GetYaxis()->GetBinCenter(iY+1)<binY[iB+1])
					{
						sum += allTrue->GetBinContent(iX+1,iY+1,iZ+1);
					}
				}
				new_allTrue->SetBinContent(iX+1,iB+1,iZ+1,sum);
			}
		}
	}

	// allReco
	for(unsigned iX=0; iX<NX; iX++)
	{
		for(unsigned int iZ=0; iZ<NZ; iZ++)
		{
			for(unsigned int iB=0; iB<18; iB++)
			{
				double sum = 0;
				for(unsigned int iY=0; iY<NY; iY++)
				{
					if(allReco->GetYaxis()->GetBinCenter(iY+1)>=binY[iB] && allReco->GetYaxis()->GetBinCenter(iY+1)<binY[iB+1])
					{
						sum += allReco->GetBinContent(iX+1,iY+1,iZ+1);
					}
				}
				new_allReco->SetBinContent(iX+1,iB+1,iZ+1,sum);
			}
		}
	}

	// matchedRecoNorm
	for(unsigned iX=0; iX<NX; iX++)
	{
		for(unsigned int iZ=0; iZ<NZ; iZ++)
		{
			for(unsigned int iB=0; iB<18; iB++)
			{
				double sum = 0;
				for(unsigned int iY=0; iY<NY; iY++)
				{
					if(matchedRecoNorm->GetYaxis()->GetBinCenter(iY+1)>=binY[iB] && matchedRecoNorm->GetYaxis()->GetBinCenter(iY+1)<binY[iB+1])
					{
						sum += matchedRecoNorm->GetBinContent(iX+1,iY+1,iZ+1);
					}
				}
				new_matchedRecoNorm->SetBinContent(iX+1,iB+1,iZ+1,sum);
			}
		}
	}

	// unmatchedRecoToPrimary
	for(unsigned iX=0; iX<NX; iX++)
	{
		for(unsigned int iZ=0; iZ<NZ; iZ++)
		{
			for(unsigned int iB=0; iB<18; iB++)
			{
				double sum = 0;
				for(unsigned int iY=0; iY<NY; iY++)
				{
					if(unmatchedRecoToPrimary->GetYaxis()->GetBinCenter(iY+1)>=binY[iB] && unmatchedRecoToPrimary->GetYaxis()->GetBinCenter(iY+1)<binY[iB+1])
					{
						sum += unmatchedRecoToPrimary->GetBinContent(iX+1,iY+1,iZ+1);
					}
				}
				new_unmatchedRecoToPrimary->SetBinContent(iX+1,iB+1,iZ+1,sum);
			}
		}
	}

	TFile* fOut = new TFile("trkEff_LOOSE.root","RECREATE");
	fOut->cd();
	new_allTrue->Write();
	new_allReco->Write();
	new_matchedRecoNorm->Write();
	new_unmatchedRecoToPrimary->Write();
	fOut->Close();

	TH1D* h0 = (TH1D*)new_allReco->ProjectionZ("h0");
	TH1D* h1 = (TH1D*)new_unmatchedRecoToPrimary->ProjectionZ("h1");
	h1->Divide(h0);
	h1->Draw();
}
