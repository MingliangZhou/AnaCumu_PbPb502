#include "Plot.h"

Plot::Plot(unsigned int iBin)
{
	initialize(iBin);
	execute(iBin);
	finalize();
}

Plot::~Plot()
{
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
		}
	}
}

void Plot::execute(unsigned int iBin)
{
	cout<<"execute..."<<endl;


	sprintf(name,"../PLOT/bin%d/canvas.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();

	// 0: c4_1sub
	for(unsigned int iS=1; iS<NS; iS++)
	{
		for(unsigned int iV=1; iV<4; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				draw_sepa(c4_1sub[1][iS][iV][iP],c4_1sub[0][iS][iV][iP],ratio_c4_1sub[iS][iV][iP],iS,iV,iP,0,iBin);
			}
		}
	}
	
	// 1: sc_1sub
	for(unsigned int iS=5; iS<6; iS++)
	{
		for(unsigned int iV=2; iV<4; iV++)
		{
			for(unsigned int iP=0; iP<1; iP++)
			{
				draw_sepa(sc_1sub[1][iS][iV][iP],sc_1sub[0][iS][iV][iP],ratio_sc_1sub[iS][iV][iP],iS,iV,iP,1,iBin);
			}
		}
	}
	
	// 2: ac_1sub
	for(unsigned int iS=5; iS<6; iS++)
	{
		for(unsigned int iV=2; iV<3; iV++)
		{
			for(unsigned int iP=0; iP<1; iP++)
			{
				draw_sepa(ac_1sub[1][iS][iV][iP],ac_1sub[0][iS][iV][iP],ratio_ac_1sub[iS][iV][iP],iS,iV,iP,2,iBin);
			}
		}
	}

	//-----------------

	vector<TGraphErrors*> gVec;
	
	// 0: c2_1sub
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c2_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_c2_1sub[iV][iP]);
			gVec.push_back(sysLw_c2_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,0,iBin);
			gVec.clear();
		}
	}
	// 1: c4_1sub
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c4_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_c4_1sub[iV][iP]);
			gVec.push_back(sysLw_c4_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,1,iBin);
			gVec.clear();
		}
	}
	// 2: c6_1sub
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c6_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_c6_1sub[iV][iP]);
			gVec.push_back(sysLw_c6_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,2,iBin);
			gVec.clear();
		}
	}
	// 3: nc4_1sub
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nc4_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_nc4_1sub[iV][iP]);
			gVec.push_back(sysLw_nc4_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,3,iBin);
			gVec.clear();
		}
	}
	// 4: nc6_1sub
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nc6_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_nc6_1sub[iV][iP]);
			gVec.push_back(sysLw_nc6_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,4,iBin);
			gVec.clear();
		}
	}
	// 5: sc_1sub
	for(unsigned int iV=2; iV<4; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_sc_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_sc_1sub[iV][iP]);
			gVec.push_back(sysLw_sc_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,5,iBin);
			gVec.clear();
		}
	}
	// 6: nsc_1sub
	for(unsigned int iV=2; iV<4; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nsc_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_nsc_1sub[iV][iP]);
			gVec.push_back(sysLw_nsc_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,6,iBin);
			gVec.clear();
		}
	}
	// 7: ac_1sub
	for(unsigned int iV=2; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_ac_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_ac_1sub[iV][iP]);
			gVec.push_back(sysLw_ac_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,7,iBin);
			gVec.clear();
		}
	}
	// 8: nac_1sub
	for(unsigned int iV=2; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nac_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_nac_1sub[iV][iP]);
			gVec.push_back(sysLw_nac_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,8,iBin);
			gVec.clear();
		}
	}
	// 9: isGauss_1sub
	for(unsigned int iV=2; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_isGauss_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_isGauss_1sub[iV][iP]);
			gVec.push_back(sysLw_isGauss_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,9,iBin);
			gVec.clear();
		}
	}
	// 10: isPower_1sub
	for(unsigned int iV=2; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_isPower_1sub[iS][iV][iP]);
			gVec.push_back(sysUp_isPower_1sub[iV][iP]);
			gVec.push_back(sysLw_isPower_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,10,iBin);
			gVec.clear();
		}
	}
	// 11: c2_3sub
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c2_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_c2_3sub[iV][iP]);
			gVec.push_back(sysLw_c2_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,11,iBin);
			gVec.clear();
		}
	}
	// 12: c4_3sub
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_c4_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_c4_3sub[iV][iP]);
			gVec.push_back(sysLw_c4_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,12,iBin);
			gVec.clear();
		}
	}
	// 13: nc4_3sub
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nc4_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_nc4_3sub[iV][iP]);
			gVec.push_back(sysLw_nc4_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,13,iBin);
			gVec.clear();
		}
	}
	// 14: sc_3sub
	for(unsigned int iV=2; iV<4; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_sc_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_sc_3sub[iV][iP]);
			gVec.push_back(sysLw_sc_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,14,iBin);
			gVec.clear();
		}
	}
	// 15: nsc_3sub
	for(unsigned int iV=2; iV<4; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nsc_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_nsc_3sub[iV][iP]);
			gVec.push_back(sysLw_nsc_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,15,iBin);
			gVec.clear();
		}
	}
	// 16: ac_3sub
	for(unsigned int iV=2; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_ac_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_ac_3sub[iV][iP]);
			gVec.push_back(sysLw_ac_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,16,iBin);
			gVec.clear();
		}
	}
	// 17: nac_3sub
	for(unsigned int iV=2; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			for(unsigned int iS=0; iS<NS; iS++) gVec.push_back(ratio_nac_3sub[iS][iV][iP]);
			gVec.push_back(sysUp_nac_3sub[iV][iP]);
			gVec.push_back(sysLw_nac_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,17,iBin);
			gVec.clear();
		}
	}


	fOut->Close();
}

void Plot::initialize(unsigned int iBin)
{
	cout<<"initialize..."<<endl;

	sprintf(name,"../OUTPUT/Sys_bin%d.root",iBin);
	TFile* fIn = new TFile(name,"READ");
	for(unsigned int iF=0; iF<2; iF++)
	{
		for(unsigned int iS=0; iS<NS; iS++)
		{
			for(unsigned int iV=0; iV<NV; iV++)
			{
				for(unsigned int iP=0; iP<NP; iP++)
				{
					readHist_FSVP(fIn,c4_1sub[iF][iS][iV][iP],"c4_1sub",iF,iS,iV,iP);
					readHist_FSVP(fIn,sc_1sub[iF][iS][iV][iP],"sc_3sub",iF,iS,iV,iP);
					readHist_FSVP(fIn,ac_1sub[iF][iS][iV][iP],"ac_3sub",iF,iS,iV,iP);
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
				readHist_SVP(fIn,ratio_c2_1sub[iS][iV][iP],"ratio_c2_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_c4_1sub[iS][iV][iP],"ratio_c4_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_c6_1sub[iS][iV][iP],"ratio_c6_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nc4_1sub[iS][iV][iP],"ratio_nc4_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nc6_1sub[iS][iV][iP],"ratio_nc6_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_sc_1sub[iS][iV][iP],"ratio_sc_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nsc_1sub[iS][iV][iP],"ratio_nsc_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_ac_1sub[iS][iV][iP],"ratio_ac_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nac_1sub[iS][iV][iP],"ratio_nac_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_isGauss_1sub[iS][iV][iP],"ratio_isGauss_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_isPower_1sub[iS][iV][iP],"ratio_isPower_1sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_c2_3sub[iS][iV][iP],"ratio_c2_3sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_c4_3sub[iS][iV][iP],"ratio_c4_3sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nc4_3sub[iS][iV][iP],"ratio_nc4_3sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_sc_3sub[iS][iV][iP],"ratio_sc_3sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nsc_3sub[iS][iV][iP],"ratio_nsc_3sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_ac_3sub[iS][iV][iP],"ratio_ac_3sub",iS,iV,iP);
				readHist_SVP(fIn,ratio_nac_3sub[iS][iV][iP],"ratio_nac_3sub",iS,iV,iP);
			}
		}
	}
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			readHist_VP(fIn,sysUp_c2_1sub[iV][iP],"sysUp_c2_1sub",iV,iP);
			readHist_VP(fIn,sysUp_c4_1sub[iV][iP],"sysUp_c4_1sub",iV,iP);
			readHist_VP(fIn,sysUp_c6_1sub[iV][iP],"sysUp_c6_1sub",iV,iP);
			readHist_VP(fIn,sysUp_nc4_1sub[iV][iP],"sysUp_nc4_1sub",iV,iP);
			readHist_VP(fIn,sysUp_nc6_1sub[iV][iP],"sysUp_nc6_1sub",iV,iP);
			readHist_VP(fIn,sysUp_sc_1sub[iV][iP],"sysUp_sc_1sub",iV,iP);
			readHist_VP(fIn,sysUp_nsc_1sub[iV][iP],"sysUp_nsc_1sub",iV,iP);
			readHist_VP(fIn,sysUp_ac_1sub[iV][iP],"sysUp_ac_1sub",iV,iP);
			readHist_VP(fIn,sysUp_nac_1sub[iV][iP],"sysUp_nac_1sub",iV,iP);
			readHist_VP(fIn,sysUp_isGauss_1sub[iV][iP],"sysUp_isGauss_1sub",iV,iP);
			readHist_VP(fIn,sysUp_isPower_1sub[iV][iP],"sysUp_isPower_1sub",iV,iP);
			readHist_VP(fIn,sysUp_c2_3sub[iV][iP],"sysUp_c2_3sub",iV,iP);
			readHist_VP(fIn,sysUp_c4_3sub[iV][iP],"sysUp_c4_3sub",iV,iP);
			readHist_VP(fIn,sysUp_nc4_3sub[iV][iP],"sysUp_nc4_3sub",iV,iP);
			readHist_VP(fIn,sysUp_sc_3sub[iV][iP],"sysUp_sc_3sub",iV,iP);
			readHist_VP(fIn,sysUp_nsc_3sub[iV][iP],"sysUp_nsc_3sub",iV,iP);
			readHist_VP(fIn,sysUp_ac_3sub[iV][iP],"sysUp_ac_3sub",iV,iP);
			readHist_VP(fIn,sysUp_nac_3sub[iV][iP],"sysUp_nac_3sub",iV,iP);

			readHist_VP(fIn,sysLw_c2_1sub[iV][iP],"sysLw_c2_1sub",iV,iP);
			readHist_VP(fIn,sysLw_c4_1sub[iV][iP],"sysLw_c4_1sub",iV,iP);
			readHist_VP(fIn,sysLw_c6_1sub[iV][iP],"sysLw_c6_1sub",iV,iP);
			readHist_VP(fIn,sysLw_nc4_1sub[iV][iP],"sysLw_nc4_1sub",iV,iP);
			readHist_VP(fIn,sysLw_nc6_1sub[iV][iP],"sysLw_nc6_1sub",iV,iP);
			readHist_VP(fIn,sysLw_sc_1sub[iV][iP],"sysLw_sc_1sub",iV,iP);
			readHist_VP(fIn,sysLw_nsc_1sub[iV][iP],"sysLw_nsc_1sub",iV,iP);
			readHist_VP(fIn,sysLw_ac_1sub[iV][iP],"sysLw_ac_1sub",iV,iP);
			readHist_VP(fIn,sysLw_nac_1sub[iV][iP],"sysLw_nac_1sub",iV,iP);
			readHist_VP(fIn,sysLw_isGauss_1sub[iV][iP],"sysLw_isGauss_1sub",iV,iP);
			readHist_VP(fIn,sysLw_isPower_1sub[iV][iP],"sysLw_isPower_1sub",iV,iP);
			readHist_VP(fIn,sysLw_c2_3sub[iV][iP],"sysLw_c2_3sub",iV,iP);
			readHist_VP(fIn,sysLw_c4_3sub[iV][iP],"sysLw_c4_3sub",iV,iP);
			readHist_VP(fIn,sysLw_nc4_3sub[iV][iP],"sysLw_nc4_3sub",iV,iP);
			readHist_VP(fIn,sysLw_sc_3sub[iV][iP],"sysLw_sc_3sub",iV,iP);
			readHist_VP(fIn,sysLw_nsc_3sub[iV][iP],"sysLw_nsc_3sub",iV,iP);
			readHist_VP(fIn,sysLw_ac_3sub[iV][iP],"sysLw_ac_3sub",iV,iP);
			readHist_VP(fIn,sysLw_nac_3sub[iV][iP],"sysLw_nac_3sub",iV,iP);
		}
	}
}

void Plot::finalize()
{
	cout<<"finalize..."<<endl;
}

void Plot::draw_graph(vector<TGraphErrors*> vIn, int iV, int iP, int iOpt, unsigned int iBin)
{
	int NG = vIn.size();
	TGraphErrors* gIn[10];
	for(int iG=0; iG<NG; iG++)
	{
		gIn[iG] = (TGraphErrors*)vIn.at(iG)->Clone("gIn");
		if(iG==0) styleGraph(gIn[iG],0);
		else styleGraph(gIn[iG],iG-1);
	}
	for(int iG=1; iG<NG; iG++)
	{
		for(int iB=0; iB<gIn[0]->GetN(); iB++) gIn[iG]->SetPointError(iB,0,0);
	}
	gIn[NG-2]->SetLineColor(kBlack);  gIn[NG-2]->SetLineWidth(4);
	gIn[NG-1]->SetLineColor(kBlack);  gIn[NG-1]->SetLineWidth(4);

	TLatex* tex = new TLatex();
	tex->SetTextFont(42);
	tex->SetTextSize(0.045);
	tex->SetTextAlign(12);
	tex->SetNDC(1);
	TLine* lin = new TLine();
	lin->SetLineColor(1);
	lin->SetLineStyle(2);
	lin->SetLineWidth(2);
	TLegend* leg = new TLegend(0.15,0.125,0.95,0.25);
	leg->SetTextSize(0.05);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(3);
	leg->AddEntry(gIn[0],"#font[42]{stat. err}","f");
	leg->AddEntry(gIn[1],"#font[42]{lower eff.}","p");
	leg->AddEntry(gIn[2],"#font[42]{higher eff.}","p");
	leg->AddEntry(gIn[3],"#font[42]{tight sel.}","p");
	leg->AddEntry(gIn[4],"#font[42]{pileup}","p");
	leg->AddEntry(gIn[5],"#font[42]{MC closure}","p");
	leg->AddEntry(gIn[6],"#font[42]{flattening}","p");
	leg->AddEntry(gIn[7],"#font[42]{Combined}","L");

	double xMin =  0;
	double xMax =  80;
	double yMin =  1;
	double yMax = -1;
	for(int iG=0; iG<NG; iG++) getYrange(gIn[iG],yMin,yMax);
	double diff = yMax-yMin;
	yMax += 0.5*diff;
	yMin -= 0.5*diff;

	TH1D* hAxis = new TH1D("hAxis","",1000,xMin,xMax);
	for(int i=0; i<1000; i++) hAxis->SetBinContent(i+1,1E9);
	styleGraph(hAxis,0);

	TCanvas* cOut = new TCanvas("cOut","",400,400);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	cOut->cd();
	gPad->SetTicks(1,1);
	gPad->SetTopMargin(0.075);
	gPad->SetBottomMargin(0.1);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.025);
	hAxis->GetXaxis()->SetTitle("Centrality / %");
	hAxis->GetXaxis()->SetTitleOffset(1.15);
	hAxis->GetXaxis()->SetRangeUser(xMin,xMax);
	if(iOpt==0) sprintf(name,"c_{%d}{2}",iV);
	if(iOpt==1) sprintf(name,"c_{%d}{4}",iV);
	if(iOpt==2) sprintf(name,"c_{%d}{6}",iV);
	if(iOpt==3) sprintf(name,"#hat{c}_{%d}{4}",iV);
	if(iOpt==4) sprintf(name,"#hat{c}_{%d}{6}",iV);
	if(iOpt==5)
	{
		if(iV==0) sprintf(name,"SC(1,2)");
		if(iV==1) sprintf(name,"SC(1,3)");
		if(iV==2) sprintf(name,"SC(2,3)");
		if(iV==3) sprintf(name,"SC(2,4)");
	}
	if(iOpt==6)
	{
		if(iV==0) sprintf(name,"#hat{SC}(1,2)");
		if(iV==1) sprintf(name,"#hat{SC}(1,3)");
		if(iV==2) sprintf(name,"#hat{SC}(2,3)");
		if(iV==3) sprintf(name,"#hat{SC}(2,4)");
	}
	if(iOpt==7)
	{
		if(iV==0) sprintf(name,"ASC(1,1,2)");
		if(iV==1) sprintf(name,"ASC(1,2,3)");
		if(iV==2) sprintf(name,"ASC(2,2,4)");
	}
	if(iOpt==8)
	{
		if(iV==0) sprintf(name,"#hat{ASC}(1,1,2)");
		if(iV==1) sprintf(name,"#hat{ASC}(1,2,3)");
		if(iV==2) sprintf(name,"#hat{ASC}(2,2,4)");
	}
	if(iOpt==9) sprintf(name,"v_{%d} Gaussian?",iV);
	if(iOpt==10) sprintf(name,"v_{%d} Power-law?",iV);
	if(iOpt==11) sprintf(name,"c_{%d}{2}",iV);
	if(iOpt==12) sprintf(name,"c_{%d}{4}",iV);
	if(iOpt==13) sprintf(name,"#hat{c}_{%d}{4}",iV);
	if(iOpt==14)
	{
		if(iV==0) sprintf(name,"SC(1,2)");
		if(iV==1) sprintf(name,"SC(1,3)");
		if(iV==2) sprintf(name,"SC(2,3)");
		if(iV==3) sprintf(name,"SC(2,4)");
	}
	if(iOpt==15)
	{
		if(iV==0) sprintf(name,"#hat{SC}(1,2)");
		if(iV==1) sprintf(name,"#hat{SC}(1,3)");
		if(iV==2) sprintf(name,"#hat{SC}(2,3)");
		if(iV==3) sprintf(name,"#hat{SC}(2,4)");
	}
	if(iOpt==16)
	{
		if(iV==0) sprintf(name,"ASC(1,1,2)");
		if(iV==1) sprintf(name,"ASC(1,2,3)");
		if(iV==2) sprintf(name,"ASC(2,2,4)");
	}
	if(iOpt==17)
	{
		if(iV==0) sprintf(name,"#hat{ASC}(1,1,2)");
		if(iV==1) sprintf(name,"#hat{ASC}(1,2,3)");
		if(iV==2) sprintf(name,"#hat{ASC}(2,2,4)");
	}
	hAxis->GetYaxis()->SetTitle(name);
	hAxis->GetYaxis()->SetTitleOffset(1.6);
	hAxis->GetYaxis()->SetRangeUser(yMin,yMax);
	hAxis->Draw();
	gIn[0]->SetFillColor(1);
	gIn[0]->SetFillStyle(3004);
	gIn[0]->Draw("E3");
	for(int iG=1; iG<NG-2; iG++) gIn[iG]->Draw("PL");
	gIn[NG-2]->Draw("L");
	gIn[NG-1]->Draw("L");
	tex->DrawLatex(0.175,0.875,"#font[72]{ATLAS} #font[62]{Internal}");
	tex->DrawLatex(0.175,0.82,"#font[42]{Pb+Pb #sqrt{s_{NN}}=5.02 TeV}");
	if(iOpt<=10) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV   standard}",minPt[iP],maxPt[iP]);
	else sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV   3-subevent}",minPt[iP],maxPt[iP]);
	tex->DrawLatex(0.175,0.76,name);
	lin->DrawLine(xMin,0,xMax,0);
	leg->Draw();

	if(iOpt==0) sprintf(name,"../PLOT/bin%d/sys_c2_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==1) sprintf(name,"../PLOT/bin%d/sys_c4_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==2) sprintf(name,"../PLOT/bin%d/sys_c6_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==3) sprintf(name,"../PLOT/bin%d/sys_nc4_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==4) sprintf(name,"../PLOT/bin%d/sys_nc6_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==5) sprintf(name,"../PLOT/bin%d/sys_sc_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==6) sprintf(name,"../PLOT/bin%d/sys_nsc_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==7) sprintf(name,"../PLOT/bin%d/sys_ac_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==8) sprintf(name,"../PLOT/bin%d/sys_nac_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==9) sprintf(name,"../PLOT/bin%d/sys_isGauss_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==10) sprintf(name,"../PLOT/bin%d/sys_isPower_1sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==11) sprintf(name,"../PLOT/bin%d/sys_c2_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==12) sprintf(name,"../PLOT/bin%d/sys_c4_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==13) sprintf(name,"../PLOT/bin%d/sys_nc4_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==14) sprintf(name,"../PLOT/bin%d/sys_sc_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==15) sprintf(name,"../PLOT/bin%d/sys_nsc_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==16) sprintf(name,"../PLOT/bin%d/sys_ac_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==17) sprintf(name,"../PLOT/bin%d/sys_nac_3sub_Har%d_Pt%d.pdf",iBin,iV,iP);
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"sys_c2_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==1) sprintf(name,"sys_c4_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==2) sprintf(name,"sys_c6_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==3) sprintf(name,"sys_nc4_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==4) sprintf(name,"sys_nc6_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==5) sprintf(name,"sys_sc_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==6) sprintf(name,"sys_nsc_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==7) sprintf(name,"sys_ac_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==8) sprintf(name,"sys_nac_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==9) sprintf(name,"sys_isGauss_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==10) sprintf(name,"sys_isPower_1sub_Har%d_Pt%d",iV,iP);
	if(iOpt==11) sprintf(name,"sys_c2_3sub_Har%d_Pt%d",iV,iP);
	if(iOpt==12) sprintf(name,"sys_c4_3sub_Har%d_Pt%d",iV,iP);
	if(iOpt==13) sprintf(name,"sys_nc4_3sub_Har%d_Pt%d",iV,iP);
	if(iOpt==14) sprintf(name,"sys_sc_3sub_Har%d_Pt%d",iV,iP);
	if(iOpt==15) sprintf(name,"sys_nsc_3sub_Har%d_Pt%d",iV,iP);
	if(iOpt==16) sprintf(name,"sys_ac_3sub_Har%d_Pt%d",iV,iP);
	if(iOpt==17) sprintf(name,"sys_nac_3sub_Har%d_Pt%d",iV,iP);
	cOut->SetName(name);
	cOut->Write();

	delete hAxis;
	delete cOut;
}

void Plot::draw_sepa(TGraphErrors* g0, TGraphErrors* g1, TGraphErrors* gR, int iS, int iV, int iP, int iOpt, unsigned int iBin)
{
	styleGraph(g0,0); styleGraph(g1,1); styleGraph(gR,0);
	for(int iB=0; iB<g1->GetN(); iB++)
	{
		double x; double y;
		g1->GetPoint(iB,x,y);
		x += 2;
		g1->SetPoint(iB,x,y);
	}

	TLatex* tex = new TLatex();
	tex->SetTextFont(42);
	tex->SetTextSize(0.045);
	tex->SetTextAlign(12);
	tex->SetNDC(1);
	TLine* lin = new TLine();
	lin->SetLineColor(1);
	lin->SetLineStyle(2);
	lin->SetLineWidth(2);
	TLegend* leg = new TLegend(0.225,0.025,0.95,0.125);
	leg->SetTextSize(0.05);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(2);
	if(iS==0)
	{
		leg->AddEntry(g0,"#font[42]{default}","p");
		leg->AddEntry(g1,"#font[42]{default}","p");
	}
	if(iS==1)
	{
		leg->AddEntry(g0,"#font[42]{default}","p");
		leg->AddEntry(g1,"#font[42]{low track eff.}","p");
	}
	if(iS==2)
	{
		leg->AddEntry(g0,"#font[42]{default}","p");
		leg->AddEntry(g1,"#font[42]{high track eff.}","p");
	}
	if(iS==3)
	{
		leg->AddEntry(g0,"#font[42]{default}","p");
		leg->AddEntry(g1,"#font[42]{tight track sel.}","p");
	}
	if(iS==4)
	{
		leg->AddEntry(g0,"#font[42]{default}","p");
		leg->AddEntry(g1,"#font[42]{loose pileup rej.}","p");
	}
	if(iS==5)
	{
		leg->AddEntry(g0,"#font[42]{recon}","p");
		leg->AddEntry(g1,"#font[42]{truth}","p");
	}
	if(iS==6)
	{
		leg->AddEntry(g0,"#font[42]{default}","p");
		leg->AddEntry(g1,"#font[42]{no flat.}","p");
	}

	double xMin =  0;
	double xMax =  70;
	double yMin =  1;
	double yMax = -1;
	getYrange(g0,yMin,yMax);
	getYrange(g1,yMin,yMax);
	double diff = yMax-yMin;
	yMax += 0.5*diff;
	yMin -= 0.5*diff;

	TH1D* hAxis = new TH1D("hAxis","",1000,xMin,xMax);
	for(int i=0; i<1000; i++) hAxis->SetBinContent(i+1,1E9);
	styleGraph(hAxis,0);

	TCanvas* cOut = new TCanvas("cOut","",300,500);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	cOut->cd();
	TPad* pad1 = new TPad("pad1","",0,0.3,1,1);
	pad1->Draw();
	pad1->cd();
	pad1->SetTicks(1,1);
	pad1->SetTopMargin(0.075);
	pad1->SetBottomMargin(0.01);
	pad1->SetLeftMargin(0.2);
	pad1->SetRightMargin(0.025);
	hAxis->GetXaxis()->SetTitle("Centrality / %");
	hAxis->GetXaxis()->SetTitleOffset(5);
	hAxis->GetXaxis()->SetRangeUser(xMin,xMax);
	if(iOpt==0) sprintf(name,"c_{%d}{4}",iV);
	if(iOpt==1)
	{
		if(iV==2) sprintf(name,"sc_{2,3}{4}");
		if(iV==3) sprintf(name,"sc_{2,4}{4}");
	}
	if(iOpt==2)
	{
		if(iV==2) sprintf(name,"ac_{2,4}{3}");
	}
	hAxis->GetYaxis()->SetTitle(name);
	hAxis->GetYaxis()->SetTitleOffset(3.3);
	hAxis->GetYaxis()->SetRangeUser(yMin,yMax);
	hAxis->DrawCopy();
	g0->Draw("P");
	g1->Draw("P");
	tex->DrawLatex(0.25,0.875,"#font[72]{ATLAS} #font[62]{Internal}");
	tex->DrawLatex(0.25,0.82,"#font[42]{Pb+Pb #sqrt{s_{NN}}=5.02 TeV}");
	if(iOpt<=10) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV   standard}",minPt[iP],maxPt[iP]);
	else sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV   3-subevent}",minPt[iP],maxPt[iP]);
	tex->DrawLatex(0.25,0.76,name);
	if(yMax>0 && yMin<0) lin->DrawLine(xMin,0,xMax,0);
	leg->Draw();
	cOut->cd();
	TPad* pad2 = new TPad("pad2","",0,0,1,0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTicks(1,1);
	pad2->SetTopMargin(0.01);
	pad2->SetBottomMargin(0.2);
	pad2->SetLeftMargin(0.2);
	pad2->SetRightMargin(0.025);
	hAxis->GetXaxis()->SetTitle("Centrality / %");
	hAxis->GetXaxis()->SetTitleOffset(2.7);
	hAxis->GetXaxis()->SetRangeUser(xMin,xMax);
	sprintf(name,"check / default - 1");
	hAxis->GetYaxis()->SetTitle(name);
	hAxis->GetYaxis()->SetTitleOffset(3.3);
	hAxis->GetYaxis()->SetRangeUser(-0.14,0.14);
	hAxis->DrawCopy();
	gR->Draw("P");
	lin->DrawLine(xMin,0,xMax,0);

	if(iOpt==0) sprintf(name,"../PLOT/bin%d/sys%d_c4_1sub_Har%d_Pt%d.pdf",iBin,iS,iV,iP);
	if(iOpt==1) sprintf(name,"../PLOT/bin%d/sys%d_sc_1sub_Har%d_Pt%d.pdf",iBin,iS,iV,iP);
	if(iOpt==2) sprintf(name,"../PLOT/bin%d/sys%d_ac_1sub_Har%d_Pt%d.pdf",iBin,iS,iV,iP);
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"sys%d_c4_1sub_Har%d_Pt%d",iS,iV,iP);
	if(iOpt==1) sprintf(name,"sys%d_sc_1sub_Har%d_Pt%d",iS,iV,iP);
	if(iOpt==2) sprintf(name,"sys%d_ac_1sub_Har%d_Pt%d",iS,iV,iP);
	cOut->SetName(name);
	cOut->Write();

	delete hAxis;
	delete cOut;
	//delete pad1;
	//delete pad2;
}

void Plot::styleGraph(TGraph* hIn, int k)
{
	hIn->SetTitle("");
	hIn->GetYaxis()->SetNdivisions(505);
	hIn->GetYaxis()->SetTitleOffset(1.);    hIn->GetXaxis()->SetTitleOffset(1.);
	hIn->GetYaxis()->SetTitleFont(43);      hIn->GetYaxis()->SetTitleSize(15);
	hIn->GetYaxis()->SetLabelFont(43);      hIn->GetYaxis()->SetLabelSize(15);
	hIn->GetXaxis()->SetTitleFont(43);      hIn->GetXaxis()->SetTitleSize(15);
	hIn->GetXaxis()->SetLabelFont(43);      hIn->GetXaxis()->SetLabelSize(15);
	hIn->SetMarkerStyle(mS[k]);
	hIn->SetMarkerColor(mC[k]);
	hIn->SetLineStyle(lS[k]);
	hIn->SetLineColor(lC[k]);
	hIn->SetMarkerSize(1.2);
	hIn->SetLineWidth(2);
}

void Plot::styleGraph(TH1D* hIn, int k)
{
	hIn->SetTitle("");
	hIn->GetYaxis()->SetTitleOffset(1.);    hIn->GetXaxis()->SetTitleOffset(1.);
	hIn->GetYaxis()->SetTitleFont(43);      hIn->GetYaxis()->SetTitleSize(15);
	hIn->GetYaxis()->SetLabelFont(43);      hIn->GetYaxis()->SetLabelSize(15);
	hIn->GetXaxis()->SetTitleFont(43);      hIn->GetXaxis()->SetTitleSize(15);
	hIn->GetXaxis()->SetLabelFont(43);      hIn->GetXaxis()->SetLabelSize(15);
	hIn->SetMarkerStyle(mS[k]);
	hIn->SetMarkerColor(mC[k]);
	hIn->SetLineStyle(lS[k]);
	hIn->SetLineColor(lC[k]);
	hIn->SetMarkerSize(1.2);
	hIn->SetLineWidth(2);
}

void Plot::getYrange(TGraph* hIn, double& yMin, double& yMax)
{
	double x; double y;
	for(int i=0; i<hIn->GetN(); i++)
	{
		hIn->GetPoint(i,x,y);
		if(x>80) continue;
		if(y<yMin) yMin = y;
		if(y>yMax) yMax = y;
	}
}

void Plot::readHist_FSVP(TFile* fIn, TGraphErrors*& hIn, const char* hName, int iF, int iS, int iV, int iP)
{
	sprintf(name,"%s_File%d_Sys%d_Har%d_Pt%d",hName,iF,iS,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
}

void Plot::readHist_SVP(TFile* fIn, TGraphErrors*& hIn, const char* hName, int iS, int iV, int iP)
{
	sprintf(name,"%s_Sys%d_Har%d_Pt%d",hName,iS,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
}

void Plot::readHist_VP(TFile* fIn, TGraphErrors*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
}
