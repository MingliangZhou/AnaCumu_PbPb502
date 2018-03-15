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
			delete c2_1sub[iV][iP];
			delete c4_1sub[iV][iP];
			delete c6_1sub[iV][iP];
			delete nc4_1sub[iV][iP];
			delete nc6_1sub[iV][iP];
			delete sc_1sub[iV][iP];
			delete nsc_1sub[iV][iP];
			delete ac_1sub[iV][iP];
			delete nac_1sub[iV][iP];
			delete isGauss_1sub[iV][iP];
			delete isPower_1sub[iV][iP];
			delete c2_3sub[iV][iP];
			delete c4_3sub[iV][iP];
			delete nc4_3sub[iV][iP];
			delete sc_3sub[iV][iP];
			delete nsc_3sub[iV][iP];
			delete ac_3sub[iV][iP];
			delete nac_3sub[iV][iP];
		}
	}
}

void Plot::execute(unsigned int iBin)
{
  cout<<"execute..."<<endl;

	sprintf(name,"../PLOT/bin%d/canvas.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();
	vector<TGraphErrors*> gVec;

	// 0: mtd_c2
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(c2_1sub[iV][iP]);
			gVec.push_back(c2_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,0,iBin);
			gVec.clear();
		}
	}

	// 1: mtd_c4
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(c4_1sub[iV][iP]);
			gVec.push_back(c4_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,1,iBin);
			gVec.clear();
		}
	}

	// 2: mtd_c6
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(c6_1sub[iV][iP]);
			gVec.push_back(c6_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,2,iBin);
			gVec.clear();
		}
	}

	// 3: mtd_nc4
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(nc4_1sub[iV][iP]);
			gVec.push_back(nc4_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,3,iBin);
			gVec.clear();
		}
	}

	// 4: mtd_nc6
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(nc6_1sub[iV][iP]);
			gVec.push_back(nc6_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,4,iBin);
			gVec.clear();
		}
	}

	// 5: mtd_sc
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(sc_1sub[iV][iP]);
			gVec.push_back(sc_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,5,iBin);
			gVec.clear();
		}
	}

	// 6: mtd_nsc
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(nsc_1sub[iV][iP]);
			gVec.push_back(nsc_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,6,iBin);
			gVec.clear();
		}
	}

	// 7: mtd_ac
	for(unsigned int iV=0; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(ac_1sub[iV][iP]);
			gVec.push_back(ac_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,7,iBin);
			gVec.clear();
		}
	}

	// 8: mtd_nac
	for(unsigned int iV=0; iV<3; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(nac_1sub[iV][iP]);
			gVec.push_back(nac_3sub[iV][iP]);
			draw_graph(gVec,iV,iP,8,iBin);
			gVec.clear();
		}
	}

	// 9: mtd_isGauss
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(isGauss_1sub[iV][iP]);
			gVec.push_back(isGauss_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,9,iBin);
			gVec.clear();
		}
	}

	// 10: mtd_isPower
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			gVec.push_back(isPower_1sub[iV][iP]);
			gVec.push_back(isPower_1sub[iV][iP]);
			draw_graph(gVec,iV,iP,10,iBin);
			gVec.clear();
		}
	}

	fOut->Close();
}

void Plot::initialize(unsigned int iBin)
{
  cout<<"initialize..."<<endl;
	
	sprintf(name,"../OUTPUT/Phase3/Phase3_bin%d.root",iBin);
	TFile* fIn = new TFile(name,"READ");
	for(unsigned int iV=0; iV<NV; iV++)
	{
		for(unsigned int iP=0; iP<NP; iP++)
		{
			readHist_VP(fIn,c2_1sub[iV][iP],"c2_1sub",iV,iP);
			readHist_VP(fIn,c4_1sub[iV][iP],"c4_1sub",iV,iP);
			readHist_VP(fIn,c6_1sub[iV][iP],"c6_1sub",iV,iP);
			readHist_VP(fIn,nc4_1sub[iV][iP],"nc4_1sub",iV,iP);
			readHist_VP(fIn,nc6_1sub[iV][iP],"nc6_1sub",iV,iP);
			readHist_VP(fIn,sc_1sub[iV][iP],"sc_1sub",iV,iP);
			readHist_VP(fIn,nsc_1sub[iV][iP],"nsc_1sub",iV,iP);
			readHist_VP(fIn,ac_1sub[iV][iP],"ac_1sub",iV,iP);
			readHist_VP(fIn,nac_1sub[iV][iP],"nac_1sub",iV,iP);
			readHist_VP(fIn,isGauss_1sub[iV][iP],"isGauss_1sub",iV,iP);
			readHist_VP(fIn,isPower_1sub[iV][iP],"isPower_1sub",iV,iP);
			readHist_VP(fIn,c2_3sub[iV][iP],"c2_3sub",iV,iP);
			readHist_VP(fIn,c4_3sub[iV][iP],"c4_3sub",iV,iP);
			readHist_VP(fIn,nc4_3sub[iV][iP],"nc4_3sub",iV,iP);
			readHist_VP(fIn,sc_3sub[iV][iP],"sc_3sub",iV,iP);
			readHist_VP(fIn,nsc_3sub[iV][iP],"nsc_3sub",iV,iP);
			readHist_VP(fIn,ac_3sub[iV][iP],"ac_3sub",iV,iP);
			readHist_VP(fIn,nac_3sub[iV][iP],"nac_3sub",iV,iP);
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
		styleGraph(gIn[iG],iG);
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
	TLegend* leg = new TLegend(0.15,0.125,0.95,0.25);
	leg->SetTextSize(0.05);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(3);
	leg->AddEntry(gIn[0],"#font[42]{standard}","p");
	leg->AddEntry(gIn[1],"#font[42]{3-subevent}","p");

	double xMin =  0;
	double xMax = 80;
	double yMin =  1;
	double yMax = -1;
	for(int iG=0; iG<NG; iG++) getYrange(gIn[iG],yMin,yMax);
	double diff = yMax-yMin;
	yMax += 0.5*diff;
	yMin -= 0.5*diff;

	TH1D* hAxis = new TH1D("hAxis","",800,xMin,xMax);
	hAxis->Fill(1,1E9);
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
	if(iOpt==9) sprintf(name,"v_{%d} Gauss check",iV);
	if(iOpt==10) sprintf(name,"v_{%d} Power check",iV);
	hAxis->GetYaxis()->SetTitle(name);
	hAxis->GetYaxis()->SetTitleOffset(1.6);
	hAxis->GetYaxis()->SetRangeUser(yMin,yMax);
	if(iOpt==9 || iOpt==10) hAxis->GetYaxis()->SetRangeUser(-0.5,2.5);
	hAxis->Draw();
	for(int iG=0; iG<NG; iG++) gIn[iG]->Draw("PL");
	tex->DrawLatex(0.175,0.875,"#font[72]{ATLAS} #font[62]{Internal}");
	tex->DrawLatex(0.175,0.82,"#font[42]{Xe+Xe #sqrt{s_{NN}}=5.44 TeV}");
	sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[iP],maxPt[iP]);
	tex->DrawLatex(0.175,0.76,name);
	if(iOpt!=9 && iOpt!=10) lin->DrawLine(xMin,0,xMax,0);
	else lin->DrawLine(xMin,1,xMax,1);
	leg->Draw();

	if(iOpt==0) sprintf(name,"../PLOT/bin%d/mtd_c2_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==1) sprintf(name,"../PLOT/bin%d/mtd_c4_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==2) sprintf(name,"../PLOT/bin%d/mtd_c6_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==3) sprintf(name,"../PLOT/bin%d/mtd_nc4_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==4) sprintf(name,"../PLOT/bin%d/mtd_nc6_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==5) sprintf(name,"../PLOT/bin%d/mtd_sc_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==6) sprintf(name,"../PLOT/bin%d/mtd_nsc_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==7) sprintf(name,"../PLOT/bin%d/mtd_ac_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==8) sprintf(name,"../PLOT/bin%d/mtd_nac_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==9) sprintf(name,"../PLOT/bin%d/mtd_isGauss_Har%d_Pt%d.pdf",iBin,iV,iP);
	if(iOpt==10) sprintf(name,"../PLOT/bin%d/mtd_isPower_Har%d_Pt%d.pdf",iBin,iV,iP);
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"mtd_c2_Har%d_Pt%d",iV,iP);
	if(iOpt==1) sprintf(name,"mtd_c4_Har%d_Pt%d",iV,iP);
	if(iOpt==2) sprintf(name,"mtd_c6_Har%d_Pt%d",iV,iP);
	if(iOpt==3) sprintf(name,"mtd_nc4_Har%d_Pt%d",iV,iP);
	if(iOpt==4) sprintf(name,"mtd_nc6_Har%d_Pt%d",iV,iP);
	if(iOpt==5) sprintf(name,"mtd_sc_Har%d_Pt%d",iV,iP);
	if(iOpt==6) sprintf(name,"mtd_nsc_Har%d_Pt%d",iV,iP);
	if(iOpt==7) sprintf(name,"mtd_ac_Har%d_Pt%d",iV,iP);
	if(iOpt==8) sprintf(name,"mtd_nac_Har%d_Pt%d",iV,iP);
	if(iOpt==9) sprintf(name,"mtd_isGauss_Har%d_Pt%d",iV,iP);
	if(iOpt==10) sprintf(name,"mtd_isPower_Har%d_Pt%d",iV,iP);
	cOut->SetName(name);
	cOut->Write();

	delete hAxis;
	delete cOut;
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
	hIn->SetLineWidth(2);
}

void Plot::getYrange(TGraph* hIn, double& yMin, double& yMax)
{
	double x; double y;
	for(int i=0; i<hIn->GetN(); i++)
	{
		hIn->GetPoint(i,x,y);
		if(x>50) continue;
		if(y<yMin) yMin = y;
		if(y>yMax) yMax = y;
	}
}

void Plot::readHist_VP(TFile* fIn, TGraphErrors*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
}
