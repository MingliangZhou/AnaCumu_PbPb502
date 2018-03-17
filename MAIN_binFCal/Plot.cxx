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

	sprintf(name,"/phenix/plhf/mzhou/AnaCumu_PbPb502/MAIN_binCent/PLOT/bin%d/canvas.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();
	vector<TGraphErrors*> gVec;

	fOut->Close();
}

void Plot::initialize(unsigned int iBin)
{
  cout<<"initialize..."<<endl;
	
	sprintf(name,"/phenix/plhf/mzhou/AnaCumu_PbPb502/MAIN_binCent/OUTPUT_Phase3/Phase3_bin%d.root",iBin);
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

void Plot::draw_graph(vector<TGraphErrors*> vIn, int iV, int iP, int iOpt)
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
	if(iOpt==0)
	{
		leg->AddEntry(gIn[0],"#font[42]{standard}","p");
		leg->AddEntry(gIn[1],"#font[42]{3-subevent}","p");
	}

	double xMin =  0;
	double xMax =  0;
	double yMin =  1;
	double yMax = -1;
	for(int iG=0; iG<NG; iG++) getYrange(gIn[iG],yMin,yMax,iOpt);
	double diff = yMax-yMin;
	yMax += 0.5*diff;
	yMin -= 0.5*diff;

	TH1D* hAxis = new TH1D("hAxis","",1,xMin,xMax);
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
	if(iOpt!=8 && iOpt!=19 && iOpt!=20 && iOpt!=27) sprintf(name,"Centrality / %%");
	else if(iOpt!=27) sprintf(name,"p_{T} / GeV");
	else sprintf(name,"v_{%d}{4}/v_{%d}{2}",iV,iV);
	hAxis->GetXaxis()->SetTitle(name);
	hAxis->GetXaxis()->SetTitleOffset(1.15);
	hAxis->GetXaxis()->SetRangeUser(xMin,xMax);
	if(iOpt==0) sprintf(name,"c_{%d}{2}",iV);
	hAxis->GetYaxis()->SetTitle(name);
	hAxis->GetYaxis()->SetTitleOffset(1.6);
	hAxis->GetYaxis()->SetRangeUser(yMin,yMax);
	if(iOpt==7 && iV!=2) hAxis->GetYaxis()->SetRangeUser(-0.1,0.2);
	if(iOpt==23 || iOpt==24) hAxis->GetYaxis()->SetRangeUser(-0.5,2.5);
	if(iOpt==27) hAxis->GetYaxis()->SetRangeUser(0.6,1.4);
	hAxis->Draw();
	for(int iG=0; iG<NG; iG++) gIn[iG]->Draw("PL");
	tex->DrawLatex(0.175,0.875,"#font[72]{ATLAS} #font[62]{Internal}");
	tex->DrawLatex(0.175,0.82,"#font[42]{Xe+Xe #sqrt{s_{NN}}=5.44 TeV}");
	if(iOpt== 0) sprintf(name,"#font[42]{%.1f<p_{T}^{RFP}<%.1f GeV}",minPtRef[iR],maxPtRef[iR]);
	tex->DrawLatex(0.175,0.76,name);
	if(iOpt<2 || iOpt==9 || iOpt==10 || iOpt==11 || iOpt==12 || iOpt==16 || iOpt==21 || iOpt==22) lin->DrawLine(xMin,0,xMax,0);
	if(iOpt==23 || iOpt==24) lin->DrawLine(xMin,1,xMax,1);
	leg->Draw();

	if(iOpt==0) sprintf(name,"PLOT/mtd_c2_Har%d_PtRef%d.pdf",iV,iR);
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"mtd_c2_Har%d_PtRef%d",iV,iR);
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

void Plot::getYrange(TGraph* hIn, double& yMin, double& yMax, int iOpt)
{
	double x; double y;
	for(int i=0; i<hIn->GetN(); i++)
	{
		hIn->GetPoint(i,x,y);
		if(iOpt!=8 && x>50) continue;
		if(iOpt==8 && x>10) continue;
		if(y<yMin) yMin = y;
		if(y>yMax) yMax = y;
	}
}

void Plot::readHist_VP(TFile* fIn, TGraphErrors*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
}
