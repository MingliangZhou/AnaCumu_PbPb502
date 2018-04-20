#include "Plot.h"

Plot::Plot()
{
	initialize();
	execute();
	finalize();
}

Plot::~Plot()
{
}

void Plot::execute()
{
  cout<<"execute..."<<endl;

	sprintf(name,"PLOT/canvas_LOOSE.root");
	TFile* fOut = new TFile(name,"RECREATE");
	fOut->cd();
	vector<TH1D*> gVec;

	for(unsigned int iP=0; iP<nPt; iP++)
	{
		gVec.push_back(hEtaEff[0][iP]);
		gVec.push_back(hEtaEff[3][iP]);
		gVec.push_back(hEtaEff[6][iP]);
		draw_graph(gVec,iP,0);
		gVec.clear();
	}

	for(unsigned int iP=0; iP<nPt; iP++)
	{
		gVec.push_back(hEtaFak[0][iP]);
		gVec.push_back(hEtaFak[3][iP]);
		gVec.push_back(hEtaFak[6][iP]);
		draw_graph(gVec,iP,1);
		gVec.clear();
	}

	fOut->Close();
}

void Plot::initialize()
{
  cout<<"initialize..."<<endl;
	
	sprintf(name,"trkEff_mon_LOOSE.root");
	TFile* fIn = new TFile(name,"READ");
	for(unsigned int iC=0; iC<nCent; iC++)
	{
		for(unsigned int iP=0; iP<nPt; iP++)
		{
			readHist_CP(fIn,hEtaEff[iC][iP],"hEtaEff",iC,iP);
			readHist_CP(fIn,hEtaFak[iC][iP],"hEtaFak",iC,iP);
		}
	}
}

void Plot::finalize()
{
	cout<<"finalize..."<<endl;
}

void Plot::draw_graph(vector<TH1D*> vIn, int iP, int iOpt)
{
	int NG = vIn.size();
	TH1D* gIn[10];
	for(int iG=0; iG<NG; iG++)
	{
		gIn[iG] = (TH1D*)vIn.at(iG)->Clone("gIn");
		styleGraph(gIn[iG],iG);
	}

	TLatex* tex = new TLatex();
	tex->SetTextSize(0.04);
	tex->SetTextAlign(12);
	tex->SetNDC(1);
	TLine* lin = new TLine();
	lin->SetLineColor(1);
	lin->SetLineStyle(2);
	lin->SetLineWidth(2);
	TLegend* leg = new TLegend(0.75,0.7,0.95,0.9);
	leg->SetTextSize(0.04);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(1);
	if(iOpt==0)
	{
		leg->AddEntry(gIn[0],"0-5%","p");
		leg->AddEntry(gIn[1],"20-30%","p");
		leg->AddEntry(gIn[2],"60-80%","p");
	}
	if(iOpt==1)
	{
		leg->AddEntry(gIn[0],"0-5%","p");
		leg->AddEntry(gIn[1],"20-30%","p");
		leg->AddEntry(gIn[2],"30-40%","p");
	}

	double xMin =  -2.5;
	double xMax =   2.5;

	TCanvas* cOut = new TCanvas("cOut","",400,400);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	cOut->cd();
	if(iOpt==0) gPad->SetLogy(0);
	if(iOpt==1) gPad->SetLogy(0);
	gPad->SetTicks(1,1);
	gPad->SetTopMargin(0.075);
	gPad->SetBottomMargin(0.1);
	gPad->SetLeftMargin(0.125);
	gPad->SetRightMargin(0.025);
	gIn[0]->GetXaxis()->SetTitle("#eta");
	gIn[0]->GetXaxis()->SetTitleOffset(1.1);
	gIn[0]->GetXaxis()->SetRangeUser(xMin,xMax);
	if(iOpt==0) gIn[0]->GetYaxis()->SetTitle("efficiency");
	if(iOpt==1) gIn[0]->GetYaxis()->SetTitle("fake rate");
	gIn[0]->GetYaxis()->SetTitleOffset(1.5);
	if(iOpt==0) gIn[0]->GetYaxis()->SetRangeUser(0.2,1);
	if(iOpt==1) gIn[0]->GetYaxis()->SetRangeUser(0.,0.3);
	gIn[0]->Draw("P");
	for(int iG=1; iG<NG; iG++) gIn[iG]->Draw("P same");
	tex->DrawLatex(0.175,0.875,"#font[72]{ATLAS} #font[62]{Simulation Internal}");
	tex->DrawLatex(0.175,0.825,"#font[42]{HIJING Pb+Pb 5.02 TeV}");
	sprintf(name,"#font[42]{HILOOSE   %0.1f<p_{T}<%0.1f GeV}",binPt[iP],binPt[iP+1]);
	tex->DrawLatex(0.175,0.765,name);
	leg->Draw();

	if(iOpt==0) sprintf(name,"PLOT/PbPb502_LOOSE_eff_Pt%d.pdf",iP);
	if(iOpt==1) sprintf(name,"PLOT/PbPb502_LOOSE_fak_Pt%d.pdf",iP);
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"PbPb502_LOOSE_eff_Pt%d",iP);
	if(iOpt==1) sprintf(name,"PbPb502_LOOSE_fak_Pt%d",iP);
	cOut->SetName(name);
	cOut->Write();

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

void Plot::readHist_CP(TFile* fIn, TH1D*& hIn, const char* hName, int iC, int iP)
{
	sprintf(name,"%s_Cent%d_Pt%d",hName,iC,iP);
	hIn = (TH1D*)fIn->Get(name);
}
