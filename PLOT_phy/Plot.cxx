#include "Plot.h"

Plot::Plot()
{
	initialize();
	execute();
	finalize();
}

Plot::~Plot()
{
	for(unsigned int iF=0; iF<NF; iF++)
	{
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				delete sts_c2_1sub[iF][iV][iP];
				delete sts_c4_1sub[iF][iV][iP];
				delete sts_c6_1sub[iF][iV][iP];
				delete sts_nc4_1sub[iF][iV][iP];
				delete sts_nc6_1sub[iF][iV][iP];
				delete sts_sc_1sub[iF][iV][iP];
				delete sts_nsc_1sub[iF][iV][iP];
				delete sts_ac_1sub[iF][iV][iP];
				delete sts_nac_1sub[iF][iV][iP];
				delete sts_isGauss_1sub[iF][iV][iP];
				delete sts_isPower_1sub[iF][iV][iP];
				delete sts_c2_3sub[iF][iV][iP];
				delete sts_c4_3sub[iF][iV][iP];
				delete sts_nc4_3sub[iF][iV][iP];
				delete sts_sc_3sub[iF][iV][iP];
				delete sts_nsc_3sub[iF][iV][iP];
				delete sts_ac_3sub[iF][iV][iP];
				delete sts_nac_3sub[iF][iV][iP];

				delete sys_c2_1sub[iF][iV][iP];
				delete sys_c4_1sub[iF][iV][iP];
				delete sys_c6_1sub[iF][iV][iP];
				delete sys_nc4_1sub[iF][iV][iP];
				delete sys_nc6_1sub[iF][iV][iP];
				delete sys_sc_1sub[iF][iV][iP];
				delete sys_nsc_1sub[iF][iV][iP];
				delete sys_ac_1sub[iF][iV][iP];
				delete sys_nac_1sub[iF][iV][iP];
				delete sys_isGauss_1sub[iF][iV][iP];
				delete sys_isPower_1sub[iF][iV][iP];
				delete sys_c2_3sub[iF][iV][iP];
				delete sys_c4_3sub[iF][iV][iP];
				delete sys_nc4_3sub[iF][iV][iP];
				delete sys_sc_3sub[iF][iV][iP];
				delete sys_nsc_3sub[iF][iV][iP];
				delete sys_ac_3sub[iF][iV][iP];
				delete sys_nac_3sub[iF][iV][iP];
			}
		}
	}
}

void Plot::execute()
{
  cout<<"execute..."<<endl;

	TFile* fOut = new TFile("PLOT/canvas.root","RECREATE");
	fOut->cd();
	vector<TGraphErrors*> vSts;
	vector<TGraphAsymmErrors*> vSys;

	// 0: mtd_c2
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSts.push_back(sts_c2_1sub[0][iV][0]);
		vSts.push_back(sts_c2_3sub[0][iV][0]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSts.push_back(sts_c2_1sub[0][iV][6]);
		vSts.push_back(sts_c2_3sub[0][iV][6]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSys.push_back(sys_c2_1sub[0][iV][0]);
		vSys.push_back(sys_c2_3sub[0][iV][0]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSys.push_back(sys_c2_1sub[0][iV][6]);
		vSys.push_back(sys_c2_3sub[0][iV][6]);
	}
	draw_graph_4by2(vSts,vSys,0);
	vSts.clear();
	vSys.clear();

	// 1: mtd_c4
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSts.push_back(sts_c4_1sub[0][iV][0]);
		vSts.push_back(sts_c4_3sub[0][iV][0]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSts.push_back(sts_c4_1sub[0][iV][6]);
		vSts.push_back(sts_c4_3sub[0][iV][6]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSys.push_back(sys_c4_1sub[0][iV][0]);
		vSys.push_back(sys_c4_3sub[0][iV][0]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSys.push_back(sys_c4_1sub[0][iV][6]);
		vSys.push_back(sys_c4_3sub[0][iV][6]);
	}
	draw_graph_4by2(vSts,vSys,1);
	vSts.clear();
	vSys.clear();

	// 2: mtd_nc4
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSts.push_back(sts_nc4_1sub[0][iV][0]);
		vSts.push_back(sts_nc4_3sub[0][iV][0]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSts.push_back(sts_nc4_1sub[0][iV][6]);
		vSts.push_back(sts_nc4_3sub[0][iV][6]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSys.push_back(sys_nc4_1sub[0][iV][0]);
		vSys.push_back(sys_nc4_3sub[0][iV][0]);
	}
	for(unsigned int iV=1; iV<NV; iV++)
	{
		vSys.push_back(sys_nc4_1sub[0][iV][6]);
		vSys.push_back(sys_nc4_3sub[0][iV][6]);
	}
	draw_graph_4by2(vSts,vSys,2);
	vSts.clear();
	vSys.clear();

	// 3: mtd_sc
	vSts.push_back(sts_sc_1sub[0][2][0]);
	vSts.push_back(sts_sc_3sub[0][2][0]);
	vSts.push_back(sts_sc_1sub[0][3][0]);
	vSts.push_back(sts_sc_3sub[0][3][0]);
	vSts.push_back(sts_ac_1sub[0][0][0]);
	vSts.push_back(sts_ac_3sub[0][0][0]);
	vSts.push_back(sts_ac_1sub[0][2][0]);
	vSts.push_back(sts_ac_3sub[0][2][0]);
	vSts.push_back(sts_sc_1sub[0][2][6]);
	vSts.push_back(sts_sc_3sub[0][2][6]);
	vSts.push_back(sts_sc_1sub[0][3][6]);
	vSts.push_back(sts_sc_3sub[0][3][6]);
	vSts.push_back(sts_ac_1sub[0][0][6]);
	vSts.push_back(sts_ac_3sub[0][0][6]);
	vSts.push_back(sts_ac_1sub[0][2][6]);
	vSts.push_back(sts_ac_3sub[0][2][6]);

	vSys.push_back(sys_sc_1sub[0][2][0]);
	vSys.push_back(sys_sc_3sub[0][2][0]);
	vSys.push_back(sys_sc_1sub[0][3][0]);
	vSys.push_back(sys_sc_3sub[0][3][0]);
	vSys.push_back(sys_ac_1sub[0][0][0]);
	vSys.push_back(sys_ac_3sub[0][0][0]);
	vSys.push_back(sys_ac_1sub[0][2][0]);
	vSys.push_back(sys_ac_3sub[0][2][0]);
	vSys.push_back(sys_sc_1sub[0][2][6]);
	vSys.push_back(sys_sc_3sub[0][2][6]);
	vSys.push_back(sys_sc_1sub[0][3][6]);
	vSys.push_back(sys_sc_3sub[0][3][6]);
	vSys.push_back(sys_ac_1sub[0][0][6]);
	vSys.push_back(sys_ac_3sub[0][0][6]);
	vSys.push_back(sys_ac_1sub[0][2][6]);
	vSys.push_back(sys_ac_3sub[0][2][6]);
	draw_graph_4by2(vSts,vSys,3);
	vSts.clear();
	vSys.clear();

	const int selPt[4] = {0,2,4,6};
	// 0: phy_c2
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_c2_3sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_c2_3sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,0);
		vSts.clear();
		vSys.clear();
	}

	// 1: phy_c4
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_c4_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_c4_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,1);
		vSts.clear();
		vSys.clear();
	}

	// 2: phy_nc4
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_nc4_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_nc4_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,2);
		vSts.clear();
		vSys.clear();
	}

	// 3: phy_c6
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_c6_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_c6_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,3);
		vSts.clear();
		vSys.clear();
	}

	// 4: phy_nc6
	for(unsigned int iV=1; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_nc6_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_nc6_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,4);
		vSts.clear();
		vSys.clear();
	}

	// 5: phy_sc
	for(unsigned int iV=0; iV<4; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_sc_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_sc_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,5);
		vSts.clear();
		vSys.clear();
	}

	// 6: phy_nsc
	for(unsigned int iV=0; iV<4; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_nsc_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_nsc_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,6);
		vSts.clear();
		vSys.clear();
	}

	// 7: phy_ac
	for(unsigned int iV=0; iV<3; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_ac_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_ac_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,7);
		vSts.clear();
		vSys.clear();
	}

	// 8: phy_nac
	for(unsigned int iV=0; iV<3; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_nac_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_nac_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,8);
		vSts.clear();
		vSys.clear();
	}

	// 9: phy_isGauss
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_isGauss_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_isGauss_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,9);
		vSts.clear();
		vSys.clear();
	}

	// 10: phy_isPower
	for(unsigned int iV=2; iV<NV; iV++)
	{
		for(unsigned int iF=0; iF<NF; iF++)
		{
			for(unsigned int iP=0; iP<4; iP++)
			{
				vSts.push_back(sts_isPower_1sub[iF][iV][selPt[iP]]);
				vSys.push_back(sys_isPower_1sub[iF][iV][selPt[iP]]);
			}
		}
		draw_graph_4by3(vSts,vSys,iV,10);
		vSts.clear();
		vSys.clear();
	}

	fOut->Close();
}

void Plot::initialize()
{
  cout<<"initialize..."<<endl;
	
	TFile* fIn[NF][5];
	for(unsigned int iF=0; iF<NF; iF++)
	{
		for(unsigned int iB=0; iB<5; iB++)
		{
			if(iF==0) sprintf(name,"../MAIN_binCent/OUTPUT/Phase4/hist_PbPb502_binCent_bin%d.root",iB);
			if(iF==1) sprintf(name,"../MAIN_binFCal/OUTPUT/Phase4/hist_PbPb502_binFCal_bin%d.root",iB);
			if(iF==2) sprintf(name,"../MAIN_binNch/OUTPUT/Phase4/hist_PbPb502_binNch_bin%d.root",iB);
			fIn[iF][iB] = new TFile(name,"READ");
		}
		for(unsigned int iV=0; iV<NV; iV++)
		{
			for(unsigned int iP=0; iP<NP; iP++)
			{
				int iB;
				int def = 2;
				iB = def;
				readHist_VP(fIn[iF][iB],sts_c2_1sub[iF][iV][iP],"sts_c2_1sub",iV,iP);

				if(iV==1) iB = 4;
				else iB = def;
				readHist_VP(fIn[iF][iB],sts_c4_1sub[iF][iV][iP],"sts_c4_1sub",iV,iP);
				iB = def;

				readHist_VP(fIn[iF][iB],sts_c6_1sub[iF][iV][iP],"sts_c6_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_nc4_1sub[iF][iV][iP],"sts_nc4_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_nc6_1sub[iF][iV][iP],"sts_nc6_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_sc_1sub[iF][iV][iP],"sts_sc_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_nsc_1sub[iF][iV][iP],"sts_nsc_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_ac_1sub[iF][iV][iP],"sts_ac_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_nac_1sub[iF][iV][iP],"sts_nac_1sub",iV,iP);

				iB = 3;
				readHist_VP(fIn[iF][iB],sts_isGauss_1sub[iF][iV][iP],"sts_isGauss_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_isPower_1sub[iF][iV][iP],"sts_isPower_1sub",iV,iP);
				iB = def;

				readHist_VP(fIn[iF][iB],sts_c2_3sub[iF][iV][iP],"sts_c2_3sub",iV,iP);

				if(iV==1) iB = 4;
				else iB = def;
				readHist_VP(fIn[iF][iB],sts_c4_3sub[iF][iV][iP],"sts_c4_3sub",iV,iP);
				iB = def;

				readHist_VP(fIn[iF][iB],sts_nc4_3sub[iF][iV][iP],"sts_nc4_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_sc_3sub[iF][iV][iP],"sts_sc_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_nsc_3sub[iF][iV][iP],"sts_nsc_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_ac_3sub[iF][iV][iP],"sts_ac_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sts_nac_3sub[iF][iV][iP],"sts_nac_3sub",iV,iP);
				
				// ---

				iB = def;
				readHist_VP(fIn[iF][iB],sys_c2_1sub[iF][iV][iP],"sys_c2_1sub",iV,iP);

				if(iV==1) iB = 4;
				else iB = def;
				readHist_VP(fIn[iF][iB],sys_c4_1sub[iF][iV][iP],"sys_c4_1sub",iV,iP);
				iB = def;

				readHist_VP(fIn[iF][iB],sys_c6_1sub[iF][iV][iP],"sys_c6_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_nc4_1sub[iF][iV][iP],"sys_nc4_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_nc6_1sub[iF][iV][iP],"sys_nc6_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_sc_1sub[iF][iV][iP],"sys_sc_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_nsc_1sub[iF][iV][iP],"sys_nsc_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_ac_1sub[iF][iV][iP],"sys_ac_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_nac_1sub[iF][iV][iP],"sys_nac_1sub",iV,iP);

				iB = 3;
				readHist_VP(fIn[iF][iB],sys_isGauss_1sub[iF][iV][iP],"sys_isGauss_1sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_isPower_1sub[iF][iV][iP],"sys_isPower_1sub",iV,iP);
				iB = def;

				readHist_VP(fIn[iF][iB],sys_c2_3sub[iF][iV][iP],"sys_c2_3sub",iV,iP);

				if(iV==1) iB = 4;
				else iB = def;
				readHist_VP(fIn[iF][iB],sys_c4_3sub[iF][iV][iP],"sys_c4_3sub",iV,iP);
				iB = def;

				readHist_VP(fIn[iF][iB],sys_nc4_3sub[iF][iV][iP],"sys_nc4_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_sc_3sub[iF][iV][iP],"sys_sc_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_nsc_3sub[iF][iV][iP],"sys_nsc_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_ac_3sub[iF][iV][iP],"sys_ac_3sub",iV,iP);
				readHist_VP(fIn[iF][iB],sys_nac_3sub[iF][iV][iP],"sys_nac_3sub",iV,iP);
			}
		}
	}
}

void Plot::finalize()
{
	cout<<"finalize..."<<endl;
}

void Plot::draw_graph_4by2(vector<TGraphErrors*> vSts, vector<TGraphAsymmErrors*> vSys, int iOpt)
{
	int NG = vSts.size();
	TGraphErrors* gSts[20];
	TGraphAsymmErrors* gSys[20];
	for(int iG=0; iG<NG; iG++)
	{
		gSts[iG] = (TGraphErrors*)vSts.at(iG)->Clone("gSts");
		gSys[iG] = (TGraphAsymmErrors*)vSys.at(iG)->Clone("gSys");
		styleGraph(gSts[iG],iG%2);
		styleGraph(gSys[iG],iG%2);
		gSys[iG]->SetFillColor(mC[iG%2]);
		gSys[iG]->SetFillStyle(3001);
		for(int i=0; i<gSys[iG]->GetN(); i++)
		{
			gSys[iG]->SetPointEXhigh(i,1);
			gSys[iG]->SetPointEXlow(i,1);
		}
	}

	TLatex* tex = new TLatex();
	tex->SetTextFont(42);
	tex->SetTextSize(0.045);
	tex->SetTextAlign(12);
	tex->SetNDC(1);
	TLine* lin = new TLine();
	lin->SetLineColor(kGray);
	lin->SetLineStyle(2);
	lin->SetLineWidth(1);
	TLegend* leg = new TLegend(0.25,0.125,0.95,0.25);
	leg->SetTextSize(0.05);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(2);
	leg->AddEntry(gSts[0],"#font[42]{standard}","p");
	leg->AddEntry(gSts[1],"#font[42]{3-subevent}","p");

	double xMin =  0;
	double xMax = 80;
	double yMin[NG/2];
	double yMax[NG/2];
	for(int iG=0; iG<NG/2; iG++)
	{
		yMin[iG] = 1;
		yMax[iG] = -1;
		getYrange(gSts[iG*2],yMin[iG],yMax[iG],0);
		getYrange(gSts[iG*2+1],yMin[iG],yMax[iG],0);
		double diff = yMax[iG]-yMin[iG];
		yMax[iG] += 0.5*diff;
		yMin[iG] -= 0.5*diff;
	}

	TH1D* hAxis = new TH1D("hAxis","",1000,xMin,xMax);
	for(int i=0; i<1000; i++) hAxis->SetBinContent(i+1,1E9);
	styleGraph(hAxis,0);

	TCanvas* cOut = new TCanvas("cOut","",1200,600);
	cOut->Divide(4,2,1E-4,1E-4);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	for(int iC=0; iC<8; iC++)
	{
		cOut->cd(iC+1);
		gPad->SetTicks(1,1);
		gPad->SetTopMargin(0.075);
		gPad->SetBottomMargin(0.1);
		gPad->SetLeftMargin(0.2);
		gPad->SetRightMargin(0.025);
		hAxis->GetXaxis()->SetTitle("Centrality / %");
		hAxis->GetXaxis()->SetTitleOffset(1.75);
		hAxis->GetXaxis()->SetRangeUser(xMin,xMax);
		if(iOpt==0) sprintf(name,"c_{%d}{2}",iC%4+1);
		if(iOpt==1) sprintf(name,"c_{%d}{4}",iC%4+1);
		if(iOpt==2) sprintf(name,"#hat{c}_{%d}{4}",iC%4+1);
		if(iOpt==3)
		{
			if(iC%4==0) sprintf(name,"SC(2,3)");
			if(iC%4==1) sprintf(name,"SC(2,4)");
			if(iC%4==2) sprintf(name,"ASC(1,1,2)");
			if(iC%4==3) sprintf(name,"ASC(2,2,4)");
		}
		hAxis->GetYaxis()->SetTitle(name);
		hAxis->GetYaxis()->SetTitleOffset(4);
		hAxis->GetYaxis()->SetRangeUser(yMin[iC],yMax[iC]);
		hAxis->DrawCopy();
		gSys[2*iC]->Draw("2");
		gSys[2*iC+1]->Draw("2");
		gSts[2*iC]->Draw("P");
		gSts[2*iC+1]->Draw("P");
		tex->DrawLatex(0.25,0.875,"#font[72]{ATLAS} #font[62]{Internal}");
		tex->DrawLatex(0.25,0.82,"#font[42]{Pb+Pb #sqrt{s_{NN}}=5.02 TeV}");
		if(iC/4==0) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[0],maxPt[0]);
		else sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[6],maxPt[6]);
		tex->DrawLatex(0.25,0.76,name);
		if(yMax[iC]>0 && yMin[iC]<0) lin->DrawLine(xMin,0,xMax,0);
		leg->Draw();
	}

	if(iOpt==0) sprintf(name,"PLOT/mtd_c2.pdf");
	if(iOpt==1) sprintf(name,"PLOT/mtd_c4.pdf");
	if(iOpt==2) sprintf(name,"PLOT/mtd_nc4.pdf");
	if(iOpt==3) sprintf(name,"PLOT/mtd_sc.pdf");
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"mtd_c2");
	if(iOpt==1) sprintf(name,"mtd_c4");
	if(iOpt==2) sprintf(name,"mtd_nc4");
	if(iOpt==3) sprintf(name,"mtd_sc");
	cOut->SetName(name);
	cOut->Write();

	delete hAxis;
	delete cOut;
}

void Plot::draw_graph_4by3(vector<TGraphErrors*> vSts, vector<TGraphAsymmErrors*> vSys, int iV, int iOpt)
{
	int NG = vSts.size();
	TGraphErrors* gSts[20];
	TGraphAsymmErrors* gSys[20];
	for(int iG=0; iG<NG; iG++)
	{
		gSts[iG] = (TGraphErrors*)vSts.at(iG)->Clone("gSts");
		gSys[iG] = (TGraphAsymmErrors*)vSys.at(iG)->Clone("gSys");
		styleGraph(gSts[iG],0);
		styleGraph(gSys[iG],0);
		gSys[iG]->SetFillColor(mC[0]);
		gSys[iG]->SetFillStyle(3001);
		for(int i=0; i<gSys[iG]->GetN(); i++)
		{
			if(iG/4==0)
			{
				gSys[iG]->SetPointEXhigh(i,1);
				gSys[iG]->SetPointEXlow(i,1);
			}
			else if(iG/4==1)
			{
				gSys[iG]->SetPointEXhigh(i,0.1);
				gSys[iG]->SetPointEXlow(i,0.1);
			}
			else
			{
				gSys[iG]->SetPointEXhigh(i,50);
				gSys[iG]->SetPointEXlow(i,50);
			}
		}
	}

	TLatex* tex = new TLatex();
	tex->SetTextFont(42);
	tex->SetTextSize(0.045);
	tex->SetTextAlign(12);
	tex->SetNDC(1);
	TLine* lin = new TLine();
	lin->SetLineColor(kGray);
	lin->SetLineStyle(2);
	lin->SetLineWidth(1);
	TLegend* leg = new TLegend(0.25,0.125,0.95,0.25);
	leg->SetTextSize(0.05);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetNColumns(2);
	if(iOpt==0) leg->AddEntry(gSts[0],"#font[42]{3-subevent}","p");
	else leg->AddEntry(gSts[0],"#font[42]{standard}","p");

	double xMin =  0;
	double xMax = 80;
	double yMin[NG];
	double yMax[NG];
	for(int iG=0; iG<NG; iG++)
	{
		yMin[iG] = 1;
		yMax[iG] = -1;
		getYrange(gSts[iG],yMin[iG],yMax[iG],iG/4);
		double diff = yMax[iG]-yMin[iG];
		yMax[iG] += 0.5*diff;
		yMin[iG] -= 0.5*diff;
	}

	TH1D* hAxis = new TH1D("hAxis","",8000,0,4000);
	for(int i=0; i<8000; i++) hAxis->SetBinContent(i+1,1E9);
	styleGraph(hAxis,0);

	TCanvas* cOut = new TCanvas("cOut","",1200,900);
	cOut->Divide(4,3,1E-4,1E-4);
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.0001);
	for(int iC=0; iC<12; iC++)
	{
		cOut->cd(iC+1);
		gPad->SetTicks(1,1);
		gPad->SetTopMargin(0.075);
		gPad->SetBottomMargin(0.125);
		gPad->SetLeftMargin(0.2);
		gPad->SetRightMargin(0.025);
		if(iC/4==0) hAxis->GetXaxis()->SetTitle("Centrality / %");
		if(iC/4==1) hAxis->GetXaxis()->SetTitle("FCal E_{T} / TeV");
		if(iC/4==2) hAxis->GetXaxis()->SetTitle("N_{ch}^{rec}");
		hAxis->GetXaxis()->SetTitleOffset(3.25);
		if(iC/4==0) {xMin = 0; xMax = 80;}
		if(iC/4==1) {xMin = 0; xMax = 4.9;}
		if(iC/4==2) {xMin = 0; xMax = 3500;}
		hAxis->GetXaxis()->SetRangeUser(xMin,xMax);
		if(iOpt==0) sprintf(name,"c_{%d}{2}",iV);
		if(iOpt==1) sprintf(name,"c_{%d}{4}",iV);
		if(iOpt==2) sprintf(name,"#hat{c}_{%d}{4}",iV);
		if(iOpt==3) sprintf(name,"c_{%d}{6}",iV);
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
		hAxis->GetYaxis()->SetTitle(name);
		hAxis->GetYaxis()->SetTitleOffset(6);
		hAxis->GetYaxis()->SetRangeUser(yMin[iC],yMax[iC]);
		if(iOpt==9 || iOpt==10) hAxis->GetYaxis()->SetRangeUser(0.75,1.25);
		hAxis->DrawCopy();
		//if(iC/4==0) reverseXaxis(hAxis,xMin,xMax,yMin[iC],1);
		//else reverseXaxis(hAxis,xMin,xMax,yMin[iC],0);
		gSys[iC]->Draw("2");
		gSts[iC]->Draw("P");
		tex->DrawLatex(0.25,0.875,"#font[72]{ATLAS} #font[62]{Internal}");
		tex->DrawLatex(0.25,0.82,"#font[42]{Pb+Pb #sqrt{s_{NN}}=5.02 TeV}");
		if(iC%4==0) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[0],maxPt[0]);
		if(iC%4==1) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[2],maxPt[2]);
		if(iC%4==2) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[4],maxPt[4]);
		if(iC%4==3) sprintf(name,"#font[42]{%.1f<p_{T}<%.1f GeV}",minPt[6],maxPt[6]);
		tex->DrawLatex(0.25,0.76,name);
		if(iOpt!=9 && iOpt!=10)
		{
			if(yMax[iC]>0 && yMin[iC]<0) lin->DrawLine(xMin,0,xMax,0);
		}
		else lin->DrawLine(xMin,1,xMax,1);
		leg->Draw();
	}

	if(iOpt==0) sprintf(name,"PLOT/phy_c2_Har%d.pdf",iV);
	if(iOpt==1) sprintf(name,"PLOT/phy_c4_Har%d.pdf",iV);
	if(iOpt==2) sprintf(name,"PLOT/phy_nc4_Har%d.pdf",iV);
	if(iOpt==3) sprintf(name,"PLOT/phy_c6_Har%d.pdf",iV);
	if(iOpt==4) sprintf(name,"PLOT/phy_nc6_Har%d.pdf",iV);
	if(iOpt==5) sprintf(name,"PLOT/phy_sc_Har%d.pdf",iV);
	if(iOpt==6) sprintf(name,"PLOT/phy_nsc_Har%d.pdf",iV);
	if(iOpt==7) sprintf(name,"PLOT/phy_ac_Har%d.pdf",iV);
	if(iOpt==8) sprintf(name,"PLOT/phy_nac_Har%d.pdf",iV);
	if(iOpt==9) sprintf(name,"PLOT/phy_isGauss_Har%d.pdf",iV);
	if(iOpt==10) sprintf(name,"PLOT/phy_isPower_Har%d.pdf",iV);
	cOut->Print(name);
	if(iOpt==0) sprintf(name,"phy_c2_Har%d",iV);
	if(iOpt==1) sprintf(name,"phy_c4_Har%d",iV);
	if(iOpt==2) sprintf(name,"phy_nc4_Har%d",iV);
	if(iOpt==3) sprintf(name,"phy_c6_Har%d",iV);
	if(iOpt==4) sprintf(name,"phy_nc6_Har%d",iV);
	if(iOpt==5) sprintf(name,"phy_sc_Har%d",iV);
	if(iOpt==6) sprintf(name,"phy_nsc_Har%d",iV);
	if(iOpt==7) sprintf(name,"phy_ac_Har%d",iV);
	if(iOpt==8) sprintf(name,"phy_nac_Har%d",iV);
	if(iOpt==9) sprintf(name,"phy_isGauss_Har%d",iV);
	if(iOpt==10) sprintf(name,"phy_isPower_Har%d",iV);
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
	hIn->SetLineWidth(1);
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
		if(iOpt==0 && x>60) continue;
		if(iOpt==1 && (x>4.2 || x<=0.5)) continue;
		if(iOpt==2 && (x>3300 || x<=250)) continue;
		if(y<yMin) yMin = y;
		if(y>yMax) yMax = y;
	}
}

void Plot::readHist_VP(TFile* fIn, TGraphErrors*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphErrors*)fIn->Get(name);
}

void Plot::readHist_VP(TFile* fIn, TGraphAsymmErrors*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = (TGraphAsymmErrors*)fIn->Get(name);
}

void Plot::reverseXaxis(TH1D* hIn, double xMin, double xMax, double y, bool isRev)
{
	double xMinNew = 0;
	double xMaxNew = 0;
	if(isRev)
	{
		xMinNew = xMax;
		xMaxNew = xMin;
	}
	else
	{
		xMinNew = xMin;
		xMaxNew = xMax;
	}
	hIn->GetXaxis()->SetLabelOffset(999);
	hIn->GetXaxis()->SetTickLength(0);
	gPad->Update();
	TGaxis *newAxis = new TGaxis(xMinNew, y, xMaxNew, y, xMin, xMax, 510, "-");
	newAxis->SetLabelFont(43);
	newAxis->SetLabelSize(15);
	newAxis->SetLabelOffset(-0.03);
	newAxis->Draw();
}
