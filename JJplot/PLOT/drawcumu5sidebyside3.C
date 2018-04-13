#include <iostream>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TLegend.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TView3D.h>
#include <TTree.h>
#include <TLatex.h>
#include <TROOT.h>
#include <stdio.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"


void Divide(TCanvas* can,int x,int y,float marx, float mary){
	double xcoor1[10][10], xcoor2[10][10], ycoor1[10][10], ycoor2[10][10];
	Double_t xlow,ylow,xup,yup;
	double ratx[]={1,1,1,1,1,1,1,1};
	double raty[]={1,1,1,1,1,1,1,1};
	double fracx[10], fracy[10];

	double xsli=0,ysli=0;//for boundary
	//define the slice size;
	ratx[0]   +=marx;
	raty[y-1] +=mary;
	for(int i=0;i<x;i++){
		xsli+=ratx[i];
	}
	for(int i=0;i<y;i++){
		ysli+=raty[i];
	}
	fracx[0]=0;  fracy[0]=1;
	for(int i=1;i<=x;i++){
		fracx[i]=fracx[i-1]+ratx[i-1]/xsli;
	}
	for(int i=1;i<=y;i++){
		fracy[i]=fracy[i-1]-raty[i-1]/ysli;
	}
	//rescale
	double scal=0.995;
	for(int i=0;i<=x;i++){
		fracx[i]= fracx[i]*scal+(1-scal)*(0.5-fracx[i]);
	}
	for(int i=0;i<=y;i++){
		fracy[i]= fracy[i]*scal+(1-scal)*(0.5-fracy[i]);
	}
	can->cd();
	can->Divide(x,y);
	int count=1;
	for(int i=0;i<y;i++){
		for(int j=0;j<x;j++){
			can->cd(count);
			count++;
			xlow = fracx[j];      xup = fracx[j+1];
			ylow = fracy[i+1];    yup = fracy[i];
			xcoor1[i][j] = xlow;      xcoor2[i][j] = xup;
			ycoor1[i][j] = ylow;      ycoor2[i][j] = yup;
			//cout<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
			gPad->SetPad(xlow,ylow,xup,yup);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			gPad->SetTopMargin(0.06);
			gPad->SetBottomMargin(0.01);
			if(j==0){
				gPad->SetLeftMargin((marx+0.15)/ratx[0]);
			}
			if(i==y-1){
				gPad->SetBottomMargin(mary/raty[y-1]);
			}
		}
	}
}

void Divide1(TCanvas* can,int x,int y,float marx, float mary){
	double xcoor1[10][10], xcoor2[10][10], ycoor1[10][10], ycoor2[10][10];
	Double_t xlow,ylow,xup,yup;
	double ratx[]={1,1,1,1,1,1,1,1};
	double raty[]={1,0.5,1,0.5,1,0.5,1,0.3};
	double fracx[10], fracy[10];

	double xsli=0,ysli=0;//for boundary
	//define the slice size;
	ratx[0]   +=marx;
	raty[y-1] +=mary;
	for(int i=0;i<x;i++){
		xsli+=ratx[i];
	}
	for(int i=0;i<y;i++){
		ysli+=raty[i];
	}
	fracx[0]=0;  fracy[0]=1;
	for(int i=1;i<=x;i++){
		fracx[i]=fracx[i-1]+ratx[i-1]/xsli;
	}
	for(int i=1;i<=y;i++){
		fracy[i]=fracy[i-1]-raty[i-1]/ysli;
	}
	//rescale
	double scal=0.995;
	for(int i=0;i<=x;i++){
		fracx[i]= fracx[i]*scal+(1-scal)*(0.5-fracx[i]);
	}
	for(int i=0;i<=y;i++){
		fracy[i]= fracy[i]*scal+(1-scal)*(0.5-fracy[i]);
	}
	can->cd();
	can->Divide(x,y);
	int count=1;
	for(int i=0;i<y;i++){
		for(int j=0;j<x;j++){
			can->cd(count);
			count++;
			xlow = fracx[j];      xup = fracx[j+1];
			ylow = fracy[i+1];    yup = fracy[i];
			xcoor1[i][j] = xlow;      xcoor2[i][j] = xup;
			ycoor1[i][j] = ylow;      ycoor2[i][j] = yup;
			//cout<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
			gPad->SetPad(xlow,ylow,xup,yup);
			gPad->SetLeftMargin(0.15);
			gPad->SetRightMargin(0.01);
			gPad->SetTopMargin(0.01);
			gPad->SetBottomMargin(0.01);
			if(j==0){
				gPad->SetLeftMargin((marx+0.15)/ratx[0]);
			}
			if(i==y-1){
				gPad->SetBottomMargin(mary/raty[y-1]);
			}
		}
	}
}
void setstyle(TGraph*h){
	h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
	h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
	h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
	h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
	h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
	//h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
	//h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}
void setstyle1(TGraph*h){
	h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
	h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(15);
	h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(13);
	h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(15);
	h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(13);
	//h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
	//h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}


void setstyle(TH1*h){
	h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
	h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
	h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
	h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
	h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
	//h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
	//h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}
void setstyle1(TH1*h){
	h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
	h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(15);
	h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(13);
	h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(15);
	h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(13);
	//h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
	//h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}
void setgrerror(TGraph*h1,double val){
	int N= h1->GetN();
	double *y=h1->GetY(); 
	if(h1->InheritsFrom("TGraphErrors")){
		double* ey1u = ((TGraphErrors*)h1)->GetEY(); 
		for(int i=0;i<N;i++){
			ey1u[i]=y[i]*val;
		}
	}else{
		double* ey1u = ((TGraphErrors*)h1)->GetEYhigh(); 
		double* ey1l = ((TGraphErrors*)h1)->GetEYlow(); 
		for(int i=0;i<N;i++){
			ey1u[i]=y[i]*val;
			ey1l[i]=y[i]*val;
		}
	}
}
TGraph* grratio(TGraph*h1,TGraph*h2, int opt=1){
	char n[200];  
	sprintf(n,"R%s",h1->GetName());
	TGraph *h3 = (TGraph*)h1->Clone(); h3->SetName(n);h3->SetTitle(n);

	double *x1 = h1->GetX(); 
	double *y1 = h1->GetY();  double *ty2 = h2->GetY();  double *y3 = h3->GetY();
	double *ey1u,*ey2u ,*ey3u,*ey1l,*ey2l ,*ey3l;
	double y2[1000];
	int type=0;
	if(h1->InheritsFrom("TGraphErrors")){
		ey1u = ((TGraphErrors*)h1)->GetEY(); 
		ey2u = ((TGraphErrors*)h2)->GetEY(); 
		ey3u =  ((TGraphErrors*)h3)->GetEY(); 
	}else{
		type=1;
		ey1u  = ((TGraphAsymmErrors*)h1)->GetEYhigh(); 
		ey2u = ((TGraphAsymmErrors*)h2)->GetEYhigh();// 
		ey3u  = ((TGraphAsymmErrors*)h3)->GetEYhigh(); 

		ey1l  = ((TGraphAsymmErrors*)h1)->GetEYlow(); 
		ey2l = ((TGraphAsymmErrors*)h2)->GetEYlow(); //
		ey3l  = ((TGraphAsymmErrors*)h3)->GetEYlow(); 
	}
	if(opt==0){
		int N= h1->GetN();
		for(int i=0;i<N;i++){
			y2[i] = h2->Eval(x1[i]);//
			y3[i] = ey3u[i]=0;
			if(y2[i]){
				y3[i] = y1[i]/y2[i];
				ey3u[i] = ey1u[i]/y2[i];
				if(type) ey3l[i] = ey1l[i]/y2[i];
			}
		}
	}else{//reverse
		int N= h1->GetN();
		for(int i=0;i<N;i++){
			y2[i] = h2->Eval(x1[i]);//
			y3[i] = ey3u[i]=0;
			if(y1[i]){
				y3[i] = y2[i]/y1[i];
				ey3u[i] = ey1u[i]/y2[i];
				if(type) ey3l[i] = ey1l[i]/y2[i];
			}
		}
	}
	return h3;
}
TGraphAsymmErrors *grconvert(TGraph*gr){
	int N = gr->GetN();
	double *xin = gr->GetX();
	double *yin = gr->GetY();
	double *yinel,*yineh;

	if(gr->InheritsFrom("TGraphErrors")){
		yinel = ((TGraphErrors*)gr)->GetEY();
		yineh = ((TGraphErrors*)gr)->GetEY();
	}else{
		yinel = ((TGraphAsymmErrors*)gr)->GetEYlow();
		yineh = ((TGraphAsymmErrors*)gr)->GetEYhigh();
	}
	double y[100];
	double ye1[100];
	double ye2[100];
	double val1,val2;
	double val3,val4;
	for(int i=0;i<N;i++){
		y[i] = pow(fabs(-yin[i]),0.25);
		if(yin[i]>0) y[i] *=-1;
		val1 = -yin[i]+yinel[i];
		val2 = -yin[i]-yineh[i];
		val3 = pow(fabs(val1),0.25); if(val1<0) val3*=-1;
		val4 = pow(fabs(val2),0.25); if(val2<0) val4*=-1;
		ye1[i] = val3-y[i];
		ye2[i] = y[i]-val4;    
	}
	return new TGraphAsymmErrors(N,xin,y,0,0,ye2,ye1);
}
void scalgr(TGraph*gr,double scal){

	int N = gr->GetN();
	double *X = gr->GetX();
	double *Y = gr->GetY();
	double *Ye1,*Ye2;
	double ty=1;
	if(gr->InheritsFrom("TGraphErrors")){
		Ye1 = ((TGraphErrors*) gr)->GetEY();
	}else{
		ty=0;
		Ye1 = ((TGraphAsymmErrors*) gr)->GetEYlow();
		Ye2 = ((TGraphAsymmErrors*) gr)->GetEYhigh();
	}
	double valtmp;
	for(int i=0;i<N;i++){
		if(scal>0){
			Y[i]*=scal;    Ye1[i]*=scal;    if(ty==0) Ye2[i]*=scal;
		}else{
			Y[i]*=scal;    Ye1[i]*=-scal;  if(ty==0) Ye2[i]*=-scal;
			if(ty==0){
				valtmp = Ye1[i]; Ye1[i]=Ye2[i]; Ye2[i]=valtmp;
			}
		}
	}
}
double Knee[3]={4.1,2800,1};
void grcent(TGraph*gr, int type){  
	int N = gr->GetN();
	double *X = gr->GetX();
	switch (type) {
		case 0:
			for(int i=0;i<N;i++) X[i]/=Knee[0]; break;//ET/4.1
		case 1:
			for(int i=0;i<N;i++) X[i]/=Knee[1]; break;//Nch/2800
		case 2:
			for(int i=0;i<N;i++) X[i]=1-X[i]/100.; break; //1-Cent
			//for(int i=0;i<N;i++) X[i]=X[i]; break; //1-Cent // FIX
		case 3:
			for(int i=0;i<N;i++) X[i]/=Knee[1]; break;//ET map Nch then divide by 2800
		case 4:
			for(int i=0;i<N;i++) X[i]/=Knee[0]; break;//Nch map ET then divide by 4.1
		case 5: 
			; //Cent
	}
}
void grmap(TGraph*gr1, TGraph*gr2){//mapping
	int N=gr1->GetN();
	double *X=gr1->GetX();
	for(int ib=0;ib<N;ib++){
		X[ib]=gr2->Eval(X[ib]);
	}
}

void comb(TGraph*gr,TGraph*gr1){  
	int N = gr->GetN();
	double *X = gr->GetX();
	double *Y = gr->GetY();
	double *Ye1,*Ye2;
	int N1 = gr1->GetN();
	double *X1 = gr1->GetX();
	double *Y1 = gr1->GetY();
	double *Y1e1,*Y1e2;
	//int st1=15,st2=6,np1 =25,np2=10,nrm=15;
	int st1=10,st2=4,np1 =30,np2=12,nrm=18;

	if(N==16) {
		//st1=6; st2=3; np1 =10; np2=5;nrm=5;
		st1=4; st2=2; np1 =12; np2=6;nrm=6;
	}
	for(int ip=0;ip<nrm;ip++){
		gr->RemovePoint(st1);
	}
	double ty=1;
	if(gr->InheritsFrom("TGraphErrors")){
		Ye1 = ((TGraphErrors*) gr)->GetEY();
		Y1e1 = ((TGraphErrors*) gr1)->GetEY();
	}else{
		ty=0;
		Ye1 = ((TGraphAsymmErrors*) gr)->GetEYlow();
		Ye2 = ((TGraphAsymmErrors*) gr)->GetEYhigh();
		Y1e1 = ((TGraphAsymmErrors*) gr1)->GetEYlow();
		Y1e2 = ((TGraphAsymmErrors*) gr1)->GetEYhigh();
	}

	for(int i=0;i<np2;i++){
		X[i+st1] = X1[i+st2];
		Y[i+st1] = Y1[i+st2];
		Ye1[i+st1] =  Y1e1[i+st2];
		if(ty==0){
			Ye2[i+st1] =  Y1e2[i+st2];
		}
	}
}

TGraphErrors* h2gr(TH1*h){
	int N = h->GetNbinsX();
	double X[20000],Y[20000],YE[20000];
	for(int ib=0;ib<N;ib++){
		X[ib]  = h->GetBinCenter(ib+1);
		Y[ib]  = h->GetBinContent(ib+1);
		YE[ib] = h->GetBinError(ib+1);
	}
	char n[200];
	sprintf(n,"g%s",h->GetName()); TGraphErrors*gr = new TGraphErrors(N,X,Y,0,YE);
	gr->SetName(n);  return  gr;
}
double m1=0,m2=0;
double ran1=100, ran2=3500;
void GetMax(TGraph*gr){
	m1 =-0.5; m2=0.5;
	int N = gr->GetN();
	double *X = gr->GetX();
	double *Y = gr->GetY();
	for(int i=0;i<N;i++){
		if(X[i]<ran1||X[i]>ran2) continue;

		double val = Y[i];
		if(fabs(val)<5){
			if(Y[i]>m1) m1 = Y[i];
			if(Y[i]<m2) m2 = Y[i];
		}
		//cout<<Y[i]<<endl;
	}
}
void grpow(TGraph*gr,double n){
	int N = gr->GetN();
	double *X = gr->GetX();
	double *Y = gr->GetY();
	double *Y1, *Y2;
	int type=0;
	if(gr->InheritsFrom("TGraphErrors")){
		Y1 = ((TGraphErrors*) gr)->GetEY();
	}else{
		type=1;
		Y1 = ((TGraphAsymmErrors*) gr)->GetEYlow();
		Y2 = ((TGraphAsymmErrors*) gr)->GetEYhigh();
	}
	for(int ip=0;ip<N;ip++){

		double val = Y[ip], val1=Y[ip]; int sign;
		sign=1; if(val<0) sign=-1;  
		Y[ip] = pow(fabs(val),n)*sign;
		if(Y1[ip]<0) cout<<"HI"<<endl;
		val=val1-Y1[ip];      sign=1; if(val<0) sign=-1; Y1[ip] = Y[ip]-pow(fabs(val),n)*sign;

		if(type){
			val=val1+Y2[ip];    sign=1; if(val<0) sign=-1; Y2[ip] = pow(fabs(val),n)*sign-Y[ip];
		}
	}
}
TGraphAsymmErrors* toasygr(TGraphErrors*gr){
	int N = gr->GetN();
	double *X = gr->GetX();
	double *Y = gr->GetY();
	double *YE = gr->GetEY();
	TGraphAsymmErrors*gtmp = new TGraphAsymmErrors(N,X,Y,0,0,YE,YE);
	gtmp->SetName(gr->GetName());
	return gtmp;
}


TGraphAsymmErrors* grcorr(TGraph*g1,TGraph*g2){
	int N= g1->GetN();  if(N!=g2->GetN()) cout<<"Bad"<<endl;

	char n[200];    sprintf(n,"%s_cor",g1->GetName());  TGraphAsymmErrors* g3 = (TGraphAsymmErrors*)g1->Clone(n);

	double *x   = g3->GetX(); 
	double *xl  = g3->GetEXlow(); 
	double *xu  = g3->GetEXhigh(); 

	double *x1 = g2->GetY(); 
	double *x1l = ((TGraphAsymmErrors*)g2)->GetEYlow(); //
	double *x1u = ((TGraphAsymmErrors*)g2)->GetEYhigh();// 
	for(int i=0;i<N;i++){
		x[i] = x1[i];    xl[i]= x1l[i];    xu[i]= x1u[i];
	}
	return g3;
	//= new TGraphAsymmErrors(N,y2,y1,y2l,y2u,y1l,y1u); 
}
double shvals[]={0,-0.1,0.1,-0.2,0.2,-0.3,0.3,-0.4,0.4};
void shiftgr(TGraph*gr,int ipt){
	double sh  = shvals[ipt]/2;
	int N = gr->GetN();
	double *X = gr->GetX();
	double del;
	for(int i=0;i<N;i++){
		if(i<N-1) del= X[i+1]-X[i];
		else del= X[i]-X[i-1];
		X[i]+=sh*del;    
	}
}

enum{
	NBIN=5,
	NHAR=4,
	NPT=6,
	NPTS=4,
	NPTS1=2,
	NPTS2=2,
	NT2=7,//NT1 + converted 
	NT1=3,//ET, Nch, Cent
	NT=2,//ET, Nch
};

char *labss[]={"#Sigma E_{T}","N_{ch}^{rec}"};

char *otypes[]={"FCal","Nch","Cent","FCalb","Nchb"};
char *atypes[][6]={
	{"#Sigma E_{T} [TeV]","N_{ch}^{rec}", "1 - Centrality",  "#Sigma E_{T} [TeV]","N_{ch}^{rec}", "Centrality [%]"},
	{"#LT#Sigma E_{T}#GT/(4.1TeV)","#LTN_{ch}^{rec}#GT/2800", "1 - Centrality","#LTN_{ch}^{rec}#GT/2800","#LT#Sigma E_{T}#GT/(4.1TeV)","Centrality [%]"},
	//  {"#Sigma E_{T}/(4.1TeV)","N_{ch}^{rec}/2800", "1 - Centrality","#LTN_{ch}^{rec}#GT/2800","#LT#Sigma E_{T}#GT/(4.1TeV)","Centrality [%]"},
};


char *tnames[]={
	"1sub","3sub",
};
//TF1*fun1 = new TFile("fun1","[0]/x
char *ptname1[]={  "0.5<p_{T}<5 GeV",  "1.0<p_{T}<5 GeV",  "1.2<p_{T}<5 GeV",  "1.5<p_{T}<5 GeV",  "1.7<p_{T}<5 GeV",  "2.0<p_{T}<5 GeV"};
char *ptname2[][NHAR]={
	{"c_{1}{4, #SigmaE_{T}}", "c_{2}{4, #SigmaE_{T}}", "c_{3}{4, #SigmaE_{T}}", "c_{4}{4, #SigmaE_{T}}"},
	{"c_{1}{4, N_{ch}^{rec}}","c_{2}{4, N_{ch}^{rec}}","c_{3}{4, N_{ch}^{rec}}", "c_{4}{4, N_{ch}^{rec}}"},
};
char *ptname3[][NHAR]={
	{"c_{1}{6, #SigmaE_{T}}", "c_{2}{6, #SigmaE_{T}}", "c_{3}{6, #SigmaE_{T}}", "c_{4}{6, #SigmaE_{T}}"},
	{"c_{1}{6, N_{ch}^{rec}}","c_{2}{6, N_{ch}^{rec}}","c_{3}{6, N_{ch}^{rec}}","c_{4}{6, N_{ch}^{rec}}"},
};
char *scnames[NHAR]={
	"sc_{1,2}","sc_{1,3}","sc_{2,3}","sc_{2,4}"
};
char *acnames[NHAR]={
	"ac_{1,1,2}","ac_{x}","ac_{2}","ac_{y}"
};
char *testnames[]={"isGauss_1sub","isPower_1sub"};
//char *testnames[]={"isGauss","isPower"}; old

char *ctypes[]={"Standard method", "Three-subevent method"};
int ptbins[]={0,2,3,4,5,6};
int pts[]={0,1,3,5};
int pts1[]={3,5};
int pts2[]={0,5};
double mul1[]={1,3};
char name[200];

TGraph*hc2[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hc4[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hc6[NBIN+1][NHAR][NT2][NPT][2];

TGraph*hr4[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hr6[NBIN+1][NHAR][NT2][NPT][2];

TGraph*hv2[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hv4[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hv6[NBIN+1][NHAR][NT2][NPT][2];

TGraph*hrat1[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hrat2[NBIN+1][NHAR][NT2][NPT][2];

//3sub results
TGraph*hc4b[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hr4b[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hv4b[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hrat1b[NBIN+1][NHAR][NT2][NPT][2];

//sc and ac
TGraph*hsc[NBIN+1][NHAR][NT2][NPT][2];//sc(1,2),sc(1,3), sc(2,3), sc(2,4)
TGraph*hac[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hnsc[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hnac[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hscb[NBIN+1][NHAR][NT2][NPT][2];//sc(1,2),sc(1,3), sc(2,3), sc(2,4)
TGraph*hacb[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hnscb[NBIN+1][NHAR][NT2][NPT][2];
TGraph*hnacb[NBIN+1][NHAR][NT2][NPT][2];

TGraph*hrr[NBIN+1][NHAR][10][NPT][2];//ratio between the two types
TGraph*hcor[NBIN+1][10][2][NPT][2];//correlation

TGraph*g1[10][2],*g2[10][2];
TLatex lab1,lab2;
TLatex text,text1;
TLine line,line1,line2;  
double perc[][3]={
	{4.26/4.1,   4.08/4.1,  3.92/4.1},
	{2993./2800,2872./2800,2762./2800},
};

TH1*hcoora[NT1+3],*hcoorb[NT1+3];
//char *tname = otypes[type], *aname = atypes[type]; cout<<tname<<" "<<aname<<endl;
double lims[] ={5.1, 3800,  1.02, 5.1, 3600,  1.02}; // FIX
double lims0[]={0,   0,     0.2,  0,   0,     0};
double lims1[]={3.6, 2400,  0.2,  3.6, 2400,  80};
double lims2[]={4.9, 3800,  1.02, 4.9, 3600,  80};
double lims3[]={0.35,0.35,  0.7,  0.35,0.35,  30};
double lims4[]={0.6,  0.6,  0.85, 0.6,  0.6,  15};
TLegend*l;
TCanvas*c1;
int sty[]={20,24,25,27,28,30,31};
int col[] ={1,2,4,6,3,434,8};
int fcol[]={16,623,591,607,407,422,8};
int lsty[]={1,1,1,1,1,1,1,1,1,1,1,1};
int lsty2[]={2,2,2,2,2,2,2,2,2};
//for comparision between two curves
int bsty[]={20,24};
int bcol[] ={1,2};
int bfcol[]={16,623};

void draw1();
void drawcumu5sidebyside3(int type=1){//Cent or FCal
	lab1.SetNDC(1);
	lab2.SetNDC(1);
	lab1.SetTextSize(16);
	lab2.SetTextSize(15);
	lab1.SetTextFont(43);
	lab2.SetTextFont(43);

	if(type==0){
		ran1=0.4;ran2=5.1;
		Knee[0]=Knee[1]=Knee[2]=1;
	}else{
		ran1=0.4/Knee[0];ran2=5.1/Knee[0];
		for(int it1=0;it1<NT1+3;it1++){
			//lims[it1]/=Knee[it1%NT1];
			lims1[it1]/=Knee[it1%NT1];
			lims2[it1]/=Knee[it1%NT1];
		}
		lims[0]=1.25;    lims[1]=1.32; lims[2]=1.02; lims[3]=1.32;lims[4]=1.25; lims[5]=80;
	}
	//ET,NCh,Cent, ET mapped to Ntrk, NCh mapped to ET, Cent%
	for(int it1=0;it1<NT1+3;it1++){
		sprintf(name,"hcoora%d",it1);
		hcoora[it1] = new TH1F(name,name,600,lims0[it1],lims[it1]);
		hcoora[it1]->GetXaxis()->SetTitle(atypes[type][it1]);
		setstyle(hcoora[it1]);     hcoora[it1]->SetLineStyle(3);
		sprintf(name,"hcoorb%d",it1);
		hcoorb[it1] = new TH1F(name,name,600,lims1[it1],lims2[it1]);
		hcoorb[it1]->GetXaxis()->SetTitle(atypes[type][it1]);
		setstyle(hcoorb[it1]);     hcoorb[it1]->SetLineStyle(3);
	}
	text.SetNDC(1);  text.SetTextFont(43);  text.SetTextSize(16);
	text1.SetNDC(1);  text1.SetTextFont(43);  text1.SetTextSize(17);
	line.SetLineStyle(3);
	line1.SetLineWidth(2); line2.SetLineWidth(2);
	line1.SetLineStyle(7); line1.SetLineColor(4);
	line2.SetLineStyle(7); line2.SetLineColor(2);
	TFile*f,*f1,*f2;
	int pc=1;

	TGraphErrors*gmap[10];//mapping between different x-axes.
	f2 =new TFile("../data/paper/hist_cvt.root"); // FIX
	gmap[0]= (TGraphErrors*)f2->Get("gCvt_NchRec_FCal");
	gmap[1]= (TGraphErrors*)f2->Get("gCvt_FCal_NchRec");

	int aa=1;
	for(int ib=0;ib<NBIN;ib++){
		for(int it1=0;it1<NT1;it1++){
			if(ib<4){
				if(it1==0) sprintf(name,"../data/paper/hist_PbPb502_binFCal_bin%d.root",ib); 
				else if (it1==1) sprintf(name,"../data/paper/hist_PbPb502_binNch_bin%d.root",ib+1); 
				else sprintf(name,"../data/paper/hist_PbPb502_binCent_bin%d.root",ib+1); 
			}else{
				if(it1==0) sprintf(name,"../data/paper/hist_PbPb502_binFCal_bin%d.root",4); //corsest
				else if (it1==1) sprintf(name,"../data/paper/hist_PbPb502_binNch_bin%d.root",0); //finest
				else sprintf(name,"../data/paper/hist_PbPb502_binCent_bin%d.root",0); //fineset
			}
			cout<<name<<endl;      f1 = new TFile(name);
			for(int ihar=0;ihar<NHAR;ihar++){//v1-v4
				for(int it=0;it<NPT;it++){
					sprintf(name,"sts_c2_%s_Har%d_Pt%d",tnames[1],ihar+1,ptbins[it]);
					hc2[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_c2_%s_Har%d_Pt%d",tnames[1],ihar+1,ptbins[it]);
					hc2[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_c4_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hc4[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_c4_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hc4[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_c6_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hc6[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_c6_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hc6[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);

					sprintf(name,"sts_nc4_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hr4[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nc4_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hr4[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_nc6_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hr6[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nc6_%s_Har%d_Pt%d",tnames[0],ihar+1,ptbins[it]);
					hr6[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);

					//read in 3subevent results
					sprintf(name,"sts_c4_%s_Har%d_Pt%d",tnames[1],ihar+1,ptbins[it]);
					hc4b[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_c4_%s_Har%d_Pt%d",tnames[1],ihar+1,ptbins[it]);
					hc4b[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_nc4_%s_Har%d_Pt%d",tnames[1],ihar+1,ptbins[it]);
					hr4b[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nc4_%s_Har%d_Pt%d",tnames[1],ihar+1,ptbins[it]);
					hr4b[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					//v2
					sprintf(name,"%s_v2",hc2[ib][ihar][it1][it][0]->GetName());
					hv2[ib][ihar][it1][it][0] = (TGraph*)hc2[ib][ihar][it1][it][0]->Clone(name);
					sprintf(name,"%s_v2",hc2[ib][ihar][it1][it][1]->GetName());
					hv2[ib][ihar][it1][it][1] = (TGraph*)hc2[ib][ihar][it1][it][1]->Clone(name);
					grpow(hv2[ib][ihar][it1][it][0],1./2);
					grpow(hv2[ib][ihar][it1][it][1],1./2);
					//v4
					sprintf(name,"%s_v4",hc4[ib][ihar][it1][it][0]->GetName());
					hv4[ib][ihar][it1][it][0] = (TGraph*)hc4[ib][ihar][it1][it][0]->Clone(name);
					sprintf(name,"%s_v4",hc4[ib][ihar][it1][it][1]->GetName());
					hv4[ib][ihar][it1][it][1] = (TGraph*)hc4[ib][ihar][it1][it][1]->Clone(name);
					grpow(hv4[ib][ihar][it1][it][0],1./4); scalgr(hv4[ib][ihar][it1][it][0],-1);
					grpow(hv4[ib][ihar][it1][it][1],1./4); scalgr(hv4[ib][ihar][it1][it][1],-1);

					sprintf(name,"%s_v4",hc4b[ib][ihar][it1][it][0]->GetName());
					hv4b[ib][ihar][it1][it][0] = (TGraph*)hc4b[ib][ihar][it1][it][0]->Clone(name);
					sprintf(name,"%s_v4",hc4b[ib][ihar][it1][it][1]->GetName());
					hv4b[ib][ihar][it1][it][1] = (TGraph*)hc4b[ib][ihar][it1][it][1]->Clone(name);
					grpow(hv4b[ib][ihar][it1][it][0],1./4); scalgr(hv4b[ib][ihar][it1][it][0],-1);
					grpow(hv4b[ib][ihar][it1][it][1],1./4); scalgr(hv4b[ib][ihar][it1][it][1],-1);
					//v6
					sprintf(name,"%s_v6",hc6[ib][ihar][it1][it][0]->GetName());
					hv6[ib][ihar][it1][it][0] = (TGraph*)hc6[ib][ihar][it1][it][0]->Clone(name);
					sprintf(name,"%s_v6",hc6[ib][ihar][it1][it][1]->GetName());
					hv6[ib][ihar][it1][it][1] = (TGraph*)hc6[ib][ihar][it1][it][1]->Clone(name);
					grpow(hv6[ib][ihar][it1][it][0],1./6); scalgr(hv6[ib][ihar][it1][it][0],pow(4,-1./6));
					grpow(hv6[ib][ihar][it1][it][1],1./6); scalgr(hv6[ib][ihar][it1][it][1],pow(4,-1./6));

					//v{4}/v{2}
					sprintf(name,"%s_rat4",hr4[ib][ihar][it1][it][0]->GetName());
					hrat1[ib][ihar][it1][it][0] = (TGraph*)hr4[ib][ihar][it1][it][0]->Clone(name);
					sprintf(name,"%s_rat4",hr4[ib][ihar][it1][it][1]->GetName());
					hrat1[ib][ihar][it1][it][1] = (TGraph*)hr4[ib][ihar][it1][it][1]->Clone(name);
					grpow(hrat1[ib][ihar][it1][it][0],1./4); scalgr(hrat1[ib][ihar][it1][it][0],-1);
					grpow(hrat1[ib][ihar][it1][it][1],1./4); scalgr(hrat1[ib][ihar][it1][it][1],-1);

					sprintf(name,"%s_rat4",hr4b[ib][ihar][it1][it][0]->GetName());
					hrat1b[ib][ihar][it1][it][0] = (TGraph*)hr4b[ib][ihar][it1][it][0]->Clone(name);
					sprintf(name,"%s_rat4",hr4b[ib][ihar][it1][it][1]->GetName());
					hrat1b[ib][ihar][it1][it][1] = (TGraph*)hr4b[ib][ihar][it1][it][1]->Clone(name);
					grpow(hrat1b[ib][ihar][it1][it][0],1./4); scalgr(hrat1b[ib][ihar][it1][it][0],-1);
					grpow(hrat1b[ib][ihar][it1][it][1],1./4); scalgr(hrat1b[ib][ihar][it1][it][1],-1);

					//v{6}/v{4}
					sprintf(name,"sts_%s_Har%d_Pt%d",testnames[0],ihar+1,ptbins[it]);
					hrat2[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_%s_Har%d_Pt%d",testnames[0],ihar+1,ptbins[it]);
					hrat2[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					grpow(hrat2[ib][ihar][it1][it][0],1./6);
					grpow(hrat2[ib][ihar][it1][it][1],1./6);
					//sc and ac
					sprintf(name,"sts_sc_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hsc[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_sc_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hsc[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_nsc_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hnsc[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nsc_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hnsc[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_ac_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hac[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_ac_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hac[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_nac_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hnac[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nac_%s_Har%d_Pt%d",tnames[0],ihar,ptbins[it]);
					hnac[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);

					sprintf(name,"sts_sc_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hscb[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_sc_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hscb[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_nsc_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hnscb[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nsc_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hnscb[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_ac_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hacb[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_ac_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hacb[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
					sprintf(name,"sts_nac_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hnacb[ib][ihar][it1][it][0] =  toasygr((TGraphErrors*) (f1->Get(name)));
					sprintf(name,"sys_nac_%s_Har%d_Pt%d",tnames[1],ihar,ptbins[it]);
					hnacb[ib][ihar][it1][it][1] =  (TGraph*) f1->Get(name);
				}
			}
		}

		//mapping
		for(int it1=NT1;it1<NT1+3;it1++){
			for(int ihar=0;ihar<NHAR;ihar++){//v1-v4
				for(int it=0;it<NPT;it++){
					for(int ie=0;ie<2;ie++){
						TGraph*gtmp[30],*gtmp1[30];
						int itt=1; if(it1==3) itt=0; else if(it1==4) itt=1; else itt=2;
						gtmp1[0]=hc2[ib][ihar][itt][it][ie];
						gtmp1[1]=hc4[ib][ihar][itt][it][ie];
						gtmp1[2]=hc6[ib][ihar][itt][it][ie];
						gtmp1[3]=hr4[ib][ihar][itt][it][ie];
						gtmp1[4]=hr6[ib][ihar][itt][it][ie];
						gtmp1[5]=hv2[ib][ihar][itt][it][ie];
						gtmp1[6]=hv4[ib][ihar][itt][it][ie];
						gtmp1[7]=hv6[ib][ihar][itt][it][ie];
						gtmp1[8]=hrat1[ib][ihar][itt][it][ie];
						gtmp1[9]=hrat2[ib][ihar][itt][it][ie];
						gtmp1[10]=hc4b[ib][ihar][itt][it][ie];
						gtmp1[11]=hr4b[ib][ihar][itt][it][ie];
						gtmp1[12]=hv4b[ib][ihar][itt][it][ie];
						gtmp1[13]=hrat1b[ib][ihar][itt][it][ie];
						//sc&ac
						gtmp1[14]=hsc[ib][ihar][itt][it][ie];
						gtmp1[15]=hnsc[ib][ihar][itt][it][ie];
						gtmp1[16]=hac[ib][ihar][itt][it][ie];
						gtmp1[17]=hnac[ib][ihar][itt][it][ie];
						gtmp1[18]=hscb[ib][ihar][itt][it][ie];
						gtmp1[19]=hnscb[ib][ihar][itt][it][ie];
						gtmp1[20]=hacb[ib][ihar][itt][it][ie];
						gtmp1[21]=hnacb[ib][ihar][itt][it][ie];
						for(int iw=0;iw<22;iw++){//clone with different name
							sprintf(name,"%s_map",gtmp1[iw]->GetName()); gtmp[iw] = (TGraph*)gtmp1[iw]->Clone(name);
						}
						for(int iw=0;iw<22;iw++){//map axis
							if(it1==3){//ET to Ntrk
								grmap(gtmp[iw],gmap[1]);
							}else if(it1==4){//Ntrk to ET
								grmap(gtmp[iw],gmap[0]);
							}else{//Cent Cent
							}
						}
						hc2[ib][ihar][it1][it][ie]   =gtmp[0] ;
						hc4[ib][ihar][it1][it][ie]   =gtmp[1] ;
						hc6[ib][ihar][it1][it][ie]   =gtmp[2] ;
						hr4[ib][ihar][it1][it][ie]   =gtmp[3] ;
						hr6[ib][ihar][it1][it][ie]   =gtmp[4] ;
						hv2[ib][ihar][it1][it][ie]   =gtmp[5] ;
						hv4[ib][ihar][it1][it][ie]   =gtmp[6] ;
						hv6[ib][ihar][it1][it][ie]   =gtmp[7] ;
						hrat1[ib][ihar][it1][it][ie] =gtmp[8] ;
						hrat2[ib][ihar][it1][it][ie] =gtmp[9] ;
						hc4b[ib][ihar][it1][it][ie]  =gtmp[10];
						hr4b[ib][ihar][it1][it][ie]  =gtmp[11];
						hv4b[ib][ihar][it1][it][ie]  =gtmp[12];
						hrat1b[ib][ihar][it1][it][ie]=gtmp[13];

						hsc[ib][ihar][it1][it][ie]   =gtmp[14];
						hnsc[ib][ihar][it1][it][ie]  =gtmp[15];
						hac[ib][ihar][it1][it][ie]   =gtmp[16];
						hnac[ib][ihar][it1][it][ie]  =gtmp[17];
						hscb[ib][ihar][it1][it][ie]  =gtmp[18];
						hnscb[ib][ihar][it1][it][ie] =gtmp[19];
						hacb[ib][ihar][it1][it][ie]  =gtmp[20];
						hnacb[ib][ihar][it1][it][ie] =gtmp[21];

						if(ihar==2&&ie==0&&it==0&&it1==NT1){//Check ET to Nch
							cout<<ib<<endl;
							cout<<gtmp[4]->GetN()<<endl;
						}
					}
				}
			}
		}
		//	  if(type||it1==2){//always do 1-centrality 
		for(int it1=0;it1<NT1+3;it1++){//rescal the axis
			for(int ihar=0;ihar<NHAR;ihar++){//v1-v4
				for(int it=0;it<NPT;it++){
					for(int ie=0;ie<2;ie++){
						grcent(hc2[ib][ihar][it1][it][ie],it1); 
						grcent(hc4[ib][ihar][it1][it][ie],it1);
						grcent(hc6[ib][ihar][it1][it][ie],it1);
						grcent(hr4[ib][ihar][it1][it][ie],it1);
						grcent(hr6[ib][ihar][it1][it][ie],it1);
						grcent(hv2[ib][ihar][it1][it][ie],it1);
						grcent(hv4[ib][ihar][it1][it][ie],it1);
						grcent(hv6[ib][ihar][it1][it][ie],it1);
						grcent(hrat1[ib][ihar][it1][it][ie],it1);
						grcent(hrat2[ib][ihar][it1][it][ie],it1);
						grcent(hc4b[ib][ihar][it1][it][ie],it1);
						grcent(hr4b[ib][ihar][it1][it][ie],it1);
						grcent(hv4b[ib][ihar][it1][it][ie],it1);
						grcent(hrat1b[ib][ihar][it1][it][ie],it1);
						//sc & ac
						grcent(hsc[ib][ihar][it1][it][ie],it1);
						grcent(hnsc[ib][ihar][it1][it][ie],it1);
						grcent(hac[ib][ihar][it1][it][ie],it1);
						grcent(hnac[ib][ihar][it1][it][ie],it1);
						grcent(hscb[ib][ihar][it1][it][ie],it1);
						grcent(hnscb[ib][ihar][it1][it][ie],it1);
						grcent(hacb[ib][ihar][it1][it][ie],it1);
						grcent(hnacb[ib][ihar][it1][it][ie],it1);
					}
				}
			}
		}
	}
	cout<<"hi"<<endl;

	/*
		 for(int ib=0;ib<NBIN;ib++){
		 for(int ihar=0;ihar<NHAR;ihar++){//v1-v4
		 for(int it=0;it<NPT;it++){
		 int it1=1;	    int NP=0;
		 if(ihar==1&&ib==0) NP=5;
		 int N = hc4[ib][ihar][it1][it][0]->GetN();
		 for(int ibin=0;ibin<NP;ibin++){
		 for(int ie=0;ie<2;ie++){
		 hc2[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hc4[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hc6[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hr4[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hr6[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv2[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv4[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv6[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hrat1[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hrat2[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hc4b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hr4b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv4b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hrat1b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);

		 hsc[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnsc[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hac[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnac[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hscb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnscb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hacb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnacb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 }
		 }
		 N-=NP;
		 if(ib>0){
		 for(int ibin=0;ibin<NP;ibin++){
		 for(int ie=0;ie<2;ie++){
		 hc2[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hc4[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hc6[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hr4[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hr6[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv2[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv4[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv6[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hrat1[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hrat2[ib][ihar][it1][it][ie]->RemovePoint(N-NP);

		 hc4b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hr4b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hv4b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hrat1b[ib][ihar][it1][it][ie]->RemovePoint(N-NP);

		 hsc[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnsc[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hac[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnac[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hscb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnscb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hacb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 hnacb[ib][ihar][it1][it][ie]->RemovePoint(N-NP);
		 }
		 }
		 N-=NP;
		 }
		 }
		 }
		 }
		 */

	//CMS vn6/vn4 results:
	TFile *fcms = new TFile("../data/paper/cms.root");
	TGraph*grcms = (TGraph*) fcms->Get("CMS_ratio_v6_v4");  grcent(grcms,2); grcms->RemovePoint(0);

	grcms->SetMarkerSize(1);  grcms->SetFillColor(590); grcms->SetLineColor(590); grcms->SetMarkerStyle(25);  grcms->SetMarkerColor(4);
	TGraph* grcms1 = (TGraph*) grcms->Clone();
	grcms1->SetMarkerSize(1);  grcms1->SetMarkerColor(4);

	grcms->Print();
	cout<<"done1"<<endl;;
	for(int ib=0;ib<NBIN;ib++){
		for(int ihar=0;ihar<NHAR;ihar++){
			for(int it1=0;it1<NT1+3;it1++){
				for(int it=0;it<NPT;it++){
					for(int ie=0;ie<2;ie++){
						TGraph*gtmp[30];
						gtmp[0]=hc2[ib][ihar][it1][it][ie];
						gtmp[1]=hc4[ib][ihar][it1][it][ie];
						gtmp[2]=hc6[ib][ihar][it1][it][ie];
						gtmp[3]=hr4[ib][ihar][it1][it][ie];
						gtmp[4]=hr6[ib][ihar][it1][it][ie];
						gtmp[5]=hv2[ib][ihar][it1][it][ie];
						gtmp[6]=hv4[ib][ihar][it1][it][ie];
						gtmp[7]=hv6[ib][ihar][it1][it][ie];
						gtmp[8]=hrat1[ib][ihar][it1][it][ie];
						gtmp[9]=hrat2[ib][ihar][it1][it][ie];
						gtmp[10]=hc4b[ib][ihar][it1][it][ie];
						gtmp[11]=hr4b[ib][ihar][it1][it][ie];
						gtmp[12]=hv4b[ib][ihar][it1][it][ie];
						gtmp[13]=hrat1b[ib][ihar][it1][it][ie];
						//sc&ac
						gtmp[14]=hsc[ib][ihar][it1][it][ie];
						gtmp[15]=hnsc[ib][ihar][it1][it][ie];
						gtmp[16]=hac[ib][ihar][it1][it][ie];
						gtmp[17]=hnac[ib][ihar][it1][it][ie];
						gtmp[18]=hscb[ib][ihar][it1][it][ie];
						gtmp[19]=hnscb[ib][ihar][it1][it][ie];
						gtmp[20]=hacb[ib][ihar][it1][it][ie];
						gtmp[21]=hnacb[ib][ihar][it1][it][ie];

						for(int iw=0;iw<22;iw++){
							if(ie==0){
								gtmp[iw]->SetMarkerStyle(sty[it]); gtmp[iw]->SetMarkerColor(col[it]); gtmp[iw]->SetLineColor(col[it]);
							}else{
								gtmp[iw]->SetMarkerSize(0); gtmp[iw]->SetMarkerColor(fcol[it]); gtmp[iw]->SetLineColor(fcol[it]);
								gtmp[iw]->SetLineWidth(5);
							}
						}
					}
				}
			}
		}
	}
	//create cumulant correlator
	for(int ib=0;ib<NBIN;ib++){
		for(int it1=0;it1<2;it1++){
			for(int it=0;it<NPT;it++){
				for(int ie=0;ie<2;ie++){
					hcor[ib][0][it1][it][ie] = grcorr(hr6[ib][1][it1][it][ie],  hr4[ib][1][it1][it][ie]);//nc6 vs nc4
					hcor[ib][1][it1][it][ie] = grcorr(hr6[ib][2][it1][it][ie],  hr4[ib][2][it1][it][ie]);//
					hcor[ib][2][it1][it][ie] = grcorr(hr6[ib][3][it1][it][ie],  hr4[ib][3][it1][it][ie]);//
					hcor[ib][3][it1][it][ie] = grcorr(hrat2[ib][1][it1][it][ie],hrat1[ib][1][it1][it][ie]);//v2 vs v2
					hcor[ib][4][it1][it][ie] = grcorr(hrat2[ib][2][it1][it][ie],hrat1[ib][2][it1][it][ie]);//v3 vs v3
					hcor[ib][5][it1][it][ie] = grcorr(hrat2[ib][3][it1][it][ie],hrat1[ib][3][it1][it][ie]);//v4 vs v4

					hcor[ib][6][it1][it][ie] = grcorr(hrat1[ib][2][it1][it][ie],hrat1[ib][1][it1][it][ie]);//v3 vs v2
					hcor[ib][7][it1][it][ie] = grcorr(hrat1[ib][3][it1][it][ie],hrat1[ib][1][it1][it][ie]);//v4 vs v2
				}
			}
		}
	}
	//ratio
	for(int ib=0;ib<NBIN;ib++){
		for(int ihar=0;ihar<NHAR;ihar++){//v1-v4
			//cout<<ihar<<endl;
			for(int it=0;it<NPT;it++){
				//cn{2} ratios  Nch/Et
				hrr[ib][ihar][0][it][0] = grratio(hv2[ib][ihar][4][it][0],hv2[ib][ihar][0][it][0],0);
				hrr[ib][ihar][0][it][1] = grratio(hv2[ib][ihar][4][it][1],hv2[ib][ihar][0][it][1],0);
				setgrerror(hrr[ib][ihar][0][it][1],0.005);
				hrr[ib][ihar][1][it][0] = grratio(hv2[ib][ihar][3][it][0],hv2[ib][ihar][1][it][0],1);
				hrr[ib][ihar][1][it][1] = grratio(hv2[ib][ihar][3][it][1],hv2[ib][ihar][1][it][1],1);
				setgrerror(hrr[ib][ihar][1][it][1],0.005);
				//cn{4} ratios Nch/Et
				hrr[ib][ihar][2][it][0] = grratio(hv4[ib][ihar][4][it][0],hv4[ib][ihar][0][it][0],0);
				hrr[ib][ihar][2][it][1] = grratio(hv4[ib][ihar][4][it][1],hv4[ib][ihar][0][it][1],0);
			}
		}
	}
	sprintf(name,"5TeVPbPbsbys3.root");  f=new TFile(name,"recreate");
	cout<<"done2"<<endl;
	//draw1();
	cout<<"next"<<endl;
	//paper 
	int pads1[]={0,4,2,1,3};
	double mmax1[4]={0},mmin1[4]={1};
	for(int itt=0;itt<5;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		sprintf(name,"comp_vn2_%s",otypes[it1]);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
		for(int ihar=1;ihar<NHAR;ihar++){
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			TH1*htmp[6];
			double mmax=0,mmin=1;
			c1->cd(ihar);      gPad->SetTopMargin(0.01);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][0] =  (TGraph*)hv2[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hv2[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.1; if(ihar==2) th=0.05;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				if(itt==0) { ran1=0.2/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=200/Knee[1];ran2=3000/Knee[1];}else {ran1=0.1;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			if(itt==0){//backup
				mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
			}
			gPad->SetTicks(1,1);
			htmp[itt] = (TH1*)hcoora[it1]->DrawClone();
			if(ihar==1)  htmp[itt]->SetMaximum(mmax1[ihar]*1.35); else  htmp[itt]->SetMaximum(mmax1[ihar]*1.4);      htmp[itt]->SetMinimum(0);
			if(it1==1||it1==4)  sprintf(name,"v_{%d}{2,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{2,#Sigma E_{T}}",ihar+1);

			htmp[itt]->GetYaxis()->SetTitle(name);      htmp[itt]->GetXaxis()->SetTitleOffset(1.);

			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}
			line.DrawLine(perc[0][0],0,perc[0][0],0.1);

			l = new TLegend(0.6, 0.72, 1, 0.97);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			double shif=0.15,shix=0.07;
			if(ihar==1){
				l->Draw();
			}
			// text1.DrawLatex(0.2,0.9,"Two-subevent method");
			lab1.DrawLatex(0.2+shix,0.1+shif,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.2+shix,0.03+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
		}
		sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
	}
	//vn2 ratios
	for(int it1=0;it1<2;it1++){//Nch/Et, Et/Nch
		sprintf(name,"comp_vn2rat_%s",otypes[it1]);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
		int ib=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			if(ihar==1) ib=1; if(ihar==2) ib=1; if(ihar==3) ib=2;
			TH1*htmp[6];
			double mmax=0,mmin=1;
			c1->cd(ihar);      gPad->SetTopMargin(0.01);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS2;it++){
				g1[it][0] =  (TGraph*)hrr[ib][ihar][it1][pts2[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hrr[ib][ihar][it1][pts2[it]][1]->Clone();

				TGraph*gtmp=g1[it][0];
				gtmp->SetMarkerStyle(sty[it]); gtmp->SetMarkerColor(col[it]); gtmp->SetLineColor(col[it]);	
				gtmp=g1[it][1];
				gtmp->SetMarkerSize(0); gtmp->SetMarkerColor(fcol[it]); gtmp->SetLineColor(fcol[it]);
				gtmp->SetLineWidth(5);

				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.02; if(ihar==2) th=0.02;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				if(it1==0) { ran1=0.2/Knee[0];ran2=4.4/Knee[0]; } else {ran1=200/Knee[1];ran2=3000/Knee[1];}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			if(it1==0)  htmp[ihar] = (TH1*)hcoora[0]->DrawClone();  else  htmp[ihar] = (TH1*)hcoora[1]->DrawClone();
			if(it1==0){
				htmp[ihar]->SetMaximum(1.089);      htmp[ihar]->SetMinimum(0.96);
			}else{
				htmp[ihar]->SetMaximum(1.089);      htmp[ihar]->SetMinimum(0.96);
			}
			htmp[ihar]->GetXaxis()->SetLimits(0.03,lims[it1]);
			line.DrawLine(0.03,1,lims[it1],1);

			sprintf(name,"v_{%d}{2,N_{ch}^{rec}}/v_{%d}{2,#Sigma E_{T}}",ihar+1,ihar+1); 

			htmp[ihar]->GetYaxis()->SetTitle(name);      htmp[ihar]->GetXaxis()->SetTitleOffset(1.);
			gPad->SetLeftMargin(0.18); htmp[ihar]->GetYaxis()->SetTitleOffset(1.8);

			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS2;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS2;it++){
				g1[it][0]->Draw("P");
			}
			//l = new TLegend(0.2, 0.82, 0.6, 0.97);
			l = new TLegend(0.2, 0.72, 0.6, 0.97);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS2;it++){
				l->AddEntry(g1[it][0],ptname1[pts2[it]],"pl");  
			}
			double shif=0.15,shix=0.07;
			if(ihar==1){
				l->Draw();
			}
			lab1.DrawLatex(0.2+shix,0.1+shif,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.2+shix,0.03+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
		}
		sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
	}
	//nc{4}
	for(int itt=0;itt<5;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		sprintf(name,"comp_nc4_%s",otypes[it1]);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
		for(int ihar=1;ihar<NHAR;ihar++){
			if(ihar==1){ if(it1%3==0) ib=0; else ib=0;}
			if(ihar==2){ if(it1%3==0) ib=2; else ib=1;}
			if(ihar==3){ if(it1%3==0) ib=2; else ib=1;}
			TH1*htmp[6];
			double mmax=0,mmin=1;
			c1->cd(ihar);      gPad->SetTopMargin(0.01);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][0] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.05; if(ihar==2) th=0.05;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			if(itt==0){//backup
				mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
			}
			gPad->SetTicks(1,1);
			htmp[itt] = (TH1*)hcoora[it1]->DrawClone();
			if(ihar==1){ htmp[itt]->SetMaximum(mmax1[ihar]*2); htmp[itt]->SetMinimum(mmin1[ihar]*1.32);}
			if(ihar==2){ htmp[itt]->SetMaximum(mmax1[ihar]*6); htmp[itt]->SetMinimum(mmin1[ihar]*1.);}
			if(ihar==3){ htmp[itt]->SetMaximum(mmax1[ihar]*1.2); htmp[itt]->SetMinimum(mmin1[ihar]*1.4);}

			if(it1==1||it1==4)  sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+1);
			htmp[itt]->GetYaxis()->SetTitle(name);      htmp[itt]->GetXaxis()->SetTitleOffset(1.);
			if(itt==2){
				if(ihar==1) htmp[itt]->GetXaxis()->SetLimits(0.2,lims[itt]);
				if(ihar==2) htmp[itt]->GetXaxis()->SetLimits(0.4,lims[itt]);
				if(ihar==3) htmp[itt]->GetXaxis()->SetLimits(0.5,lims[itt]);
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}      
			double shiy=0,shix=0.;
			if(ihar==3) {shix=0.4;shiy=-0.2;}
			if(ihar==2) {shix=0.4;shiy=-0.4;}
			if(ihar==1) shiy-=0.11;
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			if(ihar==1){
				l->Draw();
			}
			shiy=shix=0;
			if(ihar==3) {shix=0.1;shiy=0.56;}
			if(ihar==2) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.28,0.92,ctypes[0]);
			lab1.DrawLatex(0.20+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.20+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
		}
		sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		//subevent
		sprintf(name,"comp_nc4sub_%s",otypes[it1]);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
		for(int ihar=1;ihar<NHAR;ihar++){
			if(ihar==1){ if(it1%3==0) ib=0; else ib=0;}
			if(ihar==3){ if(it1%3==2) ib=2; else ib=3;}
			TH1*htmp[6];
			double mmax=0,mmin=1;
			c1->cd(ihar);      gPad->SetTopMargin(0.01);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				if(ihar==2){ 
					if(it<2){if(it1%3==2) ib=2; else ib=2;}
					else {if(it1%3==2) ib=2; else ib=2;}
				}
				g1[it][0] =  (TGraph*)hr4b[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hr4b[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.05; if(ihar==2) th=0.05;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
			}
			gPad->SetTicks(1,1);
			htmp[itt] = (TH1*)hcoora[it1]->DrawClone();
			if(ihar==1){ htmp[itt]->SetMaximum(mmax1[ihar]*2); htmp[itt]->SetMinimum(mmin1[ihar]*1.32);}
			if(ihar==2){ htmp[itt]->SetMaximum(mmax1[ihar]*6); htmp[itt]->SetMinimum(mmin1[ihar]*1.);}
			if(ihar==3){ htmp[itt]->SetMaximum(mmax1[ihar]*1.2); htmp[itt]->SetMinimum(mmin1[ihar]*1.4);}

			if(it1==1||it1==4)  sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+1);
			htmp[itt]->GetYaxis()->SetTitle(name);      htmp[itt]->GetXaxis()->SetTitleOffset(1.);
			if(itt==2){
				if(ihar==1) htmp[itt]->GetXaxis()->SetLimits(0.2,lims[itt]);
				if(ihar==2) htmp[itt]->GetXaxis()->SetLimits(0.4,lims[itt]);
				if(ihar==3) htmp[itt]->GetXaxis()->SetLimits(0.5,lims[itt]);
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}      
			double shiy=0,shix=0.;
			if(ihar==3) {shix=0.4;shiy=-0.2;}
			if(ihar==2) {shix=0.4;shiy=-0.4;}
			if(ihar==1) shiy-=0.11;
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			if(ihar==1){
				l->Draw();
			}
			shiy=shix=0;
			if(ihar==3) {shix=0.1;shiy=0.56;}
			if(ihar==2) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.28,0.92,ctypes[1]);
			lab1.DrawLatex(0.20+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.20+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
		}
		sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
	}

	//vn{4}
	for(int itt=0;itt<5;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			sprintf(name,"comp_vn4_har%d_%s",ihar+1,otypes[it1]);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			if(ihar==1){ if(it1%3==0) ib=0; else ib=0;}
			if(ihar==2){ if(it1%3==0) ib=2; else ib=1;}
			if(ihar==3){ if(it1%3==0) ib=2; else ib=1;}
			TH1*htmp[6];
			double mmax=0,mmin=1;     gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				g1[it][0] =  (TGraph*)hv4[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hv4[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.03; if(ihar==2) th=0.03;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				N=g1[it][0]->GetN(); ggg = (TGraphAsymmErrors*)g1[it][1];
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			htmp[itt] = (TH1*)hcoora[it1]->DrawClone();
			if(ihar==1){ mmax=0.25;mmin=-0.1;}
			if(ihar==2){ mmax=0.12;mmin=-0.1;}
			if(ihar==3){ mmax=0.09;mmin=-0.1;}
			htmp[itt]->SetMaximum(mmax); htmp[itt]->SetMinimum(mmin);	      
			if(it1==1||it1==4)  sprintf(name,"v_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{4,#Sigma E_{T}}",ihar+1);
			htmp[itt]->GetYaxis()->SetTitle(name);      htmp[itt]->GetXaxis()->SetTitleOffset(1.);
			if(itt==2){
				if(ihar==1) htmp[itt]->GetXaxis()->SetLimits(0.2,lims[itt]);
				if(ihar==2) htmp[itt]->GetXaxis()->SetLimits(0.4,lims[itt]);
				if(ihar==3) htmp[itt]->GetXaxis()->SetLimits(0.5,lims[itt]);
			}

			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}
			l = new TLegend(0.6, 0.72, 1, 0.97);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			double shif=0.15,shix=0.07;
			if(it1==0){
				l->Draw();
			}
			text1.DrawLatex(0.2,0.9,ctypes[0]);
			lab1.DrawLatex(0.2+shix,0.09+shif,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.2+shix,0.02+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}

	//nc{6}
	for(int itt=0;itt<5;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			sprintf(name,"comp_nc6_har%d_%s",ihar+1,otypes[it1]);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			if(ihar==1){ if(it1%3<2) ib=1; else ib=0;}
			if(ihar==2){ if(it1%3==0) ib=3; else ib=3;}
			if(ihar==3){ if(it1%3==0) ib=3; else ib=3;}
			TH1*htmp[6];
			double mmax=0,mmin=1;      //      c1->cd(ihar);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);
			gPad->SetTicks(1,1);
			gPad->SetTopMargin(0.01);
			int NN=NPTS; if(ihar>1) NN=2;
			for(int it=0;it<NN;it++){
				g1[it][0] =  (TGraph*)hr6[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hr6[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.09; if(ihar>1) th=0.09;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				th=0.06; if(ihar==2) th=0.02; if(ihar==3) th=0.02;
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				N=g1[it][0]->GetN(); ggg = (TGraphAsymmErrors*)g1[it][1];
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			htmp[itt] = (TH1*)hcoora[it1]->DrawClone();
			if(it1==1||it1==4)  sprintf(name,"nc_{%d}{6,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{6,#Sigma E_{T}}",ihar+1);
			htmp[itt]->GetYaxis()->SetTitle(name);      htmp[itt]->GetXaxis()->SetTitleOffset(1.);
			if(itt==2){
				if(ihar==1) htmp[itt]->GetXaxis()->SetLimits(0.2,lims[itt]);
				if(ihar==2) htmp[itt]->GetXaxis()->SetLimits(0.4,lims[itt]);
				if(ihar==3) htmp[itt]->GetXaxis()->SetLimits(0.5,lims[itt]);
			}
			if(ihar==1) {mmax=0.69; mmin=-0.16;}
			if(ihar==2) {mmax=0.018; mmin=-0.02;}
			if(ihar==3) {mmax=0.06; mmin=-0.04;}
			htmp[itt]->SetMaximum(mmax);  htmp[itt]->SetMinimum(mmin);	
			gPad->SetTicks(1,1);
			for(int it=0;it<NN;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NN;it++){
				g1[it][0]->Draw("P");
			}

			double shiy=0,shix=0.;
			if(ihar==1){if(itt==2) shiy=-0.05; else {shix=0.4;shiy=0;}}
			if(ihar==2) {shix=0.;shiy=-0.07;}
			if(ihar==3) {shix=0;shiy=-0.07;}
			if(ihar==1) l = new TLegend(0.18+shix, 0.72+shiy, 0.5+shix, 0.97+shiy);
			else l = new TLegend(0.2+shix, 0.8+shiy, 0.55+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NN;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();
			shiy=shix=0;
			if(ihar>1) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.28,0.92,ctypes[0]);
			lab1.DrawLatex(0.20+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.20+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}
	//vn{6}
	for(int itt=0;itt<5;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		for(int ihar=1;ihar<2;ihar++){
			sprintf(name,"comp_vn6_har%d_%s",ihar+1,otypes[it1]);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			if(ihar==1){ if(it1%3==0) ib=0; else ib=0;}
			if(ihar==2){ if(it1%3==0) ib=2; else ib=1;}
			if(ihar==3){ if(it1%3==0) ib=2; else ib=1;}
			TH1*htmp[6];
			double mmax=0,mmin=1;     gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				g1[it][0] =  (TGraph*)hv6[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hv6[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.03; if(ihar==2) th=0.03;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				N=g1[it][0]->GetN(); ggg = (TGraphAsymmErrors*)g1[it][1];
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			htmp[itt] = (TH1*)hcoora[it1]->DrawClone();
			if(ihar==1){ mmax=0.25;mmin=-0.05;}
			if(ihar==2){ mmax=0.12;mmin=-0.1;}
			if(ihar==3){ mmax=0.09;mmin=-0.1;}
			htmp[itt]->SetMaximum(mmax); htmp[itt]->SetMinimum(mmin);	      
			if(it1==1||it1==4)  sprintf(name,"v_{%d}{6,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{6,#Sigma E_{T}}",ihar+1);
			htmp[itt]->GetYaxis()->SetTitle(name);      htmp[itt]->GetXaxis()->SetTitleOffset(1.);
			if(itt==2){
				if(ihar==1) htmp[itt]->GetXaxis()->SetLimits(0.2,lims[itt]);
				if(ihar==2) htmp[itt]->GetXaxis()->SetLimits(0.4,lims[itt]);
				if(ihar==3) htmp[itt]->GetXaxis()->SetLimits(0.5,lims[itt]);
			}

			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}
			l = new TLegend(0.6, 0.72, 1, 0.97);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			double shif=0.15,shix=0.07;
			if(it1==0){
				l->Draw();
			}
			text1.DrawLatex(0.2,0.9,ctypes[0]);
			lab1.DrawLatex(0.2+shix,0.09+shif,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.2+shix,0.02+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}


	//vn4/vn2
	for(int itt=0;itt<3;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			sprintf(name,"comp_vn4vn2_har%d_%s",ihar+1,otypes[it1]);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			if(ihar==1){ if(it1%3==0) ib=0; else ib=0;}
			if(ihar==2){ if(it1%3==0) ib=2; else ib=1;}
			if(ihar==3){ if(it1%3==0) ib=2; else ib=1;}
			double mmax=0,mmin=1;
			gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				g1[it][0] =  (TGraph*)hrat1[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hrat1[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.05; if(ihar==2) th=0.05;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				th=0.08; if(ihar==2) th=0.1; if(ihar==3) th=0.2;
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				ran1=0.3; ran2=0.9;
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			TH1*hh= (TH1*)hcoora[it1]->DrawClone();
			if(ihar==1){ mmax=1;mmin=0.7;}
			if(ihar==2){ mmax=1;mmin=0;}
			if(ihar==3){ mmax=1;mmin=-1.2;}
			hh->SetMaximum(mmax); hh->SetMinimum(mmin);	
			if(itt==1)  sprintf(name,"v_{%d}{4,N_{ch}^{rec}}/v_{%d}{2,N_{ch}^{rec}}",ihar+1,ihar+1); 
			else sprintf(name,"v_{%d}{4,#Sigma E_{T}}/v_{%d}{2,#Sigma E_{T}}",ihar+1,ihar+1); 
			hh->GetYaxis()->SetTitle(name);      hh->GetXaxis()->SetTitleOffset(1.);
			if(itt==2) {hh->GetXaxis()->SetLimits(0.2,lims[it1]);}
			else {hh->GetXaxis()->SetLimits(0,1);}
			if(itt==2){
				if(ihar==2) hh->GetXaxis()->SetLimits(0.4,lims[it1]);
				if(ihar==3) hh->GetXaxis()->SetLimits(0.4,lims[it1]);
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}

			double shiy=0,shix=0.4;
			if(ihar==3) {shiy=-0.4; shix=0.45;}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();
			shix=0;
			if(ihar==3) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.2,0.92,ctypes[0]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}

	//vn6/vn4

	for(int itt=0;itt<3;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		for(int ihar=1;ihar<2;ihar++){
			sprintf(name,"comp_vn6vn4_har%d_%s",ihar+1,otypes[it1]);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				if(ihar==1){
					if(itt==0) {ib=1; if(it>2) ib=2;}
					if(itt==1) {ib=1; if(it>2) ib=2;}
					if(itt==2) {ib=1; if(it>2) ib=2;}
				}
				g1[it][0] =  (TGraph*)hrat2[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hrat2[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.005;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-7)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-7)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				th=0.005;
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				ran1=0.3; ran2=0.9;
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			TH1*hh= (TH1*)hcoora[it1]->DrawClone();
			//hh->GetXaxis()->SetLimits(0.02,lims[itt]);
			if(itt==2) hh->GetXaxis()->SetLimits(0.3,1);
			mmax=1.012;mmin=0.98; if(itt==2) mmin=0.964;
			hh->SetMaximum(mmax); hh->SetMinimum(mmin);	
			if(itt==1)  sprintf(name,"v_{%d}{6,N_{ch}^{rec}}/v_{%d}{4,N_{ch}^{rec}}",ihar+1,ihar+1); 
			else sprintf(name,"v_{%d}{6,#Sigma E_{T}}/v_{%d}{4,#Sigma E_{T}}",ihar+1,ihar+1); 
			hh->GetYaxis()->SetTitle(name);      hh->GetXaxis()->SetTitleOffset(1.);

			if(itt==2) {hh->GetXaxis()->SetLimits(0.3,lims[it1]); line.DrawLine(0.3,1,lims[it1],1);}
			else {hh->GetXaxis()->SetLimits(0,1); line.DrawLine(0,1,1,1);}
			if(itt==2) {
				grcms->Draw("3");	grcms1->Draw("P");
				l = new TLegend(0.4, 0.4, 0.6, 0.45);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(13);
				l->AddEntry(grcms,"CMS unfolding 0.3<p_{T}<3 GeV");
				l->Draw();
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}
			double shiy=0,shix=0;
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();
			shix=0.2;shiy=0;

			text1.DrawLatex(0.6,0.92,ctypes[0]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}

	//nsc & nac
	for(int itt=0;itt<5;itt++){//Et, Nch and Cent.
		int it1 = pads1[itt];    int ib=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			if(ihar==1)      sprintf(name,"comp_nac_har%d_%s",ihar+1,otypes[it1]); 
			else      sprintf(name,"comp_nsc_har%d_%s",ihar,otypes[it1]); 
			c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			TH1*htmp[6];
			if(ihar==1){ if(it1==0) ib=1; else ib=1;}
			if(ihar==2){ if(it1==0) ib=1; else ib=1;}
			if(ihar==3){ if(it1==0) ib=1; else ib=1;}
			gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				if(ihar==1){
					g1[it][0] =  (TGraph*)hnac[ib][2][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnac[ib][2][it1][pts[it]][1]->Clone();
				}else{
					g1[it][0] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][1]->Clone();
				}
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];	
				double th=0.05; if(ihar==3) th=0.05;
				if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				th=0.1; if(ihar==3) th=0.15;
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
			if(itt==0){//backup
				mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
			}
			if(ihar==1) {htmp[it1]->SetMaximum(1.5*mmax1[ihar]);  htmp[it1]->SetMinimum(-0.05);}
			if(ihar==2){
				if(it1==2){htmp[it1]->SetMaximum(1.2*mmax1[ihar]);}else {htmp[it1]->SetMaximum(2.2*mmax1[ihar]);}	htmp[it1]->SetMinimum(1.5*mmin1[ihar]);
			}
			if(ihar==3) {htmp[it1]->SetMaximum(1.5*mmax1[ihar]);  htmp[it1]->SetMinimum(-0.05);}
			if(ihar==1){
				if(it1==1)  sprintf(name,"n%s{3,N_{ch}^{rec}}",acnames[2]); else  sprintf(name,"n%s{3,#Sigma E_{T}}",acnames[2]);
			}else{
				if(it1==1)  sprintf(name,"n%s{4,N_{ch}^{rec}}",scnames[ihar]); else  sprintf(name,"n%s{4,#Sigma E_{T}}",scnames[ihar]);
			}

			htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
			if(it1==2){
				if(ihar==1) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
				if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
				if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}

			double shiy=0,shix=0.;
			if(ihar==3) {
				if(it1!=2){shix=0.4;shiy=-0.35;}
				else {shix=0.;shiy=-0.54;}
			}
			if(ihar==2) {
				if(it1!=2){shix=0.4;shiy=-0.4;}
				else {shix=0.4;shiy=-0.4;}
			}
			if(ihar==1) {
				shix=0.4;shiy=0;
			}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();

			if(ihar==3) {
				if(it1!=2){shix=0.1;shiy=0.57;}
				else {shix=0.2;shiy=0.57;}
			}
			if(ihar==2) {shix=0.2;shiy=0;}
			if(ihar==1) {shix=0;shiy=0;}
			if(ihar==3&&it1==2)      text1.DrawLatex(0.35,0.9,ctypes[0]);      else text1.DrawLatex(0.3,0.9,ctypes[0]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//subevent
		for(int ihar=1;ihar<NHAR;ihar++){
			if(ihar==1)      sprintf(name,"comp_nacsub_har%d_%s",ihar+1,otypes[it1]); 
			else      sprintf(name,"comp_nscsub_har%d_%s",ihar,otypes[it1]); 
			c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			TH1*htmp[6];
			if(ihar==1){ if(it1==0) ib=1; else ib=1;}
			if(ihar==2){ if(it1==0) ib=1; else ib=1;}
			if(ihar==3){ if(it1==0) ib=1; else ib=1;}
			gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				if(ihar==1){
					g1[it][0] =  (TGraph*)hnacb[ib][2][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnacb[ib][2][it1][pts[it]][1]->Clone();
				}else{
					g1[it][0] =  (TGraph*)hnscb[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnscb[ib][ihar][it1][pts[it]][1]->Clone();
				}
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];	
				double th=0.05; if(ihar==3) th=0.05;
				if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				th=0.1; if(ihar==3) th=0.15;
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
			/*     if(itt==0){//backup
						 mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
						 }*/
			if(ihar==1) {htmp[it1]->SetMaximum(1.5*mmax1[ihar]);  htmp[it1]->SetMinimum(-0.05);}
			if(ihar==2){
				if(it1==2){htmp[it1]->SetMaximum(1.2*mmax1[ihar]);}else {htmp[it1]->SetMaximum(2.2*mmax1[ihar]);}	htmp[it1]->SetMinimum(1.5*mmin1[ihar]);
			}
			if(ihar==3) {htmp[it1]->SetMaximum(1.5*mmax1[ihar]);  htmp[it1]->SetMinimum(-0.05);}
			if(ihar==1){
				if(it1==1)  sprintf(name,"n%s{3,N_{ch}^{rec}}",acnames[2]); else  sprintf(name,"n%s{3,#Sigma E_{T}}",acnames[2]);
			}else{
				if(it1==1)  sprintf(name,"n%s{4,N_{ch}^{rec}}",scnames[ihar]); else  sprintf(name,"n%s{4,#Sigma E_{T}}",scnames[ihar]);
			}

			htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
			if(it1==2){
				if(ihar==1) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
				if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
				if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
			}
			gPad->SetTicks(1,1);
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}

			double shiy=0,shix=0.;
			if(ihar==3) {
				if(it1!=2){shix=0.44;shiy=-0.35;}
				else {shix=0.;shiy=-0.54;}
			}
			if(ihar==2) {
				if(it1!=2){shix=0.44;shiy=-0.4;}
				else {shix=0.44;shiy=-0.4;}
			}
			if(ihar==1) {
				shix=0.4;shiy=0;
			}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();

			if(ihar==3) {
				if(it1!=2){shix=0.1;shiy=0.57;}
				else {shix=0.2;shiy=0.57;}
			}
			if(ihar==2) {shix=0.2;shiy=0;}
			if(ihar==1) {shix=0;shiy=0;}
			if(ihar==3&&it1==2)      text1.DrawLatex(0.35,0.9,ctypes[1]);      else text1.DrawLatex(0.2,0.9,ctypes[1]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}
	//compare between the two bining directly
	// int pads1[]={0,4,2,1,3};
	int pads1a[]={0,3};
	int pads1b[]={4,1};
	for(int itt=0;itt<2;itt++){//Et, Nch and Cent.
		int iet  = pads1a[itt];    int ich = pads1b[itt];    int ib1=0,ib2=0;
		//nc{4} ET vs Nch
		for(int ihar=1;ihar<NHAR;ihar++){
			double mmax=0,mmin=1; TCanvas *ct[10];//look all pT      
			TH1*htmp[6];
			for(int it=0;it<NPTS;it++){
				sprintf(name,"comp_nc4comp_t%d_har%d_pt%d",itt,ihar+1,it);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15); ct[it]=c1;

				if(ihar==1) ib1=0; if(ihar==2) {ib1=2; if(it>2) ib1=3;} if(ihar==3) {ib1=2; if(it>2) ib1=3;}
				if(ihar==1) ib2=0; if(ihar==2) {ib2=2; if(it>2) ib2=3;} if(ihar==3) {ib2=2; if(it>2) ib2=3;}
				//if(ihar==1) ib1=0; if(ihar==2) ib1=2; if(ihar==3) ib1=2;
				//if(ihar==1) ib2=0; if(ihar==2) ib2=2; if(ihar==3) ib2=2;

				//if(itt==0){if(ihar==1) ib1=0; if(ihar==2) ib1=2; if(ihar==3) ib1=2;};
				//if(itt==0){if(ihar==1) ib2=0; if(ihar==2) ib2=2; if(ihar==3) ib2=2;};//mapped

				//if(itt==1){if(ihar==1) ib1=0; if(ihar==2) ib1=2; if(ihar==3) ib1=2;};//mapped
				//if(itt==1){if(ihar==1) ib2=0; if(ihar==2) ib2=2; if(ihar==3) ib2=2;};
				gPad->SetTopMargin(0.01);
				for(int iw=0;iw<2;iw++){
					if(iw==0){//ET
						g1[iw][0] =  (TGraph*)hr4[ib1][ihar][iet][pts[it]][0]->Clone();
						g1[iw][1] =  (TGraph*)hr4[ib1][ihar][iet][pts[it]][1]->Clone();
					}else{//Nch
						g1[iw][0] =  (TGraph*)hr4[ib2][ihar][ich][pts[it]][0]->Clone();
						g1[iw][1] =  (TGraph*)hr4[ib2][ihar][ich][pts[it]][1]->Clone();
					}
					g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
					g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
					g1[iw][1]->SetLineWidth(5);
					int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];
					double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[iw][0]->RemovePoint(N);
						g1[iw][1]->RemovePoint(N);
					}
					th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
					if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[iw][0]->RemovePoint(0);
						g1[iw][1]->RemovePoint(0);
					}
					if(itt==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else {ran1=400/Knee[1];ran2=3000/Knee[1];}
					GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				if(itt==0){//backup
					mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
				}
				gPad->SetTicks(1,1);
				htmp[it] = (TH1*)hcoora[itt]->DrawClone();
				htmp[it]->SetMaximum(mmax1[ihar]*2); htmp[it]->SetMinimum(mmin1[ihar]*1.32);
				sprintf(name,"nc_{%d}{4}",ihar+1);
				htmp[it]->GetYaxis()->SetTitle(name);      htmp[it]->GetXaxis()->SetTitleOffset(1.);
				for(int iw=0;iw<2;iw++){
					g1[iw][1]->Draw("P");
				}
				for(int iw=0;iw<2;iw++){
					g1[iw][0]->Draw("P");
				}
				double shiy=0,shix=0.4;
				if(ihar==1){shix=0.;shiy=-0.2;}
				if(ihar==2){shix=0.45;shiy=-0.4;}
				if(ihar==3){shix=0.45;shiy=-0.3;}
				l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+1); l->AddEntry(g1[0][0],name,"pl");  
				sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+1); l->AddEntry(g1[1][0],name,"pl");  
				l->Draw();
				if(ihar==3) {shix=0.2;shiy=0.5;}
				if(ihar==2) {shix=0.1;shiy=0;}
				if(ihar==1) {shix=0.;shiy=0;}

				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
				if(ihar==3){text1.DrawLatex(0.35,0.92,"Standard method");	text1.DrawLatex(0.35,0.84,ptname1[pts[it]]);}
				else {text1.DrawLatex(0.2,0.92,"Standard method");	text1.DrawLatex(0.2,0.84,ptname1[pts[it]]);}
			}
			for(int it=0;it<NPTS;it++){
				if(ihar==1){ htmp[it]->SetMaximum(mmax1[ihar]*1.2); htmp[it]->SetMinimum(mmin1[ihar]*1.32);}
				if(ihar==2){ htmp[it]->SetMaximum(mmax1[ihar]*2); htmp[it]->SetMinimum(mmin1[ihar]*1.);}
				if(ihar==3){ htmp[it]->SetMaximum(mmax1[ihar]*1.2); htmp[it]->SetMinimum(mmin1[ihar]*1.4);}
				c1 = ct[it]; sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
	}
	for(int itt=0;itt<2;itt++){//Et, Nch and Cent.
		int iet  = pads1a[itt];    int ich = pads1b[itt];    int ib1=0,ib2=0;
		//nc{6} ET vs Nch
		for(int ihar=1;ihar<2;ihar++){
			double mmax=0,mmin=1; TCanvas *ct[10];//look all pT      
			TH1*htmp[6];
			for(int it=0;it<NPTS;it++){
				sprintf(name,"comp_nc6comp_t%d_har%d_pt%d",itt,ihar+1,it);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15); ct[it]=c1;
				ib1=0;ib2=0; if(it>2){ ib1=1;ib2=1;}
				//ib1=1;ib2=1;
				gPad->SetTopMargin(0.01);
				for(int iw=0;iw<2;iw++){
					if(iw==0){//ET
						g1[iw][0] =  (TGraph*)hr6[ib1][ihar][iet][pts[it]][0]->Clone();
						g1[iw][1] =  (TGraph*)hr6[ib1][ihar][iet][pts[it]][1]->Clone();
					}else{//Nch
						g1[iw][0] =  (TGraph*)hr6[ib2][ihar][ich][pts[it]][0]->Clone();
						g1[iw][1] =  (TGraph*)hr6[ib2][ihar][ich][pts[it]][1]->Clone();
					}
					g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
					g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
					g1[iw][1]->SetLineWidth(5);
					int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];
					double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[iw][0]->RemovePoint(N);
						g1[iw][1]->RemovePoint(N);
					}
					th=0.05; 	  if(iw) {th=0.06;}
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[iw][0]->RemovePoint(0);
						g1[iw][1]->RemovePoint(0);
					}
					if(itt==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else {ran1=400/Knee[1];ran2=3000/Knee[1];}
					GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				if(itt==0){//backup
					mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
				}
				gPad->SetTicks(1,1);
				htmp[it] = (TH1*)hcoora[itt]->DrawClone();
				htmp[it]->SetMaximum(mmax1[ihar]*2); htmp[it]->SetMinimum(mmin1[ihar]*1.32);
				sprintf(name,"nc_{%d}{6}",ihar+1);
				htmp[it]->GetYaxis()->SetTitle(name);      htmp[it]->GetXaxis()->SetTitleOffset(1.);
				for(int iw=0;iw<2;iw++){
					g1[iw][1]->Draw("P");
				}
				for(int iw=0;iw<2;iw++){
					g1[iw][0]->Draw("P");
				}
				double shiy=0.1,shix=0.47;
				l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				sprintf(name,"nc_{%d}{6,#Sigma E_{T}}",ihar+1); l->AddEntry(g1[0][0],name,"pl");  
				sprintf(name,"nc_{%d}{6,N_{ch}^{rec}}",ihar+1); l->AddEntry(g1[1][0],name,"pl");  
				l->Draw();
				if(ihar==3) {shix=0.2;shiy=0.5;}
				if(ihar==2) {shix=0.1;shiy=0;}
				if(ihar==1) {shix=0.;shiy=0;}

				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
				text1.DrawLatex(0.2,0.92,"Standard method");	text1.DrawLatex(0.2,0.85,ptname1[pts[it]]);
			}
			for(int it=0;it<NPTS;it++){
				if(ihar==1){ htmp[it]->SetMaximum(mmax1[ihar]*1.32); htmp[it]->SetMinimum(mmin1[ihar]*1.32);}
				c1 = ct[it]; sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
	}
	//nsc & nac 
	for(int itt=0;itt<2;itt++){//Et, Nch and Cent.
		int iet  = pads1a[itt];    int ich = pads1b[itt];    int ib1=0,ib2=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			double mmax=0,mmin=1; TCanvas *ct[10];//look all pT      
			TH1*htmp[6];
			for(int it=0;it<NPTS;it++){
				if(ihar==1)      sprintf(name,"comp_nac_t%d_har%d_pt%d",itt,ihar+1,it); else      sprintf(name,"comp_nsc_t%d_har%d_pt%d",itt,ihar,it);
				c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);ct[it]=c1;
				if(ihar==1) ib1=1; if(ihar==2) {ib1=2; if(it>3) ib1=3;} if(ihar==3) {ib1=2; if(it>3) ib1=3;}
				if(ihar==1) ib2=1; if(ihar==2) {ib2=2; if(it>3) ib2=3;} if(ihar==3) {ib2=2; if(it>3) ib2=3;}
				gPad->SetTopMargin(0.01);
				for(int iw=0;iw<2;iw++){
					if(iw==0){//ET
						if(ihar==1){
							g1[iw][0] =  (TGraph*)hnac[ib1][2][iet][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnac[ib1][2][iet][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hnsc[ib1][ihar][iet][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnsc[ib1][ihar][iet][pts[it]][1]->Clone();
						}
					}else{
						if(ihar==1){
							g1[iw][0] =  (TGraph*)hnac[ib2][2][ich][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnac[ib2][2][ich][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hnsc[ib2][ihar][ich][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnsc[ib2][ihar][ich][pts[it]][1]->Clone();
						}
					}
					g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
					g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
					g1[iw][1]->SetLineWidth(5);

					int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];	
					double th=0.05; if(ihar==3) th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[iw][0]->RemovePoint(N);
						g1[iw][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[iw][0]->RemovePoint(0);
						g1[iw][1]->RemovePoint(0);
					}
					if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
					GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				if(itt==0){//backup
					mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
				}
				gPad->SetTicks(1,1);
				htmp[it] = (TH1*)hcoora[itt]->DrawClone();
				htmp[it]->SetMaximum(mmax1[ihar]*2); htmp[it]->SetMinimum(mmin1[ihar]*1.32);
				if(ihar==1){
					sprintf(name,"n%s{3}",acnames[2]);
				}else{
					sprintf(name,"n%s{4}",scnames[ihar]);
				}
				htmp[it]->GetYaxis()->SetTitle(name);      htmp[it]->GetXaxis()->SetTitleOffset(1.);
				for(int iw=0;iw<2;iw++){
					g1[iw][1]->Draw("P");
				}
				for(int iw=0;iw<2;iw++){
					g1[iw][0]->Draw("P");
				}
				double shiy=0,shix=0.;
				if(ihar==3) { shix=0.4;shiy=-0.35;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				if(ihar==1) { shix=0.4;shiy=0;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				if(ihar==1) sprintf(name,"nac_{2}{3,#Sigma E_{T}}"); if(ihar==2) sprintf(name,"nsc_{2,3}{4,#Sigma E_{T}}"); if(ihar==3) sprintf(name,"nsc_{2,4}{4,#Sigma E_{T}}");
				l->AddEntry(g1[0][0],name,"pl");  
				if(ihar==1) sprintf(name,"nac_{2}{3,N_{ch}^{rec}}"); if(ihar==2) sprintf(name,"nsc_{2,3}{4,N_{ch}^{rec}}"); if(ihar==3) sprintf(name,"nsc_{2,4}{4,N_{ch}^{rec}}");
				l->AddEntry(g1[1][0],name,"pl");        l->Draw();

				if(ihar==3) {
					shix=0.1;shiy=0.57;
				}
				if(ihar==2) {shix=0.2;shiy=0;}
				if(ihar==1) {shix=0;shiy=0;}
				text1.DrawLatex(0.3,0.93,ctypes[0]); 	text1.DrawLatex(0.35,0.88,ptname1[pts[it]]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it=0;it<NPTS;it++){
				if(ihar==1){ htmp[it]->SetMaximum(mmax1[ihar]*1.5); htmp[it]->SetMinimum(-0.05);}
				if(ihar==2){ htmp[it]->SetMaximum(mmax1[ihar]*1.2); htmp[it]->SetMinimum(mmin1[ihar]*1.5);}
				if(ihar==3){ htmp[it]->SetMaximum(mmax1[ihar]*1.5); htmp[it]->SetMinimum(-0.05);}
				c1 = ct[it]; sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
	}
	//nsc & nac  subevent
	for(int itt=0;itt<2;itt++){//Et, Nch and Cent.
		int iet  = pads1a[itt];    int ich = pads1b[itt];    int ib1=0,ib2=0;
		for(int ihar=1;ihar<NHAR;ihar++){
			double mmax=0,mmin=1; TCanvas *ct[10];//look all pT      
			TH1*htmp[6];
			for(int it=0;it<NPTS;it++){
				if(ihar==1)      sprintf(name,"comp_nacsub_t%d_har%d_pt%d",itt,ihar+1,it); else      sprintf(name,"comp_nscsub_t%d_har%d_pt%d",itt,ihar,it);
				c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);ct[it]=c1;
				if(ihar==1) ib1=1; if(ihar==2) {ib1=2; if(it>3) ib1=3;} if(ihar==3) {ib1=2; if(it>3) ib1=3;}
				if(ihar==1) ib2=1; if(ihar==2) {ib2=2; if(it>3) ib2=3;} if(ihar==3) {ib2=2; if(it>3) ib2=3;}
				gPad->SetTopMargin(0.01);
				for(int iw=0;iw<2;iw++){
					if(iw==0){//ET
						if(ihar==1){
							g1[iw][0] =  (TGraph*)hnacb[ib1][2][iet][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnacb[ib1][2][iet][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hnscb[ib1][ihar][iet][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnscb[ib1][ihar][iet][pts[it]][1]->Clone();
						}
					}else{
						if(ihar==1){
							g1[iw][0] =  (TGraph*)hnacb[ib2][2][ich][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnacb[ib2][2][ich][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hnscb[ib2][ihar][ich][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnscb[ib2][ihar][ich][pts[it]][1]->Clone();
						}
					}
					g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
					g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
					g1[iw][1]->SetLineWidth(5);

					int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];	
					double th=0.05; if(ihar==3) th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[iw][0]->RemovePoint(N);
						g1[iw][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[iw][0]->RemovePoint(0);
						g1[iw][1]->RemovePoint(0);
					}
					if(itt==0) { ran1=0.5/Knee[0];ran2=4.4/Knee[0]; } else if(itt==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
					GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				if(itt==0){//backup
					mmax1[ihar]=mmax;	mmin1[ihar]=mmin;
				}
				gPad->SetTicks(1,1);
				htmp[it] = (TH1*)hcoora[itt]->DrawClone();
				htmp[it]->SetMaximum(mmax1[ihar]*2); htmp[it]->SetMinimum(mmin1[ihar]*1.32);
				if(ihar==1){
					sprintf(name,"n%s{3}",acnames[2]);
				}else{
					sprintf(name,"n%s{4}",scnames[ihar]);
				}
				htmp[it]->GetYaxis()->SetTitle(name);      htmp[it]->GetXaxis()->SetTitleOffset(1.);
				for(int iw=0;iw<2;iw++){
					g1[iw][1]->Draw("P");
				}
				for(int iw=0;iw<2;iw++){
					g1[iw][0]->Draw("P");
				}
				double shiy=0,shix=0.;
				if(ihar==3) { shix=0.4;shiy=-0.35;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				if(ihar==1) { shix=0.4;shiy=0;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				if(ihar==1) sprintf(name,"nac_{2}{3,#Sigma E_{T}}"); if(ihar==2) sprintf(name,"nsc_{2,3}{4,#Sigma E_{T}}"); if(ihar==3) sprintf(name,"nsc_{2,4}{4,#Sigma E_{T}}");
				l->AddEntry(g1[0][0],name,"pl");  
				if(ihar==1) sprintf(name,"nac_{2}{3,N_{ch}^{rec}}"); if(ihar==2) sprintf(name,"nsc_{2,3}{4,N_{ch}^{rec}}"); if(ihar==3) sprintf(name,"nsc_{2,4}{4,N_{ch}^{rec}}");
				l->AddEntry(g1[1][0],name,"pl");        l->Draw();

				if(ihar==3) {
					shix=0.1;shiy=0.57;
				}
				if(ihar==2) {shix=0.2;shiy=0;}
				if(ihar==1) {shix=0;shiy=0;}
				text1.DrawLatex(0.3,0.93,ctypes[1]); 	text1.DrawLatex(0.35,0.88,ptname1[pts[it]]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it=0;it<NPTS;it++){
				if(ihar==1){ htmp[it]->SetMaximum(mmax1[ihar]*1.5); htmp[it]->SetMinimum(-0.05);}
				if(ihar==2){ htmp[it]->SetMaximum(mmax1[ihar]*1.2); htmp[it]->SetMinimum(mmin1[ihar]*1.5);}
				if(ihar==3){ htmp[it]->SetMaximum(mmax1[ihar]*1.5); htmp[it]->SetMinimum(-0.05);}
				c1 = ct[it]; sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
	}
	//correlation
	//nc6/4 vs 4/2
	for(int it1=0;it1<2;it1++){//Et, Nch
		int ib=0;
		for(int ihar=0;ihar<3;ihar++){
			sprintf(name,"comp_nc6vsnc4_har%d_%s",ihar+2,otypes[it1]);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				if(ihar==0){ ib=0;}
				if(ihar==1){ ib=3;}
				if(ihar==2){ ib=3;}
				g1[it][0] =  (TGraph*)hcor[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hcor[ib][ihar][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.1; if(ihar==2) th=0.05;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}
				th=0.05; if(ihar==1) th=0.01; if(ihar==2) th=0.02;
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				ggg = (TGraphAsymmErrors*)g1[it][1];	N=g1[it][1]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				th=0.1; if(ihar==1) th=0.1; if(ihar==2) th=0.1;
				N=g1[it][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEXlow()[ibin]>th||ggg->GetEXhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				ran1=0.1; ran2=0.2;
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			TH1*hh= (TH1*)hcoora[it1]->DrawClone();
			if(it1==0){
				sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+2); 
			}else{
				sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+2); 
			}
			hh->GetXaxis()->SetTitle(name);
			if(it1==0){
				sprintf(name,"nc_{%d}{6,#Sigma E_{T}}",ihar+2); 
			}else{
				sprintf(name,"nc_{%d}{6,N_{ch}^{rec}}",ihar+2); 
			}
			hh->GetYaxis()->SetTitle(name);
			hh->GetXaxis()->SetTitleOffset(1.);
			if(ihar==0){//
				hh->GetXaxis()->SetLimits(-0.7,0.2);      hh->SetMaximum(0.7);      hh->SetMinimum(-0.1);
			}else if(ihar==1){
				hh->GetXaxis()->SetLimits(-0.12,0.05);      hh->SetMaximum(0.03);      hh->SetMinimum(-0.023);
			}else{
				hh->GetXaxis()->SetLimits(-0.09,0.3);       hh->SetMaximum(0.05);      hh->SetMinimum(-0.05);
			}
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("P");
			}

			double shiy=0,shix=0.4;
			if(ihar==2) {shiy=-0.4; shix=0.45;}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();
			shix=0;
			if(ihar==2) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.2,0.92,ctypes[0]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}
	for(int it=0;it<NPTS;it++){
		for(int ihar=0;ihar<1;ihar++){
			int ib=0;
			sprintf(name,"comp_nc6vsnc4_har%d_pt%d",ihar+2,it);   c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			gPad->SetTopMargin(0.01);
			for(int it1=0;it1<2;it1++){//Et, Nch
				if(ihar==0){ ib=0;}
				if(ihar==1){ ib=3;}
				if(ihar==2){ ib=3;}
				g1[it1][0] =  (TGraph*)hcor[ib][ihar][it1][pts[it]][0]->Clone();
				g1[it1][1] =  (TGraph*)hcor[ib][ihar][it1][pts[it]][1]->Clone();
				g1[it1][0]->SetMarkerStyle(bsty[it1]); g1[it1][0]->SetMarkerColor(bcol[it1]); g1[it1][0]->SetLineColor(bcol[it1]);
				g1[it1][1]->SetMarkerSize(0); g1[it1][1]->SetMarkerColor(bfcol[it1]); g1[it1][1]->SetLineColor(bfcol[it1]);	g1[it1][1]->SetLineWidth(5);
				int N=g1[it1][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it1][0];
				double th=0.1; if(ihar==2) th=0.05;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it1][0]->RemovePoint(N);
					g1[it1][1]->RemovePoint(N);
				}
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it1][0]->RemovePoint(0);
					g1[it1][1]->RemovePoint(0);
				}
				th=0.05; if(ihar==1) th=0.01; if(ihar==2) th=0.02;
				N=g1[it1][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}
				ggg = (TGraphAsymmErrors*)g1[it1][1];	N=g1[it1][1]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}
				th=0.1; if(ihar==1) th=0.1; if(ihar==2) th=0.1;
				N=g1[it1][0]->GetN();
				for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEXlow()[ibin]>th||ggg->GetEXhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}
				ran1=0.1; ran2=0.2;
				GetMax(g1[it1][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			TH1*hh= (TH1*)hcoora[0]->DrawClone();
			sprintf(name,"nc_{%d}{4}",ihar+2);       hh->GetXaxis()->SetTitle(name);
			sprintf(name,"nc_{%d}{6}",ihar+2);       hh->GetYaxis()->SetTitle(name);
			hh->GetXaxis()->SetTitleOffset(1.);
			if(ihar==0){//
				hh->GetXaxis()->SetLimits(-0.7,0.2);      hh->SetMaximum(0.7);      hh->SetMinimum(-0.1);
			}else if(ihar==1){
				hh->GetXaxis()->SetLimits(-0.12,0.05);      hh->SetMaximum(0.03);      hh->SetMinimum(-0.023);
			}else{
				hh->GetXaxis()->SetLimits(-0.09,0.3);       hh->SetMaximum(0.05);      hh->SetMinimum(-0.05);
			}
			for(int it1=0;it1<2;it1++){//Et, Nch
				g1[it1][1]->Draw("P");
			}
			for(int it1=0;it1<2;it1++){//Et, Nch
				g1[it1][0]->Draw("Pl");
			}

			double shiy=0,shix=0.4;
			if(ihar==2) {shiy=-0.4; shix=0.45;}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			l->AddEntry(g1[0][0],labss[0],"pl");  
			l->AddEntry(g1[1][0],labss[1],"pl");  
			l->Draw();
			shix=0;
			if(ihar==2) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.2,0.92,ctypes[0]); 	text1.DrawLatex(0.35,0.84,ptname1[pts[it]]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}

	for(int it1=0;it1<2;it1++){//Et, Nch
		int ib=0; 
		for(int ihar=0;ihar<5;ihar++){
			if(ihar<3)      sprintf(name,"comp_v64vsv42_t%d_%s",ihar,otypes[it1]); else       sprintf(name,"comp_v42vsv42_t%d_%s",ihar-3,otypes[it1]); 
			c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			gPad->SetTopMargin(0.01);
			for(int it=0;it<NPTS;it++){
				if(ihar==0){ ib=0;if(it>2) ib=1;}
				if(ihar==1){ ib=3;}
				if(ihar==2){ ib=3;}
				if(ihar==3){ ib=2;}
				if(ihar==4){ ib=2;}
				g1[it][0] =  (TGraph*)hcor[ib][ihar+3][it1][pts[it]][0]->Clone();
				g1[it][1] =  (TGraph*)hcor[ib][ihar+3][it1][pts[it]][1]->Clone();
				int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
				double th=0.1; if(ihar==0) th=0.1; else if(ihar<3) th=2; else  th=0.05;
				if(ihar==1&&it==1) g1[it][0]->Print();
				cout<<it<<" HI "<<N<<" ";
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it][0]->RemovePoint(N);
					g1[it][1]->RemovePoint(N);
				}
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it][0]->RemovePoint(0);
					g1[it][1]->RemovePoint(0);
				}

				cout<<" HI "<<N<<" ";
				if(ihar==0) th=0.005; else if(ihar<3) th=2; else th=0.5; //*****
				N=g1[it][0]->GetN();
				if(N>1)for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				ggg = (TGraphAsymmErrors*)g1[it][1];	N=g1[it][1]->GetN();
				cout<<" HI "<<N<<" ";
				if(N>1) for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				if(ihar==0) th=0.03; else if(ihar<3) th=2; else th=0.03; //*****
				ggg = (TGraphAsymmErrors*)g1[it][0];N=g1[it][0]->GetN();
				cout<<" HI "<<N<<endl;
				if(N>1) for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEXlow()[ibin]>th||ggg->GetEXhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}
				ggg = (TGraphAsymmErrors*)g1[it][1];N=g1[it][0]->GetN();
				if(N>1) for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEXlow()[ibin]>th||ggg->GetEXhigh()[ibin]>th){
						g1[it][0]->RemovePoint(ibin);
						g1[it][1]->RemovePoint(ibin);
					}
				}

				ran1=-1; ran2=1;
				GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			TH1*hh= (TH1*)hcoora[it1]->DrawClone();
			if(ihar<3)sprintf(name,"v_{%d}{4,%s}/v_{%d}{2,%s}",ihar+2,labss[it1],ihar+2,labss[it1]); else sprintf(name,"v_{2}{4,%s}/v_{2}{2,%s}",labss[it1],labss[it1]);
			hh->GetXaxis()->SetTitle(name);

			if(ihar<3)sprintf(name,"v_{%d}{6,%s}/v_{%d}{4,%s}",ihar+2,labss[it1],ihar+2,labss[it1]); else sprintf(name,"v_{%d}{4,%s}/v_{%d}{2,%s}",ihar,labss[it1],ihar,labss[it1]);
			hh->GetYaxis()->SetTitle(name);

			hh->GetXaxis()->SetTitleOffset(1.);
			if(ihar==0){//
				hh->GetXaxis()->SetLimits(0.68,0.94);     hh->SetMaximum(1.04);      hh->SetMinimum(0.955);
			}else if(ihar==1){
				hh->GetXaxis()->SetLimits(-0.1,0.3);      hh->SetMaximum(0.03);      hh->SetMinimum(-0.03);
				hh->GetXaxis()->SetLimits(-0.3,0.7);      hh->SetMaximum(5);      hh->SetMinimum(-5);
			}else if(ihar==2){
				hh->GetXaxis()->SetLimits(-1,0.7);      hh->SetMaximum(0.04);      hh->SetMinimum(-0.03);
				hh->GetXaxis()->SetLimits(-1,1);      hh->SetMaximum(5);      hh->SetMinimum(-5);
			}else if(ihar==3){
				hh->GetXaxis()->SetLimits(0.5,0.96);      hh->SetMaximum(0.79);      hh->SetMinimum(0.25);
			}else{
				hh->GetXaxis()->SetLimits(0.6,1);      hh->SetMaximum(1.2);      hh->SetMinimum(-1);
			}
			for(int it=0;it<NPTS;it++){
				g1[it][1]->Draw("P");
			}
			for(int it=0;it<NPTS;it++){
				g1[it][0]->Draw("Pl");
			}

			double shiy=0,shix=0.4;
			if(ihar==2) {shiy=-0.4; shix=0;}
			if(ihar==1) {shiy=-0.4; shix=0.2;}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			for(int it=0;it<NPTS;it++){
				l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
			}
			l->Draw();
			shix=0;
			if(ihar==2) {shix=0.2;shiy=0;}

			text1.DrawLatex(0.2,0.92,ctypes[0]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}

	for(int ihar=0;ihar<1;ihar++){
		for(int it=0;it<NPTS;it++){
			int ib=0; 
			sprintf(name,"comp_v64vsv42_pt%d_har%d",it,ihar+2);      c1 = new TCanvas(name,name,400,400);Divide(c1,1,1,0.,0.15);
			double mmax=0,mmin=1;
			gPad->SetTopMargin(0.01);
			for(int it1=0;it1<2;it1++){//Et, Nch
				if(ihar==0){ ib=0;if(it>2) ib=1;}
				if(ihar==1){ ib=3;}
				if(ihar==2){ ib=3;}
				if(ihar==3){ ib=2;}
				if(ihar==4){ ib=2;}
				g1[it1][0] =  (TGraph*)hcor[ib][ihar+3][it1][pts[it]][0]->Clone();
				g1[it1][1] =  (TGraph*)hcor[ib][ihar+3][it1][pts[it]][1]->Clone();
				g1[it1][0]->SetMarkerStyle(bsty[it1]); g1[it1][0]->SetMarkerColor(bcol[it1]); g1[it1][0]->SetLineColor(bcol[it1]);
				g1[it1][1]->SetMarkerSize(0); g1[it1][1]->SetMarkerColor(bfcol[it1]); g1[it1][1]->SetLineColor(bfcol[it1]);	g1[it1][1]->SetLineWidth(5);

				int N=g1[it1][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it1][0];
				double th=0.1; if(ihar==0) th=0.1;
				while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
					N--;
					g1[it1][0]->RemovePoint(N);
					g1[it1][1]->RemovePoint(N);
				}
				while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
					g1[it1][0]->RemovePoint(0);
					g1[it1][1]->RemovePoint(0);
				}

				if(ihar==0) th=0.005; 
				N=g1[it1][0]->GetN();
				if(N>1)for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}
				ggg = (TGraphAsymmErrors*)g1[it1][1];	N=g1[it1][1]->GetN();
				if(N>1) for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEYlow()[ibin]>th||ggg->GetEYhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}
				if(ihar==0) th=0.03; 
				ggg = (TGraphAsymmErrors*)g1[it1][0];N=g1[it1][0]->GetN();
				if(N>1) for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEXlow()[ibin]>th||ggg->GetEXhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}
				ggg = (TGraphAsymmErrors*)g1[it1][1];N=g1[it1][0]->GetN();
				if(N>1) for(int ibin=N-1;ibin>=0;ibin--){
					if(ggg->GetEXlow()[ibin]>th||ggg->GetEXhigh()[ibin]>th){
						g1[it1][0]->RemovePoint(ibin);
						g1[it1][1]->RemovePoint(ibin);
					}
				}

				ran1=-1; ran2=1;
				GetMax(g1[it1][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
			}
			gPad->SetTicks(1,1);
			TH1*hh= (TH1*)hcoora[0]->DrawClone();
			sprintf(name,"v_{%d}{4}/v_{%d}{2}",ihar+2,ihar+2);      hh->GetXaxis()->SetTitle(name);
			sprintf(name,"v_{%d}{6}/v_{%d}{4}",ihar+2,ihar+2);      hh->GetYaxis()->SetTitle(name);

			hh->GetXaxis()->SetTitleOffset(1.);
			hh->GetXaxis()->SetLimits(0.68,0.94);     hh->SetMaximum(1.04);      hh->SetMinimum(0.955);
			for(int it1=0;it1<2;it1++){
				g1[it1][1]->Draw("P");
			}
			for(int it1=0;it1<2;it1++){
				g1[it1][0]->Draw("Pl");
			}

			double shiy=0,shix=0.4;
			if(ihar==2) {shiy=-0.4; shix=0;}
			if(ihar==1) {shiy=-0.4; shix=0.2;}
			l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
			l->SetFillStyle(0);  l->SetBorderSize(0);
			l->SetTextFont(43);  l->SetTextSize(15);
			l->AddEntry(g1[0][0],labss[0],"pl");  
			l->AddEntry(g1[1][0],labss[1],"pl");  
			l->Draw();
			shix=0;
			text1.DrawLatex(0.2,0.92,ctypes[0]);
			text1.DrawLatex(0.35,0.84,ptname1[pts[it]]);
			lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
			lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
	}
	/*
		 if(it1==0){
		 sprintf(name,"v_{%d}{4,#Sigma E_{T}}/v_{%d}{2,#Sigma E_{T}}",ihar+1,ihar+1); 
		 }else{
		 sprintf(name,"v_{%d}{4,N_{ch}^{rec}}/v_{%d}{2,N_{ch}^{rec}}",ihar+1,ihar+1); 
		 }
		 hh->GetXaxis()->SetTitle(name);
		 if(it1==0){
		 sprintf(name,"v_{%d}{6,#Sigma E_{T}}/v_{%d}{4,#Sigma E_{T}}",ihar+1,ihar+1); 
		 }else{
		 sprintf(name,"v_{%d}{6,N_{ch}^{rec}}/v_{%d}{4,N_{ch}^{rec}}",ihar+1,ihar+1); 
		 }
		 */

	f->Close();
	}
























	void  draw1(){
		//paper plots 
		double ss[][NPTS]={
			{1,0.25, 0.125, 0.03125},
			{1,1.0/3,0.250, 0.25},
			{1,1.0/4,0.0625,1./32},
			{1,0.125,0.0125,0.0125},
		};
		int ssa[][NPTS]={
			{1,4,8, 32},
			{1,3,4, 4},
			{1,4,16,32},
			{1,8,80,80},
		};
		int pads[]={1,2,0};
		//int pads[]={1,0,2};  
		//vn{2}
		for(int ihar=1;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_vn2_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];

			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				int pa = pads[it1]; 
				c1->cd(pa+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hv2[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hv2[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.1; if(ihar==2) th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"v_{%d}{2,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{2,#Sigma E_{T}}",ihar+1);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);

				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}
				if(pa==1){
					line.DrawLine(perc[0][0],0,perc[0][0],0.1);
					line.DrawLine(perc[0][1],0,perc[0][1],0.1);
				}
				if(pa==2){
					line.DrawLine(perc[1][0],0,perc[1][0],0.1);
					line.DrawLine(perc[1][1],0,perc[1][1],0.1);
				}

				l = new TLegend(0.6, 0.72, 1, 0.97);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				double shif=0.15,shix=0.07;
				if(it1==0){
					l->Draw();
				}
				//text1.DrawLatex(0.2,0.9,"Two-subevent method");

				lab1.DrawLatex(0.2+shix,0.1+shif,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.2+shix,0.03+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1) htmp[it1]->SetMaximum(mmax*1.2); else  htmp[it1]->SetMaximum(mmax*1.3); 
				htmp[it1]->SetMinimum(0);	
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//vn{4}
		for(int ihar=1;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_vn4_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hv4[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hv4[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.05; if(ihar==2) th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"v_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{4,#Sigma E_{T}}",ihar+1);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);

				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}
				l = new TLegend(0.6, 0.72, 1, 0.97);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				double shif=0.15,shix=0.07;
				if(it1==0){
					l->Draw();
				}
				text1.DrawLatex(0.2,0.9,ctypes[0]);
				lab1.DrawLatex(0.2+shix,0.09+shif,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.2+shix,0.02+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1) htmp[it1]->SetMaximum(mmax*1.2); else  htmp[it1]->SetMaximum(mmax*1.3); 
				htmp[it1]->SetMinimum(mmin*1.3);	
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//nc{4} subevent a) two panels
		for(int ihar=1;ihar<NHAR;ihar++){
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				int ib=0,ibb=0;
				sprintf(name,"papercomp_nc4sub_%s_har%d",otypes[it1],ihar+1);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[2];
				for(int iw=0;iw<2;iw++){
					c1->cd(iw+1);      gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int it=0;it<NPTS;it++){
						if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
						if(ihar==2){ if(it1==0) ib=2; else ib=2; if(it1==0) ibb=2; else ibb=2;}
						if(ihar==3){ if(it1==0) ib=2; else ib=1; if(it1==0&&it>2) ibb=3; else ibb=2;}
						if(iw==0){
							g1[it][0] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[it][0] =  (TGraph*)hr4b[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hr4b[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.4;ran2=0.99;}
						GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
					htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
					if(it1==1)  sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+1);
					htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
					if(it1==2){
						if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.4,lims[it1]);
						if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.4,lims[it1]);
					}
					for(int it=0;it<NPTS;it++){
						g1[it][1]->Draw("P");
					}
					for(int it=0;it<NPTS;it++){
						g1[it][0]->Draw("P");
					}

					double shiy=0,shix=0.;
					if(ihar==3) {shix=0.4;shiy=-0.2;}
					if(ihar==2) {shix=0.4;shiy=-0.35;}
					l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					for(int it=0;it<NPTS;it++){
						l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
					}
					if(iw==0){
						l->Draw();
					}

					if(ihar==3) {shix=0.1;shiy=0.6;}
					if(ihar==2) {shix=0.2;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					if(ihar==2) shix=0.1;
					text1.DrawLatex(0.2+shix,0.92,ctypes[iw]);
				}
				for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
					if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.5); else  htmp[iw]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin*0.9); else  htmp[iw]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}

		//c1{4} and v1{4} subevent
		for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
			int ib=0,ibb=0;int ihar=0;
			double mmax=0,mmin=1;
			TH1*htmp[2];
			sprintf(name,"papercomp_c4sub_%s_har%d",otypes[it1],ihar+1);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
			for(int iw=0;iw<2;iw++){
				c1->cd(iw+1);    //  gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					if(ihar==0){ib=3;  ibb=3; if(it1==2) {ib=2;ibb=2;}}
					if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
					if(ihar==2){ if(it1==0) ib=2; else ib=2; if(it1==0) ibb=2; else ibb=2;}
					if(ihar==3){ if(it1==0) ib=2; else ib=1; if(it1==0) ibb=2; else ibb=2;}
					//if(it>2) if(ihar==2) ibb=3;
					if(iw==0){
						g1[it][0] =  (TGraph*)hc4[ib][ihar][it1][pts[it]][0]->Clone();
						g1[it][1] =  (TGraph*)hc4[ib][ihar][it1][pts[it]][1]->Clone();
					}else{
						g1[it][0] =  (TGraph*)hc4b[ibb][ihar][it1][pts[it]][0]->Clone();
						g1[it][1] =  (TGraph*)hc4b[ibb][ihar][it1][pts[it]][1]->Clone();
					}
					scalgr(g1[it][0],ss[ihar][it]);	scalgr(g1[it][1],ss[ihar][it]);

					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.05e-6;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.05e-6;
					//if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}

					if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				gPad->SetTicks(1,1);
				htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
				htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"c_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"c_{%d}{4,#Sigma E_{T}}",ihar+1);
				htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==0) htmp[iw]->GetXaxis()->SetLimits(0.55,lims[it1]);
					if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=-0.35,shix=0.35;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.35;}
				l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					sprintf(name,"%s /%d",ptname1[pts[it]],ssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
				}
				if(iw==0){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.6;}
				if(ihar==2) {shix=0.2;shiy=0;}
				shix=0.1;  shiy=0;  
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
				shix=0.3;
				text1.DrawLatex(0.2+shix,0.88,ctypes[iw]);
			}
			for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
				if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.5); else  htmp[iw]->SetMaximum(mmax*1.5);  

				if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin*0.9); else  htmp[iw]->SetMinimum(mmin*1.4); 
				//  if(ihar==0){if(it1==0) htmp[iw]->SetMinimum(mmin*1.);else htmp[iw]->SetMinimum(mmin*1.4); }  
				if(ihar==0){if(it1==0) htmp[iw]->SetMinimum(mmin*1);if(it1==1) htmp[iw]->SetMinimum(mmin*1.4); if(it1==2) htmp[iw]->SetMinimum(mmin*0.8);}

			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
		for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
			int ib=0,ibb=0;int ihar=0;
			double mmax=0,mmin=1;
			TH1*htmp[2];
			sprintf(name,"papercomp_v4sub_%s_har%d",otypes[it1],ihar+1);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
			for(int iw=0;iw<2;iw++){
				c1->cd(iw+1);    //  gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS1;it++){
					if(ihar==0){ib=3;  ibb=3; if(it1==2) {ib=2;ibb=2;}}
					if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
					if(ihar==2){ if(it1==0) ib=2; else ib=2; if(it1==0) ibb=2; else ibb=2;}
					if(ihar==3){ if(it1==0) ib=2; else ib=1; if(it1==0) ibb=2; else ibb=2;}
					//if(it>2) if(ihar==2) ibb=3;
					if(iw==0){
						g1[it][0] =  (TGraph*)hv4[ib][ihar][it1][pts1[it]][0]->Clone();
						g1[it][1] =  (TGraph*)hv4[ib][ihar][it1][pts1[it]][1]->Clone();
					}else{
						g1[it][0] =  (TGraph*)hv4b[ibb][ihar][it1][pts1[it]][0]->Clone();
						g1[it][1] =  (TGraph*)hv4b[ibb][ihar][it1][pts1[it]][1]->Clone();
					}
					//scalgr(g1[it][0],ss[ihar][it]);	scalgr(g1[it][1],ss[ihar][it]);

					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.05;
					//if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}

					if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				gPad->SetTicks(1,1);
				htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
				htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"v_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{4,#Sigma E_{T}}",ihar+1);
				htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==0) htmp[iw]->GetXaxis()->SetLimits(0.55,lims[it1]);
					if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS1;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS1;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=-0.15,shix=0.35;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.35;}
				l = new TLegend(0.2+shix, 0.75+shiy, 0.6+shix, 0.88+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS1;it++){
					//sprintf(name,"%s /%d",ptname1[pts[it]],ssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
					l->AddEntry(g1[it][0],ptname1[pts1[it]],"pl");  
				}
				if(iw==0){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.6;}
				if(ihar==2) {shix=0.2;shiy=0;}
				shix=0.05;  shiy=0.55;  
				lab1.DrawLatex(0.2+shix,0.26+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.2+shix,0.20+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
				text1.DrawLatex(0.2+shix,0.32+shiy,ctypes[iw]);
			}
			for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
				if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.5); else  htmp[iw]->SetMaximum(mmax*1.5);  

				if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin*0.9); else  htmp[iw]->SetMinimum(mmin*1.4); 
				htmp[iw]->SetMaximum(mmax*1.2);
				//if(ihar==0){if(it1==0) htmp[iw]->SetMinimum(mmin*1);if(it1==1) htmp[iw]->SetMinimum(mmin*1.4); if(it1==2) htmp[iw]->SetMinimum(mmin*0.8);}
				htmp[iw]->SetMinimum(0); //htmp[iw]->SetMinimum(-mmax*0.7);
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//cn{4}
		for(int ihar=0;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_c4_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==0) ib=3; if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==0){if(it1<2) ib=3; else ib=3;}
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=2; else ib=2;}
				if(ihar==3){ if(it1==0) ib=2; else ib=1;}
				c1->cd(pads[it1]+1);    
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g2[it][0] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][0]->Clone();

					g1[it][0] =  (TGraph*)hc4[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hc4[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g2[it][0];
					if(ihar){
						double th=0.05; if(ihar==2) th=0.05;
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g2[it][0]->RemovePoint(N);
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g2[it][0]->RemovePoint(0);
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					scalgr(g1[it][0],ss[ihar][it]);	scalgr(g1[it][1],ss[ihar][it]);

					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"c_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"c_{%d}{4,#Sigma E_{T}}",ihar+1);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=-0.3,shix=0.36;
				if(ihar==3) {shix=0.35;shiy=-0.2;}
				if(ihar==2) {shix=0.35;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					sprintf(name,"%s /%d",ptname1[pts[it]],ssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shiy=0;shix=0.1;
				if(ihar==3) {shix=0.1;shiy=0.6;}
				if(ihar==2) {shix=0.2;shiy=0;}

				if(ihar==0)   text1.DrawLatex(0.6,0.88,ctypes[0]);
				if(ihar==1)   text1.DrawLatex(0.2,0.88,ctypes[0]);
				if(ihar==2)   text1.DrawLatex(0.2,0.88,ctypes[0]);
				if(ihar==3)   text1.DrawLatex(0.6,0.40,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5);
				else if(ihar==1) htmp[it1]->SetMaximum(mmax*1.3);  else htmp[it1]->SetMaximum(mmax*2.2);
				if(ihar==0) htmp[it1]->SetMinimum(mmin*1.2); else if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//nc{4}
		for(int ihar=1;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_nc4_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=2; else ib=1;}
				if(ihar==3){ if(it1==0) ib=2; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.05; if(ihar==2) th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+1);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				if(ihar==1 && it1==2) shiy-=0.1;
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}
				shiy=shix=0;
				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.5,0.92,ctypes[0]);
				text1.DrawLatex(0.5,0.92,ctypes[0]);
				lab1.DrawLatex(0.20+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.20+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
		//vn4/vn2
		for(int ihar=1;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_vn4vn2_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=2; else ib=2;}
				if(ihar==3){ if(it1==0) ib=2; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hrat1[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hrat1[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.05; if(ihar==2) th=0.05;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=3.9/Knee[0]; } else if(it1==1) {ran1=2000/Knee[1];ran2=2700/Knee[1];}else {ran1=0.5;ran2=0.97;}
					if(ihar==3){
						//if(it1==0) { ran1=3.;ran2=3.9; } else if(it1==1) {ran1=2000;ran2=2500;}else {ran1=0.85;ran2=0.97;}
					}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"v_{%d}{4,N_{ch}^{rec}}/v_{%d}{2,N_{ch}^{rec}}",ihar+1,ihar+1); 
				else sprintf(name,"v_{%d}{4,#Sigma E_{T}}/v_{%d}{2,#Sigma E_{T}}",ihar+1,ihar+1); 
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.4;
				//if(ihar==3) {shix=0.4;shiy=-0.2;}
				//if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shix=0;
				//if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.2,0.92,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				//if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
				htmp[it1]->SetMaximum(1);
				if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
				if(ihar==3){
					//htmp[it1]->SetMaximum(1);	htmp[it1]->SetMinimum(0);
				}
				if(ihar==2){
					htmp[it1]->SetMinimum(0);
				}
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//nc{6}
		for(int ihar=1;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_nc6_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1<2) ib=1; else ib=0;}
				if(ihar==2){ if(it1==0) ib=3; else ib=3;}
				if(ihar==3){ if(it1==0) ib=3; else ib=3;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hr6[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hr6[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.09; if(ihar==2) th=0.09;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.04; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"nc_{%d}{6,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{6,#Sigma E_{T}}",ihar+1);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.44;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shix=0.0;
				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.2,0.92,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
		//v2{6}
		for(int ihar=1;ihar<2;ihar++){
			int ib=0;
			sprintf(name,"papercomp_vn6_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1<2) ib=1; else ib=0;}
				if(ihar==2){ if(it1==0) ib=3; else ib=3;}
				if(ihar==3){ if(it1==0) ib=3; else ib=3;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hv6[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hv6[ib][ihar][it1][pts[it]][1]->Clone();
					g2[it][0] =  (TGraph*)hr6[ib][ihar][it1][pts[it]][0]->Clone();//for monitoring

					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g2[it][0];
					double th=0.03; if(ihar==2) th=0.09;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g2[it][0]->RemovePoint(N);
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.04; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g2[it][0]->RemovePoint(0);
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(ihar==1&&it1<2){
						//if(it==3&&it1==0) g1[it][0]->Print();
						int flag=1;
						while(flag){
							flag=0;N=g1[it][0]->GetN(); int M=N-5;
							if(M>0){
								double *ey = ((TGraphAsymmErrors*)g1[it][0])->GetEYlow();
								double *y = g1[it][0]->GetY();
								for(int ip=M;ip<N;ip++){
									if(y[ip]<0||(ey[ip]/y[ip])>0.5){
										cout<<ip<<" "<<y[ip]<<" "<<ey[ip]<<endl;
										M = ip; flag=1; break; 
									}
								}
								cout<<flag<<endl;
								if(flag){
									for(int ip=M;ip<N;ip++){
										g1[it][0]->RemovePoint(M);
										g1[it][1]->RemovePoint(M);
									}
								}
							}
						}
						//if(it==3&&it1==0) g1[it][0]->Print();
					}

					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.5;ran2=0.99;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"v_{%d}{6,N_{ch}^{rec}}",ihar+1); else sprintf(name,"v_{%d}{6,#Sigma E_{T}}",ihar+1);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.44;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shix=0.0;
				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.2,0.92,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
		//vn6/vn4
		for(int ihar=1;ihar<4;ihar++){
			int ib=0;
			sprintf(name,"papercomp_vn6vn4_har%d",ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1<2) ib=1; else ib=0;}
				if(ihar==2){ if(it1==0) ib=3; else ib=3;}
				if(ihar==3){ if(it1==0) ib=3; else ib=3;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){ 
					if(ihar==1){
						if(it1==0) {ib=1; if(it>2) ib=2;}
						if(it1==1) {ib=1; if(it>2) ib=2;}
						if(it1==2) {ib=1; if(it>2) ib=2;}
					}
					g1[it][0] =  (TGraph*)hrat2[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hrat2[ib][ihar][it1][pts[it]][1]->Clone();
					//g2[it][0] =  (TGraph*)hr6[ib][ihar][it1][pts[it]][0]->Clone();//for monitoring

					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					double th=0.005;// if(ihar==2) th=0.005;
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-7)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.01; //if(ihar==2) th=0.05; if(ihar==3) th=0.08;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					int flag=1;
					while(flag&&it1<2){
						flag=0;N=g1[it][0]->GetN(); int M=N-5;
						if(M>0){
							double *ey = ((TGraphAsymmErrors*)g1[it][0])->GetEYlow();//((TGraphErrors*)g1[it][0])->GetEY();
							for(int ip=M;ip<N;ip++){
								if(ey[ip]>th){
									M = ip; flag=1; break; 
								}
							}
							cout<<flag<<endl;
							if(flag){
								for(int ip=M;ip<N;ip++){
									g1[it][0]->RemovePoint(M);
									g1[it][1]->RemovePoint(M);
								}
							}
						}
					}
					if(it==1&&it1==0) g1[it][0]->Print();

					if(it1==0) { ran1=2/Knee[0];ran2=3.5/Knee[0]; } else if(it1==1) {ran1=1400/Knee[1];ran2=2000/Knee[1];}else {ran1=0.5;ran2=0.9;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"v_{%d}{6,N_{ch}^{rec}}/v_{%d}{4,N_{ch}^{rec}}",ihar+1,ihar+1); 
				else sprintf(name,"v_{%d}{6,#Sigma E_{T}}/v_{%d}{4,#Sigma E_{T}}",ihar+1,ihar+1); 
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==1) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.4,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}
				line.DrawLine(lims0[it1],1,lims2[it1],1);
				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shix=0.1;
				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.25+shix,0.31+shiy,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				htmp[it1]->SetMaximum(1.009);htmp[it1]->SetMinimum(0.98);
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}





		for(int ihar=1;ihar<2;ihar++){
			int ib=0;     
			sprintf(name,"papercomp_r6_har%d",ihar+1);   c1 = new TCanvas(name,name,900,400);Divide(c1,2,1,0.,0.15);
			if(ihar==1) ib=1; if(ihar==2) ib=3; if(ihar==3) ib=3;
			double maxb=0, minb=0;
			for(int it1=NT-1;it1>=0;it1--){//Et and FCal.
				c1->cd(pads[it1]+1);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				double mmax=0,mmin=1;
				for(int it=0;it<NPT;it++){
					g1[it][0] =  (TGraph*)hr6[ib][ihar][it1][it][0]->Clone();
					g1[it][1] =  (TGraph*)hr6[ib][ihar][it1][it][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					while((ggg->GetEYlow()[N-1]>0.1)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else{ran1=2400/Knee[1];ran2=3000/Knee[1];}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					if(it1==0){ran1=3.6/Knee[0];ran2=4.9/Knee[0];}else {ran1=2400/Knee[1];ran2=3600/Knee[1];}
				}
				gPad->SetTicks(1,1);
				TH1*hh = (TH1*)hcoora[it1]->DrawClone();
				float del = mmax-mmin;      mmax+=0.2*del;      mmin-=1.2*del;
				if(it1){maxb=mmax;minb=mmin;}else {mmax=maxb;mmin=minb;}
				hh->SetMaximum(mmax);  hh->SetMinimum(mmin);	
				if(it1==0) sprintf(name,"nc_{%d}{6,#Sigma E_{T}}",ihar+1); else sprintf(name,"nc_{%d}{6,N_{ch}^{rec}}",ihar+1); hh->GetYaxis()->SetTitle(name);
				// hh->GetXaxis()->SetLimits(ran1,ran2);
				hh->GetXaxis()->SetTitleOffset(1.);;
				for(int it=0;it<NPT;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPT;it++){
					g1[it][0]->Draw("P");
				}
				l = new TLegend(0.7, 0.15, 0.9, 0.40);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPT;it++){
					l->AddEntry(g1[it][0],ptname1[it],"pl");  
				}
				double shif=0.15,shix=0.03;
				if(it1==0){
					l->Draw();
				}
				// lab1.DrawLatex(0.2+shix,0.28+shif,ptname1[it]);
				text1.DrawLatex(0.19+shix,0.2+shif,"Standard method");
				lab1.DrawLatex(0.2+shix,0.12+shif,"#it{#bf{ATLAS}} Internal");
				lab1.DrawLatex(0.2+shix,0.04+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			//sprintf(name,"figs/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		for(int ihar=1;ihar<2;ihar++){
			int ib=0;     
			int ptts[]={0,1,3,4};
			sprintf(name,"papercomp_gauss_har%d",ihar+1);   c1 = new TCanvas(name,name,900,400);Divide(c1,2,1,0.,0.15);
			if(ihar==1) ib=1; if(ihar==2) ib=3; if(ihar==3) ib=3;
			double maxb=0, minb=0;
			for(int it1=NT-1;it1>=0;it1--){//Et and FCal.
				if(it1) ib=1; else ib=2;
				c1->cd(pads[it1]+1);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				double mmax=0,mmin=1;
				for(int it=0;it<3;it++){
					int ipt = ptts[it];
					g1[it][0] =  (TGraph*)hrat2[ib][ihar][it1][ipt][0]->Clone();
					g1[it][1] =  (TGraph*)hrat2[ib][ihar][it1][ipt][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					while((ggg->GetEYlow()[N-1]>0.01)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else{ran1=1000/Knee[1];ran2=2400/Knee[1];}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					if(it1==0){ran1=0.6/Knee[0];ran2=4.9/Knee[0];}else {ran1=100/Knee[1];ran2=3600/Knee[1];}
				}
				gPad->SetTicks(1,1);
				TH1*hh = (TH1*)hcoorb[it1]->DrawClone();
				if(it1) line.DrawLine(100,1,2800,1);    else  line.DrawLine(0,1,4.9,1);
				mmin=0.985;      mmax=1.007;
				if(it1){maxb=mmax;minb=mmin;}else {mmax=maxb;mmin=minb;}
				hh->SetMaximum(mmax);  hh->SetMinimum(mmin);	
				if(it1==0) sprintf(name,"v_{%d}{6,#Sigma E_{T}}/v_{%d}{4,#Sigma E_{T}}",ihar+1,ihar+1);
				else sprintf(name,"v_{%d}{6,N_{ch}^{rec}}/v_{%d}{4,N_{ch}^{rec}}",ihar+1,ihar+1); hh->GetYaxis()->SetTitle(name);
				hh->GetXaxis()->SetTitleOffset(1.);;
				for(int it=0;it<3;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<3;it++){
					g1[it][0]->Draw("P");
				}
				l = new TLegend(0.2, 0.65, 0.5, 0.9);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<3;it++){
					int ipt = ptts[it];
					l->AddEntry(g1[it][0],ptname1[ipt],"pl");  
				}
				double shif=0.15,shix=0.3;
				if(it1==0){
					l->Draw();
				}
				text1.DrawLatex(0.19+shix,0.2+shif,"Standard method");
				lab1.DrawLatex(0.2+shix,0.12+shif,"#it{#bf{ATLAS}} Internal");
				lab1.DrawLatex(0.2+shix,0.04+shif,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			//sprintf(name,"figs/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//sc and ac
		double sss[][NPTS]={
			{1,0.25, 0.125, 0.03125},
			{1,1.0/3,0.250, 0.25},

			{1,1.0/4,0.125,0.1},
			{1,1.0/4,0.125,0.1},
		};
		int sssa[][NPTS]={
			{1,4,8, 32},
			{1,3,4, 4},
			{1,4,8,10},
			{1,4,8,10},
		};

		//sc
		for(int ihar=2;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_sc_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==0) ib=3; if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==0){if(it1<2) ib=3; else ib=3;}
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);    
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hsc[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hsc[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					if(ihar){
						double th=0.05; if(ihar==2) th=0.05;
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g2[it][0]->RemovePoint(N);
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g2[it][0]->RemovePoint(0);
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
					scalgr(g1[it][0],sss[ihar][it]);	scalgr(g1[it][1],sss[ihar][it]);

					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"%s,#Sigma E_{T}",scnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=-0.3,shix=0.36;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.1+shix, 0.7+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					sprintf(name,"%s /%d",ptname1[pts[it]],sssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shiy=0;shix=0.1;
				if(ihar==3) {shix=0.1;shiy=0.59;}
				if(ihar==2) {shix=0.2;shiy=0;}

				if(ihar==2) {
					if(it1<2) text1.DrawLatex(0.6,0.8,ctypes[0]);
					else text1.DrawLatex(0.2,0.85,ctypes[0]);
				}else{
					text1.DrawLatex(0.65,0.85,ctypes[0]);
				}
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5);
				else if(ihar==1) htmp[it1]->SetMaximum(mmax*1.3);  else htmp[it1]->SetMaximum(mmax*2);
				if(ihar==0) htmp[it1]->SetMinimum(mmin*1.2); else if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		//nsc
		for(int ihar=2;ihar<NHAR;ihar++){
			int ib=0;
			sprintf(name,"papercomp_nsc_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];

					double th=0.05; if(ihar==3) th=0.05;
					if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=5.0/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=4000/Knee[1];}else {ran1=0.3;ran2=1;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();

				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",scnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.35;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.57;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.3,0.9,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();

			sprintf(name,"papercomp_nsczoom_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;    mmax=0;mmin=1;
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];

					double th=0.05; if(ihar==3) th=0.05;
					if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					ran1=lims3[it1]; ran2 = lims2[it1];
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",scnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				htmp[it1]->GetXaxis()->SetLimits(lims3[it1],lims[it1]);
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.3;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.57;}
				if(ihar==2) {shix=0.2;shiy=0;}

				text1.DrawLatex(0.3,0.9,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				htmp[it1]->SetMaximum(mmax*1.1); htmp[it1]->SetMinimum(mmin*2);
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();

			sprintf(name,"papercomp_nscsub_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;    mmax=0;mmin=1;
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hnscb[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnscb[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];

					double th=0.05; if(ihar==3) th=0.05;
					if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=5.0/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=4000/Knee[1];}else {ran1=0.3;ran2=1;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",scnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.35;}
				if(ihar==2) {shix=0.4;shiy=-0.4;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.57;}
				if(ihar==2) {shix=0.2;shiy=0;}

				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
				text1.DrawLatex(0.2+shix,0.92,ctypes[1]);
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.3); else  htmp[it1]->SetMinimum(mmin*1.5); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}
		//sc subevent
		for(int ihar=2;ihar<NHAR;ihar++){
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				int ib=0,ibb=0;
				sprintf(name,"papercomp_scsub_%s_har%d",otypes[it1],ihar);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[2];
				for(int iw=0;iw<2;iw++){
					c1->cd(iw+1);      //gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int it=0;it<NPTS;it++){
						if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
						if(ihar==2){ if(it1==0) ib=1; else ib=1; if(it1==0) ibb=1; else ibb=1;}
						if(ihar==3){ if(it1==0) ib=1; else ib=1; if(it1==0&&it>2) ibb=1; else ibb=1;}
						if(iw==0){
							g1[it][0] =  (TGraph*)hsc[ib][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hsc[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[it][0] =  (TGraph*)hscb[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hscb[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.32;ran2=1;}
						scalgr(g1[it][0],sss[ihar][it]);	scalgr(g1[it][1],sss[ihar][it]);
						GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
					htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
					if(it1==1)  sprintf(name,"%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"%s,#Sigma E_{T}",scnames[ihar]);
					htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
					if(it1==2){
						if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
						if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
					}
					for(int it=0;it<NPTS;it++){
						g1[it][1]->Draw("P");
					}
					for(int it=0;it<NPTS;it++){
						g1[it][0]->Draw("P");
					}

					double shiy=0,shix=0.;
					if(ihar==3) {shix=0.4;shiy=-0.2;}
					if(ihar==2) {shix=0.35;shiy=-0.3;}
					l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					for(int it=0;it<NPTS;it++){
						sprintf(name,"%s /%d",ptname1[pts[it]],sssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
					}
					if(iw==0){
						l->Draw();
					}

					if(ihar==3) {shix=0.1;shiy=0.55;}
					if(ihar==2) {shix=0.1;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					if(it1==2){
						text1.DrawLatex(0.3,0.87,ctypes[iw]);
					}else{
						if(ihar==2)text1.DrawLatex(0.5,0.7,ctypes[iw]);
						if(ihar==3)text1.DrawLatex(0.5,0.87,ctypes[iw]);
					}
				}
				if(mmin>0) mmin=0;
				for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
					if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.5); else  htmp[iw]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin*1.3); else  htmp[iw]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
		//nsc subevent
		for(int ihar=2;ihar<NHAR;ihar++){
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				int ib=0,ibb=0;
				sprintf(name,"papercomp_nscsub_%s_har%d",otypes[it1],ihar);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[2];
				for(int iw=0;iw<2;iw++){
					c1->cd(iw+1);      gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int it=0;it<NPTS;it++){
						if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
						if(ihar==2){ if(it1==0) ib=1; else ib=1; if(it1==0) ibb=1; else ibb=1;}
						if(ihar==3){ if(it1==0) ib=1; else ib=1; if(it1==0&&it>2) ibb=1; else ibb=1;}
						if(iw==0){
							g1[it][0] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[it][0] =  (TGraph*)hnscb[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hnscb[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=1;}
						GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
					htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
					if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",scnames[ihar]);
					htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
					if(it1==2){
						if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
						if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
					}
					for(int it=0;it<NPTS;it++){
						g1[it][1]->Draw("P");
					}
					for(int it=0;it<NPTS;it++){
						g1[it][0]->Draw("P");
					}

					double shiy=0,shix=0.;
					if(ihar==3) {shix=0.4;shiy=-0.2;}
					if(ihar==2) {shix=0.4;shiy=-0.35;}
					l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					for(int it=0;it<NPTS;it++){
						l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
					}
					if(iw==0){
						l->Draw();
					}

					if(ihar==3) {shix=0.1;shiy=0.6;}
					if(ihar==2) {shix=0.2;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					if(ihar==2) shix=0.1;
					text1.DrawLatex(0.2+shix,0.92,ctypes[iw]);
				}
				if(mmin>0) mmin=0;
				for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
					if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.5); else  htmp[iw]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin*0.9); else  htmp[iw]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}


		//ac
		for(int ihar=2;ihar<3;ihar++){
			int ib=0;
			sprintf(name,"papercomp_ac_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==0) ib=3; if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==0){if(it1<2) ib=3; else ib=3;}
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);    
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hac[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hac[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					if(ihar){
						double th=0.05; if(ihar==2) th=0.05;
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g2[it][0]->RemovePoint(N);
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g2[it][0]->RemovePoint(0);
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
					scalgr(g1[it][0],sss[ihar][it]);	scalgr(g1[it][1],sss[ihar][it]);

					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"%s,#Sigma E_{T}",acnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0;
				l = new TLegend(0.2+shix, 0.2+shiy, 0.6+shix, 0.4+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					sprintf(name,"%s /%d",ptname1[pts[it]],sssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
				}
				if(it1==2){
					l->Draw();
				}
				shiy=0;shix=0.1;
				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0.64;}

				text1.DrawLatex(0.6,0.74,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.3); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.3);
				else if(ihar==1) htmp[it1]->SetMaximum(mmax*1.1);  else htmp[it1]->SetMaximum(mmax*1.2);
				if(ihar==0) htmp[it1]->SetMinimum(mmin*1.); else if(ihar==1) htmp[it1]->SetMinimum(mmin*1.); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.); else  htmp[it1]->SetMinimum(mmin*1.); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();


			sprintf(name,"papercomp_acsub_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==0) ib=3; if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;    mmax=0;mmin=1;
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==0){if(it1<2) ib=3; else ib=3;}
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);    
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hacb[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hacb[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
					if(ihar){
						double th=0.05; if(ihar==2) th=0.05;
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g2[it][0]->RemovePoint(N);
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.05; if(ihar==3) th=0.08;
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g2[it][0]->RemovePoint(0);
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.3;ran2=0.99;}
					scalgr(g1[it][0],sss[ihar][it]);	scalgr(g1[it][1],sss[ihar][it]);

					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"%s,#Sigma E_{T}",acnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.2,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=-0.3,shix=0.36;
				if(ihar==3) {shix=0.3;shiy=-0.2;}
				if(ihar==2) {shix=0.3;shiy=-0.25;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					sprintf(name,"%s /%d",ptname1[pts[it]],sssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
				}
				if(it1==0){
					l->Draw();
				}
				shiy=0;shix=0.1;
				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0.2;shiy=0.64;}

				text1.DrawLatex(0.5,0.9,ctypes[1]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.3); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.3);
				else if(ihar==1) htmp[it1]->SetMaximum(mmax*1.1);  else htmp[it1]->SetMaximum(mmax*1.2);
				if(ihar==0) htmp[it1]->SetMinimum(mmin*1.); else if(ihar==1) htmp[it1]->SetMinimum(mmin*1.); else if(ihar==2) htmp[it1]->SetMinimum(mmin*1.); else  htmp[it1]->SetMinimum(mmin*1.); 
			}
		}
		//ac subevent 
		for(int ihar=2;ihar<3;ihar++){
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				int ib=0,ibb=0;
				sprintf(name,"papercomp_acsub_%s_har%d",otypes[it1],ihar);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[2];
				for(int iw=0;iw<2;iw++){
					c1->cd(iw+1);      //gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int it=0;it<NPTS;it++){
						if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
						if(ihar==2){ if(it1==0) ib=1; else ib=1; if(it1==0) ibb=1; else ibb=1;}
						if(ihar==3){ if(it1==0) ib=1; else ib=1; if(it1==0&&it>2) ibb=1; else ibb=1;}
						if(iw==0){
							g1[it][0] =  (TGraph*)hac[ib][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hac[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[it][0] =  (TGraph*)hacb[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hacb[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.32;ran2=1;}
						scalgr(g1[it][0],sss[ihar][it]);	scalgr(g1[it][1],sss[ihar][it]);
						GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
					htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
					if(it1==1)  sprintf(name,"%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"%s,#Sigma E_{T}",acnames[ihar]);
					htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
					if(it1==2){
						if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
						if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
					}
					for(int it=0;it<NPTS;it++){
						g1[it][1]->Draw("P");
					}
					for(int it=0;it<NPTS;it++){
						g1[it][0]->Draw("P");
					}

					double shiy=0,shix=0.;
					if(ihar==3) {shix=0.4;shiy=-0.2;}
					if(ihar==2) {shix=0.4;shiy=-0.;}
					l = new TLegend(0.17+shix, 0.6+shiy, 0.55+shix, 0.81+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					for(int it=0;it<NPTS;it++){
						sprintf(name,"%s /%d",ptname1[pts[it]],sssa[ihar][it]);	l->AddEntry(g1[it][0],name,"pl");  
					}
					if(iw==0){
						l->Draw();
					}

					if(ihar==3) {shix=0.1;shiy=0.6;}
					if(ihar==2) {shix=-0.03;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					if(ihar==2) shix=0.1;
					text1.DrawLatex(0.2+shix,0.87,ctypes[1]);
				}
				if(mmin>0) mmin=0;
				for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
					if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.5); else  htmp[iw]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin*1.3); else  htmp[iw]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
		//nac subevent 
		for(int ihar=2;ihar<3;ihar++){
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				int ib=0,ibb=0;
				sprintf(name,"papercomp_nacsub_%s_har%d",otypes[it1],ihar);   c1 = new TCanvas(name,name,800,400);Divide(c1,2,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[2];
				for(int iw=0;iw<2;iw++){
					c1->cd(iw+1);      gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int it=0;it<NPTS;it++){
						if(ihar==1){ if(it1==0) ib=0; else ib=0; if(it1==0) ibb=0; else ibb=0;}
						if(ihar==2){ if(it1==0) ib=1; else ib=1; if(it1==0) ibb=1; else ibb=1;}
						if(ihar==3){ if(it1==0) ib=1; else ib=1; if(it1==0&&it>2) ibb=1; else ibb=1;}
						if(iw==0){
							g1[it][0] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[it][0] =  (TGraph*)hnacb[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[it][1] =  (TGraph*)hnacb[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[it][0]->RemovePoint(N);
							g1[it][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[it][0]->RemovePoint(0);
							g1[it][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.32;ran2=1;}
						GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[iw] = (TH1*)hcoora[it1]->DrawClone();
					htmp[iw]->SetMaximum(mmax);  htmp[it1]->SetMinimum(mmin);	
					if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",acnames[ihar]);
					htmp[iw]->GetYaxis()->SetTitle(name);      htmp[iw]->GetXaxis()->SetTitleOffset(1.);
					if(it1==2){
						if(ihar==2) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
						if(ihar==3) htmp[iw]->GetXaxis()->SetLimits(0.2,lims[it1]);
					}
					for(int it=0;it<NPTS;it++){
						g1[it][1]->Draw("P");
					}
					for(int it=0;it<NPTS;it++){
						g1[it][0]->Draw("P");
					}

					double shiy=0,shix=0.;
					if(ihar==3) {shix=0.4;shiy=-0.2;}
					if(ihar==2) {shix=0.4;shiy=-0.;}
					l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					for(int it=0;it<NPTS;it++){
						l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
					}
					if(iw==0){
						l->Draw();
					}

					if(ihar==3) {shix=0.1;shiy=0.66;}
					if(ihar==2) {shix=-0.03;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					text1.DrawLatex(0.4,0.92,ctypes[iw]);
				}
				if(mmin>0) mmin=0;
				for(int iw=0;iw<2;iw++){//Et, Nch and Cent.
					if(ihar==3) htmp[iw]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[iw]->SetMaximum(mmax*1.2); else  htmp[iw]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[iw]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[iw]->SetMinimum(mmin); else  htmp[iw]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
		//nac
		for(int ihar=2;ihar<3;ihar++){
			int ib=0;
			sprintf(name,"papercomp_nac_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;
			double mmax=0,mmin=1;
			TH1*htmp[NT1];
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];

					double th=0.05; if(ihar==3) th=0.05;
					if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=5.0/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=4000/Knee[1];}else {ran1=0.3;ran2=1;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax*0.8);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",acnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.1;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0;shiy=0;}

				text1.DrawLatex(0.6,0.9,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.2); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.1); else  htmp[it1]->SetMaximum(mmax*1.1);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin); else if(ihar==2) htmp[it1]->SetMinimum(mmin); else  htmp[it1]->SetMinimum(mmin); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();

			sprintf(name,"papercomp_naczoom_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1; mmax=0;mmin=1;
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; if(it1==1) ib=0; if(it1==2) ib=0;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01); gPad->SetLogy();
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];

					double th=0.05; if(ihar==3) th=0.05;
					if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					ran1=lims4[it1]; ran2 = lims2[it1];
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax*0.8);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",acnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				htmp[it1]->GetXaxis()->SetLimits(lims4[it1],lims[it1]);
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				if(ihar==3) {shix=0.4;shiy=-0.2;}
				if(ihar==2) {shix=0.4;shiy=-0.1;}
				l = new TLegend(0.2+shix, 0.72+shiy, 0.6+shix, 0.97+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0;shiy=0;}

				text1.DrawLatex(0.6,0.9,ctypes[0]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				htmp[it1]->SetMaximum(mmax*1.1); htmp[it1]->SetMinimum(0.03);
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();

			sprintf(name,"papercomp_nacsub_har%d",ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
			if(ihar==1) ib=0; if(ihar==2) ib=1; if(ihar==3) ib=1;    mmax=0;mmin=1;
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==1){ if(it1==0) ib=0; else ib=0;}
				if(ihar==2){ if(it1==0) ib=1; else ib=1;}
				if(ihar==3){ if(it1==0) ib=1; else ib=1;}
				c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);
				gPad->SetTicks(1,1);
				for(int it=0;it<NPTS;it++){
					g1[it][0] =  (TGraph*)hnacb[ib][ihar][it1][pts[it]][0]->Clone();
					g1[it][1] =  (TGraph*)hnacb[ib][ihar][it1][pts[it]][1]->Clone();
					int N=g1[it][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[it][0];

					double th=0.05; if(ihar==3) th=0.05;
					if(it1==2) {th=0.05; if(ihar==3) th=0.3;}
					while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
						N--;
						g1[it][0]->RemovePoint(N);
						g1[it][1]->RemovePoint(N);
					}
					th=0.1; if(ihar==3) th=0.15;
					while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
						g1[it][0]->RemovePoint(0);
						g1[it][1]->RemovePoint(0);
					}
					if(it1==0) { ran1=3.6/Knee[0];ran2=5.0/Knee[0]; } else if(it1==1) {ran1=2400/Knee[1];ran2=4000/Knee[1];}else {ran1=0.3;ran2=1;}
					GetMax(g1[it][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					//ran1=lims1[it1];ran2=lims2[it1];
				}
				gPad->SetTicks(1,1);
				htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
				htmp[it1]->SetMaximum(mmax*0.8);  htmp[it1]->SetMinimum(mmin);	
				if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",acnames[ihar]);
				htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);
				if(it1==2){
					if(ihar==2) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
					if(ihar==3) htmp[it1]->GetXaxis()->SetLimits(0.3,lims[it1]);
				}
				for(int it=0;it<NPTS;it++){
					g1[it][1]->Draw("P");
				}
				for(int it=0;it<NPTS;it++){
					g1[it][0]->Draw("P");
				}

				double shiy=0,shix=0.;
				l = new TLegend(0.2+shix, 0.32+shiy, 0.6+shix, 0.52+shiy);
				l->SetFillStyle(0);  l->SetBorderSize(0);
				l->SetTextFont(43);  l->SetTextSize(15);
				for(int it=0;it<NPTS;it++){
					l->AddEntry(g1[it][0],ptname1[pts[it]],"pl");  
				}
				if(it1==2){
					l->Draw();
				}

				if(ihar==3) {shix=0.1;shiy=0.66;}
				if(ihar==2) {shix=0;shiy=0;}

				text1.DrawLatex(0.45,0.9,ctypes[1]);
				lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
				lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
			}
			for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
				if(ihar==3) htmp[it1]->SetMaximum(mmax*1.2); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.1); else  htmp[it1]->SetMaximum(mmax*1.1);  
				if(ihar==1) htmp[it1]->SetMinimum(mmin); else if(ihar==2) htmp[it1]->SetMinimum(mmin); else  htmp[it1]->SetMinimum(mmin); 
			}
			sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
		}

		double val0a[]={0.9,0.9,0.96};
		double val0b[]={0.82,0.82,0.94};
		double val2[]={1.18,1.18,1.003};
		double val3[]={1.22,1.22,1.003};
		double *val1 = val0a;
		//nc{4} subevent b) three panels
		for(int ihar=1;ihar<NHAR;ihar++){
			for(int it=0;it<NPTS;it++){
				int ib=0,ibb=0;
				sprintf(name,"papercomp_nc4sub3_pt%d_har%d",it,ihar+1);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[3];
				for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
					c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int iw=0;iw<2;iw++){
						if(it1==2) {if(ihar==1) ib=ibb=4; else if(ihar==2) ib=ibb=0; else ib=ibb=2; }
						else{
							if(ihar==1){ ib=ibb=0;}
							if(ihar==2){ if(it1==0) ib=ibb=2;   else ib=ibb=1;}
							if(ihar==3){ if(it1==0) {ib=ibb=2;} else ib=ibb=1;}
						}
						if(iw==0){
							g1[iw][0] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hr4[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hr4b[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hr4b[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
						g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
						g1[iw][1]->SetLineWidth(5);

						int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[iw][0]->RemovePoint(N);
							g1[iw][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[iw][0]->RemovePoint(0);
							g1[iw][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.4;ran2=0.99;}
						ran1 = val1[it1];ran2=val2[it1];
						GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
					if(it1==1)  sprintf(name,"nc_{%d}{4,N_{ch}^{rec}}",ihar+1); else sprintf(name,"nc_{%d}{4,#Sigma E_{T}}",ihar+1);
					htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);


					for(int iw=0;iw<2;iw++){
						g1[iw][1]->Draw("P");
					}
					for(int iw=0;iw<2;iw++){
						g1[iw][0]->Draw("P");
					}

					double shiy=-0.4,shix=0.3;
					l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					l->AddEntry(g1[0][0],ctypes[0],"pl");  
					l->AddEntry(g1[1][0],ctypes[1],"pl");  
					l->Draw();

					if(ihar==3) {shix=0.;shiy=0;}
					if(ihar==2) {shix=0.;shiy=0;}
					if(ihar==1) {shix=0.;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					text1.DrawLatex(0.3+shix,0.92,ptname1[pts[it]]);
				}
				for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
					htmp[it1]->GetXaxis()->SetLimits(val1[it1],val3[it1]);
					if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*0.9); else  htmp[it1]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}

		val1 = val0b;
		//nsc{4} subevent b) three panels
		for(int ihar=2;ihar<NHAR;ihar++){
			for(int it=0;it<NPTS;it++){
				int ib=0,ibb=0;
				sprintf(name,"papercomp_nscsub3_pt%d_har%d",it,ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[3];
				for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
					c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int iw=0;iw<2;iw++){
						if(it1==2) {if(ihar==2) {ib=0; ibb=0;} else if(ihar==3) {ib=4; ibb=0;}}
						else if(it1==0) {if(ihar==2) {ib=1; ibb=1;} else if(ihar==3) {ib=0; ibb=1;}}
						else {if(ihar==2) {ib=1; ibb=1;} else if(ihar==3) {ib=0; ibb=1;}}

						if(iw==0){
							g1[iw][0] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnsc[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hnscb[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnscb[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
						g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
						g1[iw][1]->SetLineWidth(5);

						int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[iw][0]->RemovePoint(N);
							g1[iw][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[iw][0]->RemovePoint(0);
							g1[iw][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.4;ran2=0.99;}
						ran1 = val1[it1];ran2=val2[it1];
						GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
					if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",scnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",scnames[ihar]);
					htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);

					for(int iw=0;iw<2;iw++){
						g1[iw][1]->Draw("P");
					}
					for(int iw=0;iw<2;iw++){
						g1[iw][0]->Draw("P");
					}

					double shiy=-0.4,shix=0.3;
					l = new TLegend(0.2+shix, 0.67+shiy, 0.6+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					l->AddEntry(g1[0][0],ctypes[0],"pl");  
					l->AddEntry(g1[1][0],ctypes[1],"pl");  
					l->Draw();

					if(ihar==3) {shix=0.;shiy=0;}
					if(ihar==2) {shix=0.;shiy=0;}
					if(ihar==1) {shix=0.;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					text1.DrawLatex(0.3+shix,0.92,ptname1[pts[it]]);
				}
				for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
					htmp[it1]->GetXaxis()->SetLimits(val1[it1],val3[it1]);
					if(ihar==3) htmp[it1]->SetMaximum(mmax*1.5); else if(ihar==2) htmp[it1]->SetMaximum(mmax*1.5); else  htmp[it1]->SetMaximum(mmax*1.3);  
					if(ihar==1) htmp[it1]->SetMinimum(mmin*1.3); else if(ihar==2) htmp[it1]->SetMinimum(mmin*0.9); else  htmp[it1]->SetMinimum(mmin*1.5); 
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}

		val1 = val0b;
		//nac{4} subevent b) three panels
		for(int ihar=2;ihar<3;ihar++){
			for(int it=0;it<NPTS;it++){
				int ib=0,ibb=0;
				sprintf(name,"papercomp_nacsub3_pt%d_har%d",it,ihar);   c1 = new TCanvas(name,name,1200,400);Divide(c1,3,1,0.,0.15);
				double mmax=0,mmin=1;
				TH1*htmp[3];
				for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
					c1->cd(pads[it1]+1);      gPad->SetTopMargin(0.01);
					gStyle->SetOptTitle(0);
					gStyle->SetOptStat(0);
					gPad->SetTicks(1,1);
					for(int iw=0;iw<2;iw++){
						if(it1==2)      {if(ihar==2) {ib=4; ibb=4;} else if(ihar==3) {ib=4; ibb=0;}}
						else if(it1==0) {if(ihar==2) {ib=0; ibb=0;} else if(ihar==3) {ib=0; ibb=1;}}
						else            {if(ihar==2) {ib=0; ibb=0;} else if(ihar==3) {ib=0; ibb=1;}}

						if(iw==0){
							g1[iw][0] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnac[ib][ihar][it1][pts[it]][1]->Clone();
						}else{
							g1[iw][0] =  (TGraph*)hnacb[ibb][ihar][it1][pts[it]][0]->Clone();
							g1[iw][1] =  (TGraph*)hnacb[ibb][ihar][it1][pts[it]][1]->Clone();
						}
						g1[iw][0]->SetMarkerStyle(bsty[iw]); g1[iw][0]->SetMarkerColor(bcol[iw]); g1[iw][0]->SetLineColor(bcol[iw]);
						g1[iw][1]->SetMarkerSize(0); g1[iw][1]->SetMarkerColor(bfcol[iw]); g1[iw][1]->SetLineColor(bfcol[iw]);
						g1[iw][1]->SetLineWidth(5);

						int N=g1[iw][0]->GetN();	TGraphAsymmErrors*ggg = (TGraphAsymmErrors*)g1[iw][0];
						double th=0.05; if(ihar>1) th=0.06; if(iw) {th=0.07; if(ihar==2) th=0.08; if(ihar==3) th=0.16;}
						while((ggg->GetEYlow()[N-1]>th)||(ggg->GetEYlow()[N-1]<1e-12)){//clean up some points
							N--;
							g1[iw][0]->RemovePoint(N);
							g1[iw][1]->RemovePoint(N);
						}
						th=0.05; if(ihar==2) th=0.08; if(ihar==3) th=0.08;
						if(iw) {th=0.06; if(ihar==2) th=0.08; if(ihar==3) th=0.15;}
						while((ggg->GetEYlow()[0]>th)||(ggg->GetEYlow()[0]<1e-12)){//clean up some points
							g1[iw][0]->RemovePoint(0);
							g1[iw][1]->RemovePoint(0);
						}
						if(it1==0) { ran1=0.6/Knee[0];ran2=4.4/Knee[0]; } else if(it1==1) {ran1=400/Knee[1];ran2=3000/Knee[1];}else {ran1=0.4;ran2=0.99;}
						ran1 = val1[it1];ran2=val2[it1];
						GetMax(g1[iw][0]); if(m1>mmax) mmax = m1; if(m2<mmin) mmin = m2;
					}
					gPad->SetTicks(1,1);
					htmp[it1] = (TH1*)hcoora[it1]->DrawClone();
					if(it1==1)  sprintf(name,"N%s,N_{ch}^{rec}",acnames[ihar]); else  sprintf(name,"N%s,#Sigma E_{T}",acnames[ihar]);
					htmp[it1]->GetYaxis()->SetTitle(name);      htmp[it1]->GetXaxis()->SetTitleOffset(1.);

					for(int iw=0;iw<2;iw++){
						g1[iw][1]->Draw("P");
					}
					for(int iw=0;iw<2;iw++){
						g1[iw][0]->Draw("P");
					}

					double shiy=0,shix=0.3;
					l = new TLegend(0.15+shix, 0.67+shiy, 0.5+shix, 0.88+shiy);
					l->SetFillStyle(0);  l->SetBorderSize(0);
					l->SetTextFont(43);  l->SetTextSize(15);
					l->AddEntry(g1[0][0],ctypes[0],"pl");  
					l->AddEntry(g1[1][0],ctypes[1],"pl");  
					l->Draw();

					if(ihar==3) {shix=0.;shiy=0;}
					if(ihar==2) {shix=0.;shiy=0;}
					if(ihar==1) {shix=0.;shiy=0;}

					lab1.DrawLatex(0.25+shix,0.24+shiy,"#it{#bf{ATLAS}} Internal"); //  lab2.DrawLatex(0.32, 0.42,"Preliminary");
					lab1.DrawLatex(0.25+shix,0.18+shiy,"Pb+Pb 5.02 TeV, 22-470 #mub^{-1}");
					text1.DrawLatex(0.3+shix,0.92,ptname1[pts[it]]);
				}
				for(int it1=0;it1<NT1;it1++){//Et, Nch and Cent.
					htmp[it1]->GetXaxis()->SetLimits(val1[it1],val3[it1]);
					if(ihar==2){ htmp[it1]->SetMaximum(mmax*1.2);  htmp[it1]->SetMinimum(0);}
					else{ htmp[it1]->SetMaximum(mmax*1.);  htmp[it1]->SetMinimum(0);}
				}
				sprintf(name,"figs4/%s.pdf",c1->GetName());	c1->Print(name);  c1->Write();
			}
		}
		// f->Close();
	}


	// 90,117,130,150,170,
