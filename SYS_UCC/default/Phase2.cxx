#include "Phase2.h"

Phase2::Phase2(unsigned int iBin, unsigned int iSample)
{
	initialize(iBin,iSample);
	execute();
	finalize(iBin,iSample);
}

Phase2::~Phase2()
{
	delete cumu;
}

void Phase2::execute()
{
	cout<<"execute..."<<endl;

	cumu->cal_all();
}

void Phase2::initialize(unsigned int iBin, unsigned int iSample)
{
	cout<<"initialize..."<<endl;

	if(iSample<nSample) sprintf(name,"../OUTPUT/Phase1/Phase1_Sample%d.root",iSample);
	else sprintf(name,"../OUTPUT/Phase1/Phase1_All.root");
	TFile* fIn = new TFile(name,"READ");
	cumu = new Cumu(fIn, iBin);
}

void Phase2::finalize(unsigned int iBin, unsigned int iSample)
{
	cout<<"finalize..."<<endl;

	if(iSample<nSample) sprintf(name,"../OUTPUT/Phase2/bin%d/Phase2_Sample%d.root",iBin,iSample);
	else sprintf(name,"../OUTPUT/Phase2/bin%d/Phase2_All.root",iBin);
	TFile* fOut = new TFile(name,"RECREATE");
	cumu->writeHist(fOut);
}



