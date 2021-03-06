#include "Phase1.h"

Phase1::Phase1(int iJob)
{
	initialize(iJob);
	execute(iJob);
	finalize(iJob);
}

Phase1::~Phase1()
{
	delete corr;
}

void Phase1::execute(int iJob)
{
	cout<<"execute..."<<endl;

	Event* evt = new Event(treeIn);
	for(int iEvt=0.9*treeIn->GetEntries(); iEvt<treeIn->GetEntries(); iEvt++)
	{
		if(iEvt%1000==0)
		{
			for(int i=0; i<100.*iEvt/treeIn->GetEntries(); i++) cout<<"#";
			cout<<">>>";
			for(int i=100.*iEvt/treeIn->GetEntries(); i<100-1; i++) cout<<"-";
			cout<<'\r';
		}
		treeIn->GetEntry(iEvt);

		if(!(evt->isGoodEvent())) continue;

		// speed up without mixed events
		//int tagMixZvtx = 0; int tagMixCent = 0;
		//if(!detMixZvtxCent(evt->evtZvtx(),evt->evtCls(),tagMixZvtx,tagMixCent)) continue;

		//Event* evtTmp = new Event(*evt);
		//EventPool[tagMixZvtx][tagMixCent].push_back(evtTmp);
		//if(EventPool[tagMixZvtx][tagMixCent].size()<nDepth) continue;

		corr->fill_1sub(evt);
		//corr->fill_1sub_BG(EventPool[tagMixZvtx][tagMixCent]);
		//corr->fill_3sub(evt);

		//delete EventPool[tagMixZvtx][tagMixCent].at(0);
		//EventPool[tagMixZvtx][tagMixCent].erase( EventPool[tagMixZvtx][tagMixCent].begin() );
	}
	delete evt;
}

void Phase1::initialize(int iJob)
{
	cout<<"initialize..."<<endl;

	treeIn = new TChain("HeavyIonD3PD");
	ifstream lis("../../../../MAIN_binCent/INPUT/flist_PbPb502.txt");
	int cnt = -1;
	while(!lis.eof())
	{
		cnt ++;
		string fName;
		lis >> fName;
		if(cnt<int(1141.*iJob/100) || cnt>=int(1141.*(iJob+1)/100)) continue;
		cout<<fName.c_str()<<endl;
		if(!fName.empty()) treeIn->Add(fName.c_str());
		else break;
	}

	corr = new MultiCorr();
}

void Phase1::finalize(int iJob)
{
	cout<<"finalize..."<<endl;

	sprintf(name,"../../OUTPUT/Phase1/Phase1_Job%d.root",iJob);
	TFile* fOut = new TFile(name,"RECREATE");
	corr->writeHist(fOut);
}

bool Phase1::detMixZvtxCent(double zVtx, double cent, int& tagZvtx, int& tagCent)
{
	tagZvtx = -1;
	tagCent = -1;

	tagZvtx = int((zVtx+cutZvtx)/(2*cutZvtx/nMixZvtx));
	tagCent = int(cent/(100./nMixCent));

	if(tagZvtx<0 || tagZvtx>=int(nMixZvtx)) return false;
	if(tagCent<0 || tagCent>=int(nMixCent)) return false;

	return true;
}
