#ifndef Phase1_H_
#define Phase1_H_

#include "MultiCorr.h"

TChain* treeIn;
Tool* tool;
MultiCorr* corr;

class Phase1
{
	private:
		bool detMixZvtxCent(double zVtx, double cent, int& tagZvtx, int& tagCent); // determine mixed event
		vector<Event*> EventPool[nMixZvtx][nMixCent]; // mixed event pool

	public:
		Phase1(int);
		~Phase1();
		void initialize(int);
		void execute(int);
		void finalize(int);
};

#endif
