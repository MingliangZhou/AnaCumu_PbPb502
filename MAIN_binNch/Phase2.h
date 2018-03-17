#ifndef Phase2_H_
#define Phase2_H_

#include "Cumu.h"

Cumu* cumu;

class Phase2
{
	public:
		Phase2(unsigned int iBin, unsigned int iSample);
		~Phase2();
		void initialize(unsigned int iBin, unsigned int iSample);
		void execute();
		void finalize(unsigned int iBin, unsigned int iSample);
};

#endif
