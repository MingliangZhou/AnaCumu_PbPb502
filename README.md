MAIN_binCent: main codes for centrality binning

The workflow is as flows
	* Rule.h: all the rules for the configuration
	* Tool: tools to calculate event weight, tracking efficiency, fake rates and flattening
	* Event: reading event info from the INPUT
	* MultiCorr: calculation of multi-particle correlations
	* Phase1: wrap up and submit for the Phase1 production
	* Cumu: calculation of cumulants
	* Phase2: wrap up and submit for the cumulant calculation
	* Phase3: calculation of statistical errors and finalization of the results into graphs
