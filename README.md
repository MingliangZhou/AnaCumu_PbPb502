MAIN_binCent: main codes for centrality binning
MAIN_binFCal: main codes for FCal Et binning
MAIN_binNch: main codes for NchRec binning

SYS: systematics checks for centrality binning

The workflow is as flows
	* Rule.h: all the rules for the configuration
	* Tool: tools to calculate event weight, tracking efficiency, fake rates and flattening
	* Event: reading event info from the INPUT
	* MultiCorr: calculation of multi-particle correlations
	* Phase1: wrap up and submit for the Phase1 production
	* Cumu: calculation of cumulants
	* Phase2: wrap up and submit for the cumulant calculation
	* Phase3: calculation of statistical errors and finalization of the results into graphs
	* Plot: making all the monitoring physics plots without systematics
	* Phase4: combine the systematics with physics results

RUN/compile.sh: compile the Phase1
RUN/run_Phase2.sh: compule and run Phase2 with different bin widths
RUN/run_Phase3.sh: compule and run Phase3 with different bin widths
RUN/run_Phase4.sh: compule and run Phase4 with different bin widths
RUN/SUBMIT/submit_Phase1.sh: submit Phase1 production to condor
