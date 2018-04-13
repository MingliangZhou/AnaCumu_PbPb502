#include "MultiCorr.h"

MultiCorr::MultiCorr()
{
	const double nX[4] = {4000, 6000, 4000, 7000};
	const double minX[4] = {0, 0, 0, 0};
	const double maxX[4] = {80, 6, 4000, 7000};
	for(int iT=0; iT<2; iT++) // Cent, FCal, NchRec, Nch
	{
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				sprintf(name,"cvt_Trig%d_Par%d_Par%d",iT,i,j);
				cvt[iT][i][j] = new TProfile(name,"",nX[i],minX[i],maxX[i]);
			}
		}
	}
}

MultiCorr::~MultiCorr()
{
	for(int iT=0; iT<2; iT++) // Cent, FCal, NchRec, Nch
	{
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				delete cvt[iT][i][j];
			}
		}
	}
}

void MultiCorr::fill_1sub(Event* evt)
{
	double par[4] = {0};
	par[0] = evt->evtCent();
	par[1] = evt->evtFCal();
	par[2] = evt->evtNchRec();
	par[3] = evt->evtNch();
	for(int iT=0; iT<2; iT++) // Cent, FCal, NchRec, Nch
	{
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				if(iT==0 && evt->isMB()) cvt[iT][i][j]->Fill(par[i],par[j]);
				if(iT==1 && (evt->isMB() || evt->isHMT())) cvt[iT][i][j]->Fill(par[i],par[j],evt->evtWght());
			}
		}
	}
}

void MultiCorr::fill_1sub_BG(vector<Event*> pool)
{
}

void MultiCorr::fill_3sub(Event* evt)
{
}

void MultiCorr::writeHist(TFile*& fOut)
{
  fOut->cd();
	for(int iT=0; iT<2; iT++) // Cent, FCal, NchRec, Nch
	{
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				cvt[iT][i][j]->Write();
			}
		}
	}
	fOut->Close();
}

double MultiCorr::cal_2pc_1sub(Event* evt, int iV, int iP)
{
	TComplex Q11 = evt->Q[1*iV][iP][1];
	double S21 = pow(evt->S[iP][1],2);
	double S12 = pow(evt->S[iP][2],1);
	return ((Q11*TComplex::Conjugate(Q11)).Re()-S12)/(S21-S12);
}

double MultiCorr::cal_w2pc_1sub(Event* evt, int iP)
{
	double S21 = pow(evt->S[iP][1],2);
	double S12 = pow(evt->S[iP][2],1);
	return S21-S12;
}

double MultiCorr::cal_4pc_1sub(Event* evt, int iV, int iP)
{
	TComplex Q11 = evt->Q[1*iV][iP][1];
	TComplex Q13 = evt->Q[1*iV][iP][3];
	TComplex Q22 = evt->Q[2*iV][iP][2];
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S14 = pow(evt->S[iP][4],1);
	double S21 = pow(evt->S[iP][1],2);
	double S22 = pow(evt->S[iP][2],2);
	double S41 = pow(evt->S[iP][1],4);
	return (pow((Q11*TComplex::Conjugate(Q11)).Re(),2)+(Q22*TComplex::Conjugate(Q22)).Re()-2*(Q22*TComplex::Conjugate(Q11)*TComplex::Conjugate(Q11)).Re()+8*(Q13*TComplex::Conjugate(Q11)).Re()-4*S12*(Q11*TComplex::Conjugate(Q11)).Re()+2*S22-6*S14)/(S41+8*S13*S11-6*S12*S21+3*S22-6*S14);
}

double MultiCorr::cal_w4pc_1sub(Event* evt, int iP)
{
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S14 = pow(evt->S[iP][4],1);
	double S21 = pow(evt->S[iP][1],2);
	double S22 = pow(evt->S[iP][2],2);
	double S41 = pow(evt->S[iP][1],4);
	return S41+8*S13*S11-6*S12*S21+3*S22-6*S14;
}

double MultiCorr::cal_6pc_1sub(Event* evt, int iV, int iP)
{
	TComplex Q11 = evt->Q[1*iV][iP][1];
	TComplex Q13 = evt->Q[1*iV][iP][3];
	TComplex Q15 = evt->Q[1*iV][iP][5];
	TComplex Q22 = evt->Q[2*iV][iP][2];
	TComplex Q24 = evt->Q[2*iV][iP][4];
	TComplex Q33 = evt->Q[3*iV][iP][3];
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S14 = pow(evt->S[iP][4],1);
	double S15 = pow(evt->S[iP][5],1);
	double S16 = pow(evt->S[iP][6],1);
	double S21 = pow(evt->S[iP][1],2);
	double S22 = pow(evt->S[iP][2],2);
	double S23 = pow(evt->S[iP][3],2);
	double S31 = pow(evt->S[iP][1],3);
	double S32 = pow(evt->S[iP][2],3);
	double S41 = pow(evt->S[iP][1],4);
	double S61 = pow(evt->S[iP][1],6);
	return (pow((Q11*TComplex::Conjugate(Q11)).Re(),3)-6*(Q11*TComplex::Conjugate(Q11)).Re()*(Q22*TComplex::Conjugate(Q11)*TComplex::Conjugate(Q11)).Re()+9*(Q22*TComplex::Conjugate(Q22)).Re()*(Q11*TComplex::Conjugate(Q11)).Re()+4*(Q33*TComplex::Conjugate(Q11)*TComplex::Conjugate(Q11)*TComplex::Conjugate(Q11)).Re()+18*S12*(Q22*TComplex::Conjugate(Q11)*TComplex::Conjugate(Q11)).Re()-36*(Q24*TComplex::Conjugate(Q11)*TComplex::Conjugate(Q11)).Re()-36*(Q13*Q11*TComplex::Conjugate(Q22)).Re()+18*S22*(Q11*TComplex::Conjugate(Q11)).Re()-54*S14*(Q11*TComplex::Conjugate(Q11)).Re()-72*S12*(Q13*TComplex::Conjugate(Q11)).Re()+36*(Q13*TComplex::Conjugate(Q13)).Re()+144*(Q15*TComplex::Conjugate(Q11)).Re()-9*S12*pow((Q11*TComplex::Conjugate(Q11)).Re(),2)+36*(Q11*TComplex::Conjugate(Q11)).Re()*(Q13*TComplex::Conjugate(Q11)).Re()-9*S12*(Q22*TComplex::Conjugate(Q22)).Re()+36*(Q24*TComplex::Conjugate(Q22)).Re()-12*(Q33*TComplex::Conjugate(Q22)*TComplex::Conjugate(Q11)).Re()+4*(Q33*TComplex::Conjugate(Q33)).Re()+54*S14*S12-6*S32-120*S16)/(S61-15*S12*S41+40*S13*S31+45*S22*S21-90*S14*S21-120*S13*S12*S11-15*S32+144*S15*S11+90*S14*S12+40*S23-120*S16);
}

double MultiCorr::cal_w6pc_1sub(Event* evt, int iP)
{
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S14 = pow(evt->S[iP][4],1);
	double S15 = pow(evt->S[iP][5],1);
	double S16 = pow(evt->S[iP][6],1);
	double S21 = pow(evt->S[iP][1],2);
	double S22 = pow(evt->S[iP][2],2);
	double S23 = pow(evt->S[iP][3],2);
	double S31 = pow(evt->S[iP][1],3);
	double S32 = pow(evt->S[iP][2],3);
	double S41 = pow(evt->S[iP][1],4);
	double S61 = pow(evt->S[iP][1],6);
	return S61-15*S12*S41+40*S13*S31+45*S22*S21-90*S14*S21-120*S13*S12*S11-15*S32+144*S15*S11+90*S14*S12+40*S23-120*S16;
}

double MultiCorr::cal_2pc_1sub_BG(Event* evt1, Event* evt2, int iV, int iP)
{
	TComplex Q11_1 = evt1->Q[1*iV][iP][1];
	TComplex Q11_2 = evt2->Q[1*iV][iP][1];
	double S11_1 = pow(evt1->S[iP][1],1);
	double S11_2 = pow(evt2->S[iP][1],1);
	return (Q11_1*TComplex::Conjugate(Q11_2)).Re()/(S11_1*S11_2);
}

double MultiCorr::cal_w2pc_1sub_BG(Event* evt1, Event* evt2, int iP)
{
	double S11_1 = pow(evt1->S[iP][1],1);
	double S11_2 = pow(evt2->S[iP][1],1);
	return S11_1*S11_2;
}

double MultiCorr::cal_4pc_1sub_BG(Event* evt1, Event* evt2, Event* evt3, Event* evt4, int iV, int iP)
{
	TComplex Q11_1 = evt1->Q[1*iV][iP][1];
	TComplex Q11_2 = evt2->Q[1*iV][iP][1];
	TComplex Q11_3 = evt3->Q[1*iV][iP][1];
	TComplex Q11_4 = evt4->Q[1*iV][iP][1];
	double S11_1 = pow(evt1->S[iP][1],1);
	double S11_2 = pow(evt2->S[iP][1],1);
	double S11_3 = pow(evt3->S[iP][1],1);
	double S11_4 = pow(evt4->S[iP][1],1);
	return (Q11_1*Q11_2*TComplex::Conjugate(Q11_3)*TComplex::Conjugate(Q11_4)).Re()/(S11_1*S11_2*S11_3*S11_4);
}

double MultiCorr::cal_w4pc_1sub_BG(Event* evt1, Event* evt2, Event* evt3, Event* evt4, int iP)
{
	double S11_1 = pow(evt1->S[iP][1],1);
	double S11_2 = pow(evt2->S[iP][1],1);
	double S11_3 = pow(evt3->S[iP][1],1);
	double S11_4 = pow(evt4->S[iP][1],1);
	return S11_1*S11_2*S11_3*S11_4;
}

double MultiCorr::cal_4psc_1sub(Event* evt, int iV, int iP)
{
	TComplex M11;
	TComplex M13;
	TComplex N11;
	TComplex N13;
	TComplex MN2;
	TComplex NM2;
	if(iV==0) // (1+2-1-2)
	{
		M11 = evt->Q[1][iP][1];
		M13 = evt->Q[1][iP][3];
		N11 = evt->Q[2][iP][1];
		N13 = evt->Q[2][iP][3];
		MN2 = evt->Q[3][iP][2];
		NM2 = evt->Q[1][iP][2];
	}
	if(iV==1) // (1+3-1-3)
	{
		M11 = evt->Q[1][iP][1];
		M13 = evt->Q[1][iP][3];
		N11 = evt->Q[3][iP][1];
		N13 = evt->Q[3][iP][3];
		MN2 = evt->Q[4][iP][2];
		NM2 = evt->Q[2][iP][2];
	}
	if(iV==2) // (2+3-2-3)
	{
		M11 = evt->Q[2][iP][1];
		M13 = evt->Q[2][iP][3];
		N11 = evt->Q[3][iP][1];
		N13 = evt->Q[3][iP][3];
		MN2 = evt->Q[5][iP][2];
		NM2 = evt->Q[1][iP][2];
	}
	if(iV==3) // (2+4-2-4)
	{
		M11 = evt->Q[2][iP][1];
		M13 = evt->Q[2][iP][3];
		N11 = evt->Q[4][iP][1];
		N13 = evt->Q[4][iP][3];
		MN2 = evt->Q[6][iP][2];
		NM2 = evt->Q[2][iP][2];
	}
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S14 = pow(evt->S[iP][4],1);
	double S21 = pow(evt->S[iP][1],2);
	double S22 = pow(evt->S[iP][2],2);
	double S41 = pow(evt->S[iP][1],4);
	return ((M11*TComplex::Conjugate(M11)*N11*TComplex::Conjugate(N11)).Re()-2*(MN2*TComplex::Conjugate(M11)*TComplex::Conjugate(N11)).Re()-2*(NM2*M11*TComplex::Conjugate(N11)).Re()+(MN2*TComplex::Conjugate(MN2)).Re()+(NM2*TComplex::Conjugate(NM2)).Re()+4*(M13*TComplex::Conjugate(M11)+N13*TComplex::Conjugate(N11)).Re()-S12*(M11*TComplex::Conjugate(M11)+N11*TComplex::Conjugate(N11)).Re()+S22-6*S14)/(S41+8*S13*S11-6*S12*S21+3*S22-6*S14);
}

double MultiCorr::cal_3pac_1sub(Event* evt, int iV, int iP)
{
	TComplex M11;
	TComplex M12;
	TComplex N11;
	TComplex N12;
	TComplex MN1;
	TComplex MN2;
	if(iV==0) // (1+1-2)
	{
		M11 = evt->Q[1][iP][1];
		M12 = evt->Q[1][iP][2];
		N11 = evt->Q[1][iP][1];
		N12 = evt->Q[1][iP][2];
		MN1 = evt->Q[2][iP][1];
		MN2 = evt->Q[2][iP][2];
	}
	if(iV==1) // (1+2-3)
	{
		M11 = evt->Q[1][iP][1];
		M12 = evt->Q[1][iP][2];
		N11 = evt->Q[2][iP][1];
		N12 = evt->Q[2][iP][2];
		MN1 = evt->Q[3][iP][1];
		MN2 = evt->Q[3][iP][2];
	}
	if(iV==2) // (2+2-4)
	{
		M11 = evt->Q[2][iP][1];
		M12 = evt->Q[2][iP][2];
		N11 = evt->Q[2][iP][1];
		N12 = evt->Q[2][iP][2];
		MN1 = evt->Q[4][iP][1];
		MN2 = evt->Q[4][iP][2];
	}
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S31 = pow(evt->S[iP][1],3);
	return ((M11*N11*TComplex::Conjugate(MN1)).Re()-(MN1*TComplex::Conjugate(MN2)).Re()-(M11*TComplex::Conjugate(M12)).Re()-(N11*TComplex::Conjugate(N12)).Re()+2*S13)/(S31-3*S12*S11+2*S13);
}

double MultiCorr::cal_w3pc_1sub(Event* evt, int iP)
{
	double S11 = pow(evt->S[iP][1],1);
	double S12 = pow(evt->S[iP][2],1);
	double S13 = pow(evt->S[iP][3],1);
	double S31 = pow(evt->S[iP][1],3);
	return S31-3*S12*S11+2*S13;
}

double MultiCorr::cal_2pc_3sub(Event* evt, int iV, int iP, int iMode)
{
	TComplex A11 = evt->Q3sub[1*iV][iP][0][1];
	TComplex B11;
	if(iMode==1) B11 = evt->Q3sub[1*iV][iP][1][1];
	if(iMode==2) B11 = evt->Q3sub[1*iV][iP][2][1];
	double a11 = pow(evt->nTrk3sub[iP][0][1],1);
	double b11 = 0;
	if(iMode==1) b11 = pow(evt->nTrk3sub[iP][1][1],1);
	if(iMode==2) b11 = pow(evt->nTrk3sub[iP][2][1],1);
	return (A11*TComplex::Conjugate(B11)).Re()/(a11*b11);
}

double MultiCorr::cal_w2pc_3sub(Event* evt, int iP, int iMode)
{
	double a11 = pow(evt->nTrk3sub[iP][0][1],1);
	double b11 = 0;
	if(iMode==1) b11 = pow(evt->nTrk3sub[iP][1][1],1);
	if(iMode==2) b11 = pow(evt->nTrk3sub[iP][2][1],1);
	return a11*b11;
}

double MultiCorr::cal_4pc_3sub(Event* evt, int iV, int iP)
{
	TComplex A11 = evt->Q3sub[1*iV][iP][0][1];
	TComplex A22 = evt->Q3sub[2*iV][iP][0][2];
	TComplex B11 = evt->Q3sub[1*iV][iP][1][1];
	TComplex C11 = evt->Q3sub[1*iV][iP][2][1];
	double a12 = pow(evt->nTrk3sub[iP][0][2],1);
	double a21 = pow(evt->nTrk3sub[iP][0][1],2);
	double b11 = pow(evt->nTrk3sub[iP][1][1],1);
	double c11 = pow(evt->nTrk3sub[iP][2][1],1);
	return ((A11*TComplex::Conjugate(B11)*A11*TComplex::Conjugate(C11)).Re()-(A22*TComplex::Conjugate(B11)*TComplex::Conjugate(C11)).Re())/((a21-a12)*b11*c11);
}

double MultiCorr::cal_w4pc_3sub(Event* evt, int iP)
{
	double a12 = pow(evt->nTrk3sub[iP][0][2],1);
	double a21 = pow(evt->nTrk3sub[iP][0][1],2);
	double b11 = pow(evt->nTrk3sub[iP][1][1],1);
	double c11 = pow(evt->nTrk3sub[iP][2][1],1);
	return (a21-a12)*b11*c11;
}

double MultiCorr::cal_4psc_3sub(Event* evt, int iV, int iP)
{
	TComplex M11;
	TComplex N11;
	TComplex MN2;
	TComplex B11;
	TComplex C11;
	if(iV==0) // (1+2-1-2)
	{
		M11 = evt->Q3sub[1][iP][0][1];
		N11 = evt->Q3sub[2][iP][0][1];
		MN2 = evt->Q3sub[3][iP][0][2];
		B11 = evt->Q3sub[1][iP][1][1];
		C11 = evt->Q3sub[2][iP][2][1];
	}
	if(iV==1) // (1+3-1-3)
	{
		M11 = evt->Q3sub[1][iP][0][1];
		N11 = evt->Q3sub[3][iP][0][1];
		MN2 = evt->Q3sub[4][iP][0][2];
		B11 = evt->Q3sub[1][iP][1][1];
		C11 = evt->Q3sub[3][iP][2][1];
	}
	if(iV==2) // (2+3-2-3)
	{
		M11 = evt->Q3sub[2][iP][0][1];
		N11 = evt->Q3sub[3][iP][0][1];
		MN2 = evt->Q3sub[5][iP][0][2];
		B11 = evt->Q3sub[2][iP][1][1];
		C11 = evt->Q3sub[3][iP][2][1];
	}
	if(iV==3) // (2+4-2-4)
	{
		M11 = evt->Q3sub[2][iP][0][1];
		N11 = evt->Q3sub[4][iP][0][1];
		MN2 = evt->Q3sub[6][iP][0][2];
		B11 = evt->Q3sub[2][iP][1][1];
		C11 = evt->Q3sub[4][iP][2][1];
	}
	double a12 = pow(evt->nTrk3sub[iP][0][2],1);
	double a21 = pow(evt->nTrk3sub[iP][0][1],2);
	double b11 = pow(evt->nTrk3sub[iP][1][1],1);
	double c11 = pow(evt->nTrk3sub[iP][2][1],1);
	return ((M11*TComplex::Conjugate(B11)*N11*TComplex::Conjugate(C11)).Re()-(MN2*TComplex::Conjugate(B11)*TComplex::Conjugate(C11)).Re())/((a21-a12)*b11*c11);
}

double MultiCorr::cal_3pac_3sub(Event* evt, int iV, int iP)
{
	TComplex A11;
	TComplex B11;
	TComplex C11;
	if(iV==0) // (1+1-2)
	{
		A11 = evt->Q3sub[1][iP][0][1];
		B11 = evt->Q3sub[1][iP][1][1];
		C11 = evt->Q3sub[2][iP][2][1];
	}
	if(iV==1) // (1+2-3)
	{
		A11 = evt->Q3sub[1][iP][0][1];
		B11 = evt->Q3sub[2][iP][1][1];
		C11 = evt->Q3sub[3][iP][2][1];
	}
	if(iV==2) // (2+2-4)
	{
		A11 = evt->Q3sub[2][iP][0][1];
		B11 = evt->Q3sub[2][iP][1][1];
		C11 = evt->Q3sub[4][iP][2][1];
	}
	double a11 = pow(evt->nTrk3sub[iP][0][1],1);
	double b11 = pow(evt->nTrk3sub[iP][1][1],1);
	double c11 = pow(evt->nTrk3sub[iP][2][1],1);
	return ((A11*B11*TComplex::Conjugate(C11)).Re())/(a11*b11*c11);
}

double MultiCorr::cal_w3pc_3sub(Event* evt, int iP)
{
	double a11 = pow(evt->nTrk3sub[iP][0][1],1);
	double b11 = pow(evt->nTrk3sub[iP][1][1],1);
	double c11 = pow(evt->nTrk3sub[iP][2][1],1);
	return a11*b11*c11;
}

void MultiCorr::swap(double& val1, double& val2, double& val3)
{
	double tmp = 0;
	tmp = val1;
	val1 = val2;
	val2 = val3;
	val3 = tmp;
}

void MultiCorr::swap(TComplex& val1, TComplex& val2, TComplex& val3)
{
	TComplex tmp;
	tmp = val1;
	val1 = val2;
	val2 = val3;
	val3 = tmp;
}

void MultiCorr::swap(Event*& evt1, Event*& evt2, Event*& evt3)
{
	Event* tmp = new Event(*evt1);
	*evt1 = *evt2;
	*evt2 = *evt3;
	*evt3 = *tmp;
	delete tmp;
}

void MultiCorr::initHist_VP(TH1D*& hIn, const char* hName, int iV, int iP)
{
	sprintf(name,"%s_Har%d_Pt%d",hName,iV,iP);
	hIn = new TH1D(name,"",nBin,minBin,maxBin);
}

void MultiCorr::initHist_AVP(TH1D*& hIn, const char* hName, int iA, int iV, int iP)
{
	sprintf(name,"%s_Pm%d_Har%d_Pt%d",hName,iA,iV,iP);
	hIn = new TH1D(name,"",nBin,minBin,maxBin);
}
