TH3D* initializeEfficiency(TFile* file){

  TH3D* matchedReco = (TH3D*)file->Get("matchedRecoTruePt");

  allTrue = rebinEff(allTrue,"allTrue_rebinned");
  matchedReco = rebinEff(matchedReco,"matchedReco_rebinned");
  matchedReco->Divide(allTrue);
  TH3D* efficiency = matchedReco;

  TH3D* allReco = (TH3D*)file->Get("allReco");
  TH3D* unMatchedReco = (TH3D*)file->Get("unmatchedRecoToPrimary");

  allReco = rebinEff(allReco,"allReco_rebinned");
  unMatchedReco = rebinEff(unMatchedReco,"unMatchedReco_rebinned");
  unMatchedReco->Divide(allReco);
  TH3D* fakes = unMatchedReco;

  TF1* f1 = new TF1("f1","1",-100,100);
  TH3D* weight = fakes;
  weight->Scale(-1);
  weight->Add(f1,1);
  weight->Divide(efficiency);
  weight->SetName("weight");

  return weight;

}

TH3D* rebinEff(TH3D* hist,const char* name){

  double pt_bins[] = {0.5,0.6,0.8,1.0,
                      1.2,1.5,2.0, 3.0, 4.0,
                      5.0, 8.0,
                      60};
  double eta_bins[] = {-2.5,-2.4,-2.3,   
		       -2.2,-2.1,-2.0,   
		       -1.9,-1.8,-1.7,   
		       -1.6,-1.5,-1.4, 
		       -1.3,-1.2,-1.1, 
		       -1.0,-0.9,-0.8, 
		       -0.7,-0.6,-0.5, 
		       -0.4,-0.3,-0.2, 
		       -0.1,0,0.1,0.2,
		       0.3,0.4,0.5,0.6,
		       0.7,0.8,0.9, 1.0,
		       1.1,1.2, 1.3,1.4,
		       1.5,  1.6,1.7,1.8,  
		       1.9,2.0,2.1,   2.2,2.3,2.4,2.5};

  double cent_bin[] = {
    0.063719,//80%                                                                                          
    0.097388, 0.14414,    0.207148, 0.289595,//60                                   
    0.394518, 0.525092,    0.684377,    0.87541,1.10211,//35                                    
    1.36875,    1.68058,//25                                                                             
    2.04651,//20%                                                           
    2.12711,    2.21002,    2.29572,
    2.38468,    2.47658,    2.57162,    2.66999,    2.77237,    2.87864,
    2.98931,    3.10407,    3.22397,    3.34945,    3.48077,    3.61844,
    3.7635,    3.91763,    4.08137,    4.26258 ,4.54199 ,6.};

  TH3D* dest = new TH3D(name,name,80,0,6,400,0,40,50,-2.5,2.5);
  dest->SetBins( 33,cent_bin,11,pt_bins,50,eta_bins);

  copyToReducedBinning(hist,dest);

  return dest;

}


void copyToReducedBinning(const TH3D* src, TH3D* dest ) {
  TH1::SetDefaultSumw2;
  dest->Reset();

  TH3D* newErrors = (TH3D*)dest->Clone("newErrors");

  for ( int x = 1; x <= src->GetXaxis()->GetNbins(); ++x ) {
    for ( int y = 1; y <= src->GetYaxis()->GetNbins(); ++y ) {
      for ( int z = 1; z <= src->GetZaxis()->GetNbins(); ++z ) {
        double v = src->GetBinContent(x, y, z);
        double err = src->GetBinError(x, y, z);

        double x_cent = src->GetXaxis()->GetBinCenter(x);
        double y_cent = src->GetYaxis()->GetBinCenter(y);
        double z_cent = src->GetZaxis()->GetBinCenter(z);

        int n_x = dest->GetXaxis()->FindBin(x_cent);
        int n_y = dest->GetYaxis()->FindBin(y_cent);
        int n_z = dest->GetZaxis()->FindBin(z_cent);
        double update = dest->GetBinContent( n_x, n_y, n_z )+v;
        double error_update = newErrors->GetBinContent(n_x,n_y,n_z)+(err*err);
        dest->SetBinContent( n_x, n_y, n_z, update );
        newErrors->SetBinContent(n_x,n_y,n_z,error_update);
      }
    }
  }
  //error preparation                                       

  for ( int x = 1; x <= dest->GetXaxis()->GetNbins(); ++x ) {
    for ( int y = 1; y <= dest->GetYaxis()->GetNbins(); ++y ) {
      for ( int z = 1; z <= dest->GetZaxis()->GetNbins(); ++z ) {

	double err = newErrors->GetBinContent(x,y,z);
	double update = TMath::Sqrt(err);
	dest->SetBinError(x,y,z,update);
      }
    }
  }
}

