// macro to calculate the flux prediction at JUNO detector

Bool_t doDrawPred = false;

TCanvas *c1;

TSpline5* csSpline;
TF1*      fXsec;
TGraph*   csGraph;

enum detectors{
  kFar,
  nDet
};

TString tDet[nDet]={"JUNO"};
  
enum reactors{
  kYJC1,
  kYJC2,
  kYJC3,
  kYJC4,
  kYJC5,
  kYJC6,
  kTSC1,
  kTSC2,
  kTSC3,
  kTSC4,
  kDYB,
  kHZ,
  nRxt
};

TString tRxt[nRxt] = {"YJC1","YJC2","YJC3","YJC4","YJC5","YJC6","TSC1","TSC2","TSC3","TSC4","DYB","HZ"};

TH1D* hDet[nDet][nRxt];
TH1D* hTot[nDet];

TVectorD *predSigVecs[nDet];
TVectorD *predSigVec;

std::vector<Double_t> bin_edges;

enum FissileIso{
  kU5,
  kU8,
  kPu9,
  kPu1,
  nIsotopes
};

TH1D* hMCSpF;
TH2D* hcovmat;

TString isotope[nIsotopes]={"U235","U238","Pu239","Pu241"};

Double_t Baseline[nDet][nRxt] = {{5275000,5284000,5242000,5251000,5212000,5221000,5276000,5263000,5232000,5220000,21500000,26500000}}; // in centimeters

Double_t P_th[nRxt] = {2.9,2.9,2.9,2.9,2.9,2.9,4.6,4.6,4.6,4.6,17.4,17.4}; //GW
//Double_t Duty_cycle = 0.85; //on-off ratio
Double_t Duty_cycle = 1.0; //on-off ratio
Double_t alpha_k[nIsotopes] = {0.52, 0.09, 0.33, 0.06};
Double_t Mean_E_Fission[nIsotopes] = {201.92, 205.52, 209.99, 213.6};
Double_t Average_E_Fission = 0;
Double_t DetVolume[nDet] = {23227.82}; // m^3  (20kt)
Double_t LS_density = 0.859; //ton/m^3
Double_t protons_per_ton = 7.50e28; //(20% PXE and 80% n-dodecane)
Double_t efficiency = 0.73; //IBD selection

Bool_t doPerYear = true;

const Double_t day2second = 24.0*60.0*60.0;
const Double_t year2second = 365.25*24.0*60.0*60.0;
const Double_t GJ2MeV = 6.2415096471204e21;
const Double_t PI = TMath::Pi();
const Int_t nbin = 26; //MCSpF prediction

void GetPrediction(){

  //binning definition and nominal yield values
  const int nBins_DYB = nbin;
  // Generic anti-neutrino spectrum weighted by the IBD cross section [DYB reactor paper, Table 12].
  double eBins_DYB[nBins_DYB+1]     =   {1.800,  2.125,  2.375,  2.625,  2.875,  3.125,  3.375,  3.625,  3.875,  4.125,  4.375,  4.625,  4.875,  5.125,  5.375,  5.625,  5.875,  6.125,  6.375,  6.625,  6.875,  7.125,  7.375,  7.625,  7.875,  8.125, 12};
  double spec_DYB[nBins_DYB] = {344.19, 770.96, 1080.9, 1348.4, 1528.8, 1687.0, 1746.6, 1760.6, 1719.3, 1617.6, 1466.5, 1309.3, 1203.0, 1105.4, 976.50, 852.31, 713.19, 573.90, 463.54, 368.70, 274.56, 190.00, 132.08, 92.114, 56.689, 4.0214};

  for(int i=0;i<27;i++)
    bin_edges.push_back(eBins_DYB[i]);
  
  //hSpec = new TH1D("hcovmat","",nBins_DYB,eBins_DYB);
  hMCSpF = new TH1D("hMCSpF_DYB","hMCSpF_DYB",nBins_DYB,eBins_DYB);
  
  //Loading the cov. matrix

  //Covariance matrix of antineutrino spectrum.Unit: [cm^2/fission/MeV]^2 x 10e-92
  ifstream in;
  in.open("DYB_ReacCovMatrix.dat");
  double dummy;
  
  hcovmat = new TH2D("hcovmat","",nBins_DYB,eBins_DYB,nBins_DYB,eBins_DYB);
  //hcormat = new TH2D("hcormat","",nBins_DYB,eBins_DYB,nBins_DYB,eBins_DYB);
  
  for(int i=0; i<nBins_DYB; i++) {
    //hSpec->SetBinContent(i+1, spec_DYB[i] * 1e-46);
    Double_t BinWidth = hMCSpF->GetBinWidth(i+1);
    hMCSpF->SetBinContent(i+1, spec_DYB[i] * 1e-46 * BinWidth); //scale by bin width to get absolute value
    for(int j=0; j<nBins_DYB; j++) {
      in >> dummy;
      hcovmat->SetBinContent(i+1,j+1,dummy*1e-92);
    }
    //hSpec->SetBinError(i+1, TMath::Sqrt( hcovmat->GetBinContent(i+1,i+1) ) );
    hMCSpF->SetBinError(i+1, TMath::Sqrt( hcovmat->GetBinContent(i+1,i+1) ) );
  }
  in.close();
  
  //for(int i=0; i<nBins_DYB; i++) {
  //for(int j=0; j<nBins_DYB; j++) {
      //hcormat->SetBinContent(i+1,j+1, hcovmat->GetBinContent(i+1,j+1) / TMath::Sqrt(hcovmat->GetBinContent(i+1,i+1) * hcovmat->GetBinContent(j+1,j+1) ) );
  //}
    //cout<< hSpec->GetBinCenter(i+1)  <<" " <<TMath::Sqrt(hcovmat->GetBinContent(i+1,i+1) ) / hSpec->GetBinContent(i+1)<<endl;
  //}
  //hcovmat->Draw("colz");
  //hcormat->Draw("colz");
  //hSpec->Draw();

  cout<<"\nDYB MCSpF is:"<<hMCSpF->Integral()<<endl;
  //cout<<"DYB MCSpF is:"<<hMCSpF->Integral("width")<<endl;
  
  Average_E_Fission = 0;
  
  for(int j=0;j<nIsotopes;j++){
    Average_E_Fission += Mean_E_Fission[j] * alpha_k[j];
  }
  cout<<Average_E_Fission<<endl;

  Int_t nBinsTot;

  Int_t offset = 0.0;
  
  for(int iDet = 0; iDet<nDet; iDet++){

    hTot[iDet] = new TH1D(Form("h%s",tDet[iDet].Data()),Form("%s (%4.3f kt)",tDet[iDet].Data(),LS_density*DetVolume[iDet]*1e-3),bin_edges.size()-1,&bin_edges[0]);
    hTot[iDet]->SetLineWidth(3);
    hTot[iDet]->GetXaxis()->SetTitle("E_{#bar{#nu}_{e}} (MeV)");
    hTot[iDet]->GetYaxis()->SetTitle("Events / year");

    if(iDet == 0){ // trick to initialize nbinstot here
      nBinsTot = hTot[0]->GetNbinsX(); //total bins per histo
      predSigVec = new TVectorD(nBinsTot*nDet);
    }

    predSigVecs[iDet] = new TVectorD(nBinsTot);

    for(int iRxt = 0; iRxt<nRxt; iRxt++){

      hDet[iDet][iRxt] = new TH1D(Form("h%s%s",tDet[iDet].Data(),tRxt[iRxt].Data()),Form("%s - %s",tDet[iDet].Data(),tRxt[iRxt].Data()),bin_edges.size()-1,&bin_edges[0]);
      hDet[iDet][iRxt]->SetLineWidth(3);
      if(iRxt < kTSC1)
	hDet[iDet][iRxt]->SetLineColor(kBlue+iRxt-10);
      else
	hDet[iDet][iRxt]->SetLineColor(kRed+iRxt-kTSC1-10);

      hDet[iDet][iRxt]->Add(hMCSpF);

      Double_t Norm = protons_per_ton * LS_density * DetVolume[iDet]; // total number of targets
      Norm *= efficiency;
      Norm *= P_th[iRxt] * Duty_cycle * GJ2MeV; // thermal power
      Norm /= (4.0 * PI * Baseline[iDet][iRxt] * Baseline[iDet][iRxt] * Average_E_Fission);
      if(doPerYear)
	Norm *= year2second;
      else	
	Norm *= day2second;
      
      hDet[iDet][iRxt]->Scale(Norm);

      for(int ibin=0;ibin<nBinsTot;ibin++){
	//hDet[iDet][iRxt]->SetBinContent(ibin+1, hDet[iDet][iRxt]->GetBinContent(ibin+1) * hDet[iDet][iRxt]->GetBinWidth(ibin+1) );
	hDet[iDet][iRxt]->SetBinError(ibin+1,TMath::Sqrt(hDet[iDet][iRxt]->GetBinContent(ibin+1)));
      }
      hTot[iDet]->Add(hDet[iDet][iRxt]);
    }

    for(int ibin=0;ibin<nBinsTot;ibin++){
      (*predSigVecs[iDet])[ibin] = hTot[iDet]->GetBinContent(ibin+1);
    }

    predSigVec->SetSub(offset, *(predSigVecs[iDet]) );
    offset += predSigVecs[iDet]->GetNrows();

  }

  if(doDrawPred){
  
    c1 = new TCanvas("c1","c1");
    //c1->Divide(2,2);
    cout<<endl;
    Double_t max = hTot[0]->GetMaximum()*1.1;
    for(int iDet = 0; iDet<nDet; iDet++){
      //cout<<tDet[iDet].Data()<<" = "<<hTot[iDet]->Integral()<<endl;
      c1->cd(1+iDet);
      hTot[iDet]->GetYaxis()->SetRangeUser(0,max);
      hTot[iDet]->Draw("histe");
      for(int iR=0;iR<nRxt;iR++)
	hDet[iDet][iR]->Draw("samehiste");
      gPad->BuildLegend();

      cout<<tDet[iDet].Data()<<" : "<< hTot[iDet]->Integral()<< " ( ";
      for(int iR=0;iR<nRxt;iR++)
	cout<<hDet[iDet][iR]->Integral() << " + ";
      cout<<" ) "<<endl;
    }
    cout<<endl;
  }
}

void MakeCovMat(Int_t nSims = 1000){

  TStopwatch* watch = new TStopwatch();

  GetPrediction();
  /*
  Bool_t var_baseline  = false;
  Bool_t var_MeV_fissi = false;
  Bool_t var_alpha_k   = false;
  Bool_t var_Pth       = true;
  Bool_t var_flux_pred = false;
  */
  Bool_t var_baseline  = true;
  Bool_t var_MeV_fissi = true;
  Bool_t var_alpha_k   = true;
  Bool_t var_Pth       = true;
  Bool_t var_flux_pred = true;
  
  TString tUncConfig = "";

  if(var_baseline && var_MeV_fissi && var_alpha_k && var_Pth && var_flux_pred) tUncConfig = "tot";
  else if(var_baseline)  tUncConfig = "baseline";
  else if(var_MeV_fissi) tUncConfig = "MeVperFiss";
  else if(var_alpha_k)   tUncConfig = "alphaK";
  else if(var_Pth)       tUncConfig = "Pth";
  else if(var_flux_pred) tUncConfig = "FluxPred";
    
  TGraph mVariance;
  
  //TFile fMCSpF("FluxMatrixFromDCDB.root");  
  //TMatrixD* MCSpF_cov = (TMatrixD*)fMCSpF.Get("MCSpF_cov");
  TMatrixD* MCSpF_cov = new TMatrixD(nbin,nbin);
  for(int ibin=0;ibin<nbin;ibin++)
    for(int jbin=0;jbin<nbin;jbin++){
      Double_t BinSize = hMCSpF->GetBinWidth(ibin+1) * hMCSpF->GetBinWidth(jbin+1);

      (*MCSpF_cov)(ibin,jbin) = hcovmat->GetBinContent(ibin+1,jbin+1) * BinSize;
    }
  
  TH1D* hDet_temp[nDet][nRxt];
  TH1D* hTot_temp[nDet];
  TH1D* hMCSpF_temp;
  
  for(int iDet = 0; iDet<nDet; iDet++){
    hTot_temp[iDet] = new TH1D(Form("h%s_temp",tDet[iDet].Data()),Form("%s (%4.3f kt) : var",tDet[iDet].Data(),LS_density*DetVolume[iDet]*1e-3),bin_edges.size()-1,&bin_edges[0]);
    hTot_temp[iDet]->SetLineWidth(3);
    hTot_temp[iDet]->GetXaxis()->SetTitle("E_{#bar{#nu}_{e}} (MeV)");
    hTot_temp[iDet]->GetYaxis()->SetTitle("Events / year");
    for(int iRxt = 0; iRxt<nRxt; iRxt++){
      hDet_temp[iDet][iRxt] = new TH1D(Form("h%s%s_temp",tDet[iDet].Data(),tRxt[iRxt].Data()),Form("%s - %s : var",tDet[iDet].Data(),tRxt[iRxt].Data()),bin_edges.size()-1,&bin_edges[0]);
    }
  }

  hMCSpF_temp = new TH1D("hMCSpF_temp","hMCSpF_temp",bin_edges.size()-1,&bin_edges[0]);
  
  Double_t alpha_k_temp[nIsotopes][nRxt];
  Double_t P_th_temp[nRxt];
  Double_t Baseline_temp[nDet][nRxt];
  Double_t Mean_E_Fission_temp[nIsotopes];

  Double_t err_alpha_k[nIsotopes] = {3.3e-2, 6.5e-2, 4.0e-2, 11.0e-2}; // in %
  Double_t err_Pth = 0.5e-2; //in %
  Double_t err_baseline = 1.5; //in cm
  Double_t err_MeV_fiss[nIsotopes] = {0.46, 0.96, 0.60, 0.65}; //in MeV  ###### are there correlations ???

  Double_t rho_Pth_r = 0.0;
  Double_t rho_alpha_k_r = 1.0;

  Double_t rho_alpha_k_U5U8   =  0.33;
  Double_t rho_alpha_k_U5Pu9  = -0.93;
  Double_t rho_alpha_k_U5Pu1  = -0.91;
  Double_t rho_alpha_k_U8Pu9  = -0.36;
  Double_t rho_alpha_k_U8Pu1  = -0.37;
  Double_t rho_alpha_k_Pu9Pu1 =  0.91;

  Int_t Var2Randomize = nIsotopes*nRxt + nRxt + nDet*nRxt + nIsotopes + nbin; // alpha_k + P_th + baseline + MeV per Fission + MCSpF
  
  TMatrixD* var_cov = new TMatrixD(Var2Randomize,Var2Randomize);

  //alpha_k
  (*var_cov)(0,0) = TMath::Power(err_alpha_k[0] * alpha_k[0], 2);
  (*var_cov)(1,1) = TMath::Power(err_alpha_k[1] * alpha_k[1], 2);
  (*var_cov)(2,2) = TMath::Power(err_alpha_k[2] * alpha_k[2], 2);
  (*var_cov)(3,3) = TMath::Power(err_alpha_k[3] * alpha_k[3], 2);

  (*var_cov)(0,1) = err_alpha_k[0] * alpha_k[0] * err_alpha_k[1] * alpha_k[1] * rho_alpha_k_U5U8;
  (*var_cov)(1,0) = err_alpha_k[0] * alpha_k[0] * err_alpha_k[1] * alpha_k[1] * rho_alpha_k_U5U8;

  (*var_cov)(0,2) = err_alpha_k[0] * alpha_k[0] * err_alpha_k[2] * alpha_k[2] * rho_alpha_k_U5Pu9;
  (*var_cov)(2,0) = err_alpha_k[0] * alpha_k[0] * err_alpha_k[2] * alpha_k[2] * rho_alpha_k_U5Pu9;

  (*var_cov)(0,3) = err_alpha_k[0] * alpha_k[0] * err_alpha_k[3] * alpha_k[3] * rho_alpha_k_U5Pu1;
  (*var_cov)(3,0) = err_alpha_k[0] * alpha_k[0] * err_alpha_k[3] * alpha_k[3] * rho_alpha_k_U5Pu1;

  (*var_cov)(1,2) = err_alpha_k[1] * alpha_k[1] * err_alpha_k[2] * alpha_k[2] * rho_alpha_k_U8Pu9;
  (*var_cov)(2,1) = err_alpha_k[1] * alpha_k[1] * err_alpha_k[2] * alpha_k[2] * rho_alpha_k_U8Pu9;

  (*var_cov)(1,3) = err_alpha_k[1] * alpha_k[1] * err_alpha_k[3] * alpha_k[3] * rho_alpha_k_U8Pu1;
  (*var_cov)(3,1) = err_alpha_k[1] * alpha_k[1] * err_alpha_k[3] * alpha_k[3] * rho_alpha_k_U8Pu1;

  (*var_cov)(2,3) = err_alpha_k[2] * alpha_k[2] * err_alpha_k[3] * alpha_k[3] * rho_alpha_k_Pu9Pu1;
  (*var_cov)(3,2) = err_alpha_k[2] * alpha_k[2] * err_alpha_k[3] * alpha_k[3] * rho_alpha_k_Pu9Pu1;

  for(int i=0;i<nIsotopes;i++){
    for(int j=0;j<nIsotopes;j++){
      for(int iR=0;iR<nRxt;iR++)
	(*var_cov)(i+4*iR,j+4*iR) = (*var_cov)(i,j);
    }
  }
  
  //P_th  //correlation not implemented!!!
  for(int iR=0;iR<nRxt;iR++){
    (*var_cov)(nIsotopes*nRxt+iR, nIsotopes*nRxt+iR) = TMath::Power(err_Pth * P_th[iR], 2);
  }
  
  //baseline : need to add correlation btw baseline using same detector?
  for(int iR=0;iR<nRxt;iR++){
    (*var_cov)(nIsotopes*nRxt + nRxt + iR ,nIsotopes*nRxt + nRxt + iR ) = TMath::Power(err_baseline, 2);
  }
  
  //MeV per fission

  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt    , nIsotopes*nRxt + nRxt + nDet*nRxt    ) = TMath::Power(err_MeV_fiss[0], 2);
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + 1, nIsotopes*nRxt + nRxt + nDet*nRxt + 1) = TMath::Power(err_MeV_fiss[1], 2);
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + 2, nIsotopes*nRxt + nRxt + nDet*nRxt + 2) = TMath::Power(err_MeV_fiss[2], 2);
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + 3, nIsotopes*nRxt + nRxt + nDet*nRxt + 3) = TMath::Power(err_MeV_fiss[3], 2);
  
  for(Int_t ibin = 0; ibin < nbin; ++ibin)
    for(Int_t jbin = 0; jbin < nbin; ++jbin){
      (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + nIsotopes + ibin, nIsotopes*nRxt + nRxt + nDet*nRxt + nIsotopes + jbin) = (*MCSpF_cov)(ibin,jbin);
    }


  //var_cov->Print();
  
  TDecompChol cholesky1(*var_cov);
  Bool_t success = cholesky1.Decompose();

  if(success==0)
    cout<<"Choelsky Decomposition failed "<<endl;
  else if(success==1)
    cout<<"Choelsky Decomposition Success "<<endl;

  TMatrixD lower1 = cholesky1.GetU(); // get upper triangular matrix
  //lower1.Print();
  lower1.T(); // transpose to get lower triangular matrix

  TRandom* r1[Var2Randomize];

  const Int_t nBinsTot = hDet_temp[0][0]->GetNbinsX(); //any histo, since all are the same    
  TMatrixD* simMat = new TMatrixD(nBinsTot*nDet,nBinsTot*nDet);
  TVectorD *predVars[nDet];
  TVectorD *predVar        = new TVectorD(nBinsTot*nDet);

  // sims draws
  for(Int_t isim=0; isim<nSims; isim++){

    for(Int_t i=0;i<Var2Randomize;i++){r1[i] = new TRandom3(0.0);}
    
    for(int iDet = 0; iDet<nDet; iDet++){
      hTot_temp[iDet]->Reset();
      for(int iRxt = 0; iRxt<nRxt; iRxt++){
	hDet_temp[iDet][iRxt]->Reset();
      }
    }

    hMCSpF_temp->Reset();
    
    TVectorD* random1 = new TVectorD(Var2Randomize);
    for(Int_t i=0;i<Var2Randomize;i++){  (*random1)[i] = r1[i]->Gaus(0,1);    }

    TVectorD fluc1 = lower1 * (*random1);
    Double_t flucSum[nRxt] = {0.0};

    for(Int_t j=0;j<nRxt;j++){
      for(Int_t i=0;i<nIsotopes;i++){
	flucSum[j] += fluc1(i+j*nIsotopes) + alpha_k[i]; //get sum first, to normalize alpha_k later
      }
    }
    
    for(Int_t j=0;j<nRxt;j++){
      
      if(var_Pth){
	P_th_temp[j] = fluc1(j + nIsotopes*nRxt) + P_th[j];
      }
      else{
	P_th_temp[j] = P_th[j];
      }

      for(Int_t iDet=0;iDet<nDet;iDet++){
	if(var_baseline){
	  Baseline_temp[iDet][j] = fluc1(iDet + nIsotopes*nRxt + nRxt) + Baseline[iDet][j];
	}
	else{
	  Baseline_temp[iDet][j] = Baseline[iDet][j];
	}
      }
      
      for(Int_t i=0;i<nIsotopes;i++){

	if(var_MeV_fissi)
	  Mean_E_Fission_temp[i] = fluc1(i + nIsotopes*nRxt + nRxt + nRxt*nDet) + Mean_E_Fission[i];
	else
	  Mean_E_Fission_temp[i] = Mean_E_Fission[i];
	
	if(var_alpha_k){
	  if( rho_alpha_k_r == 0)
	    alpha_k_temp[i][j] = (fluc1(i+nIsotopes*j)+alpha_k[i]) / flucSum[j];
	  else{
	    alpha_k_temp[i][0] = (fluc1(i+nIsotopes*0)+alpha_k[i]) / flucSum[0];
	    alpha_k_temp[i][j] = alpha_k_temp[i][0];
	  }
	}
	else
	  alpha_k_temp[i][j] = alpha_k[i];

	//cout<<(fluc1(i+nIsotopes*j)+alpha_k[i]) / flucSum[j] <<" "<<alpha_k[i]<<" ";
      }
      //cout<<" | ";
    }
    //cout<<endl;

    for(Int_t ibin=0;ibin<nbin;ibin++){

      Double_t binContent;
      if(var_flux_pred)
	binContent = fluc1(nIsotopes + nIsotopes*nRxt + nRxt + nRxt*nDet + ibin) + hMCSpF->GetBinContent(ibin+1);
      else
	binContent = hMCSpF->GetBinContent(ibin+1);

      hMCSpF_temp->SetBinContent(ibin+1, binContent);
	
    }

    Double_t Average_E_Fission_temp[nRxt] = {0.0};
    for(int i=0;i<nIsotopes;i++){
      for(Int_t j=0;j<nRxt;j++){
	Average_E_Fission_temp[j] += Mean_E_Fission_temp[i] * alpha_k_temp[i][j];
      }
    }
    //cout<<Average_E_Fission<<" "<<Average_E_Fission_temp[0]<<" "<<Average_E_Fission_temp[1]<<endl;

    Double_t offset = 0.0;
    
    for(int iDet = 0; iDet<nDet; iDet++){
      
      predVars[iDet] = new TVectorD(nBinsTot);
      
      for(int iRxt = 0; iRxt<nRxt; iRxt++){

	hDet_temp[iDet][iRxt]->Add(hMCSpF_temp);

	Double_t Norm = protons_per_ton * LS_density * DetVolume[iDet]; // total number of targets
	Norm *= efficiency;
	Norm *= P_th_temp[iRxt] * Duty_cycle * GJ2MeV; // thermal power
	Norm /= (4.0 * PI * Baseline_temp[iDet][iRxt] * Baseline_temp[iDet][iRxt] * Average_E_Fission_temp[iRxt]);
	if(doPerYear)
	  Norm *= year2second;
	else
	  Norm *= day2second;

	hDet_temp[iDet][iRxt]->Scale(Norm);
	
	for(int ibin=0;ibin<nBinsTot;ibin++){
	  //hDet_temp[iDet][iRxt]->SetBinContent(ibin+1, hDet_temp[iDet][iRxt]->GetBinContent(ibin+1) * hDet_temp[iDet][iRxt]->GetBinWidth(ibin+1) );
	  hDet_temp[iDet][iRxt]->SetBinError(ibin+1,TMath::Sqrt(hDet_temp[iDet][iRxt]->GetBinContent(ibin+1)));
	}
	hTot_temp[iDet]->Add(hDet_temp[iDet][iRxt]);
      }

      for(int ibin=0;ibin<nBinsTot;ibin++){
	(*predVars[iDet])[ibin] = hTot_temp[iDet]->GetBinContent(ibin+1);
      }

      predVar->SetSub(offset, *(predVars[iDet]) );
      offset += predVars[iDet]->GetNrows();

    }
    
    //fill matrix
    for(Int_t ibin = 0; ibin < nBinsTot*nDet; ++ibin)
      for(Int_t jbin = 0; jbin < nBinsTot*nDet; ++jbin){

	(*simMat)(ibin,jbin)      += ((*predVar)[ibin]      - (*predSigVec)[ibin]) * ((*predVar)[jbin]      - (*predSigVec)[jbin] );

      }
    
    delete random1;
    random1 = NULL;
    
    for(Int_t i=0;i<Var2Randomize;i++){
      delete r1[i];
      r1[i] = NULL;
    }
    for(Int_t i=0;i<nDet;i++){
      delete predVars[i];
      predVars[i] = NULL;
    }

    Double_t CovSum[nDet] = {0.0};
    for(Int_t i=0;i<nDet;i++){
      for(int jbin=0;jbin<nbin;jbin++){
	for(int ibin=0;ibin<nbin;ibin++){
	  CovSum[i] += (*simMat)(ibin+nbin*i,jbin+nbin*i);
	}
      }
    }
    
    //mVariance.SetPoint(isim, (Double_t)isim, 100.0*(TMath::Sqrt( simMat->Sum() / ((Double_t)isim+1.)) /( predSigVec->Sum() )));
    mVariance.SetPoint(isim, (Double_t)isim, 100.0*(TMath::Sqrt( CovSum[0] / ((Double_t)isim+1.)) /( hTot_temp[0]->Integral() )));
    
  }
  
  (*simMat) *= 1./nSims;
  
  Double_t CovSum[nDet] = {0.0};
  for(Int_t i=0;i<nDet;i++){
    for(int jbin=0;jbin<nbin;jbin++){
      for(int ibin=0;ibin<nbin;ibin++){
	CovSum[i] += (*simMat)(ibin+nbin*i,jbin+nbin*i);
      }
    }
  }

  
  cout<<"bin fractional error in %:"<<endl;
  for(int ibin=0;ibin<nbin;ibin++){
    cout<<ibin<<" ";
    for(Int_t i=0;i<nDet;i++)
      cout<< 100.0*TMath::Sqrt((*simMat)(ibin+nbin*i,ibin+nbin*i))/hTot_temp[i]->GetBinContent(ibin+1)<<" ";
    cout<<endl;
  }

  cout<< "\n Norm unc = " ;
  for(Int_t i=0;i<nDet;i++)
    cout<<100.0 * TMath::Sqrt(CovSum[i]) / hTot_temp[i]->Integral() << "% ";
  cout<<endl;
  
  TMatrixD* fracMat = new TMatrixD( *simMat );
  TMatrixD* corrMat = new TMatrixD( *simMat );

   for(Int_t ibin = 0; ibin < nBinsTot*nDet; ++ibin)
     for(Int_t jbin = 0; jbin < nBinsTot*nDet; ++jbin)
       {
	 if( (*predSigVec)[ibin] != 0. && (*predSigVec)[jbin] != 0. ){

           if( (*simMat)(ibin,ibin)!=0.0 &&  (*simMat)(jbin,jbin)!=0.0 ){
             (*corrMat)(ibin,jbin) /= TMath::Sqrt( (*simMat)(ibin,ibin)) * TMath::Sqrt((*simMat)(jbin,jbin));
           }

           (*fracMat)(ibin,jbin) /= (*predSigVec)[ibin]*(*predSigVec)[jbin];
	   
         }
       }

   TFile* outfile = new TFile(Form("../../data/FluxPredMat_%dsims_DayaBay_%s.root",nSims,tUncConfig.Data()),"RECREATE");

   simMat->Write("flux_cov");
   fracMat->Write("flux_frac_cov");
   corrMat->Write("flux_corr");
   mVariance.Write("variance_nsims");

   for(int iDet = 0; iDet<nDet; iDet++){
     hTot[iDet]->Write();
     for(int iR = 0; iR<nRxt; iR++){
       hDet[iDet][iR]->Write();
       hDet[iDet][iR]->Write();
     }
   }

   outfile->Close();

   // stopwatch output
   watch->Stop();
   watch->Print();

   delete watch;
   watch = NULL;

}
