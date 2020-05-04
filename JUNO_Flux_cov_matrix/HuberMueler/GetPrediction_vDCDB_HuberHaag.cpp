// macro to calculate the flux prediction at ND and FD of LiquidO-Chooz

Bool_t doDrawPred = false;
Bool_t doDrawFlux = false;

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

TH1D* hMCSpF[nIsotopes];

TString isotope[nIsotopes]={"U235","U238","Pu239","Pu241"};

Double_t Baseline[nDet][nRxt] = {{5275000,5284000,5242000,5251000,5212000,5221000,5276000,5263000,5232000,5220000,21500000,26500000}}; // in centimeters

Double_t P_th[nRxt] = {2.9,2.9,2.9,2.9,2.9,2.9,4.6,4.6,4.6,4.6,17.4,17.4}; //GW
//Double_t P_th[nRxt] = {2.9,0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001}; //GW
//Double_t Duty_cycle = 0.85; //on-off ratio
Double_t Duty_cycle = 1.0; //on-off ratio
Double_t alpha_k[nIsotopes] = {0.52, 0.09, 0.33, 0.06};
Double_t Mean_E_Fission[nIsotopes] = {201.92, 205.52, 209.99, 213.6};
Double_t DetVolume[nDet] = {23227.82}; // m^3  (20kt)
Double_t LS_density = 0.859; //ton/m^3
Double_t protons_per_ton = 7.50e28; //(LAB, yellow book?)
Double_t efficiency = 0.73; //IBD selection

Bool_t doPerYear = true;

const Double_t day2second = 24.0*60.0*60.0;
const Double_t year2second = 365.25*24.0*60.0*60.0;
const Double_t GJ2MeV = 6.2415096471204e21;
const Double_t PI = TMath::Pi();
const Int_t nbin = 25; //MCSpF prediction

Double_t integralCS[nbin];

Double_t myCS(Double_t *x, Double_t *par){

  Double_t xx = x[0];

  return csSpline->Eval(xx);
}

void GetCrossSection(){

  TString filePath("");
  filePath.Append("gm_xsec.dat");

  std::ifstream runFile;
  runFile.open( filePath.Data() );
  if( !runFile.is_open() )
    { cout << "Nope." << endl; exit(1); }

  csGraph = new TGraph(1);

  //Prefactor, with MAMBO-II value
  const Double_t K = 0.961e-43;

  std::string line;
  Int_t ipoint = 0;
  while( !runFile.eof() )
    {
      getline( runFile, line );
      //values separated by ','
      string e;
      string xs;

      e = line.substr(0,line.find(","));
      xs = line.substr(line.find(",")+1);

      csGraph->SetPoint(ipoint, atof(e.c_str()), K*atof(xs.c_str()));
      ipoint++;
    }

  csSpline = new TSpline5("csSpline",csGraph);

  fXsec = new TF1("Xsec",myCS,0.5,25.0,6);

}

void GetPrediction(){

  GetCrossSection();

  TString filename[nIsotopes];
  ifstream rIn[nIsotopes];
  
  for(int i=0;i<nIsotopes;i++){
    filename[i]= isotope[i];
    filename[i]+="_v1";
    filename[i]+=".txt";
    rIn[i].open(filename[i].Data());
  }

  Double_t lowEne[nIsotopes][nbin];
  Double_t uppEne[nIsotopes][nbin];
  Double_t EneCV[nIsotopes][nbin];
  Double_t nu_per_fis[nIsotopes][nbin];
  Double_t d_tot[nIsotopes][nbin];
  Double_t d_corr[nIsotopes][nbin];
  Double_t d_uncorr[nIsotopes][nbin];
  Double_t d_ene[nIsotopes][nbin]={0.0};
  Double_t Err[nIsotopes][nbin];

  for(int j=0;j<nIsotopes;j++){  
    for(int i=0;i<nbin;i++){
      rIn[j]>>lowEne[j][i]>>uppEne[j][i]>>nu_per_fis[j][i]>>d_tot[j][i]>>d_corr[j][i]>>d_uncorr[j][i];

      EneCV[j][i] = lowEne[j][i] + (uppEne[j][i]-lowEne[j][i])/2.0;
      Err[j][i] = d_tot[j][i]*nu_per_fis[j][i]/100.0;

    }
  }

  for(int i=0;i<nbin;i++){
    bin_edges.push_back(lowEne[0][i]);
  }
  bin_edges.push_back(uppEne[0][nbin-1]);

  for(int j=0;j<nIsotopes;j++){

    hMCSpF[j] = new TH1D(Form("h%s",isotope[j].Data()),isotope[j].Data(),bin_edges.size()-1,&bin_edges[0]);

    for(int i=0;i<nbin;i++){

      integralCS[i] = fXsec->Integral(lowEne[0][i],uppEne[0][i]);
      integralCS[i] /= (uppEne[0][i] - lowEne[0][i]);
    
      hMCSpF[j]->SetBinContent(i+1,nu_per_fis[j][i] * integralCS[i]);
      hMCSpF[j]->SetBinError(i+1,Err[j][i] * integralCS[i]);
      //hMCSpF[j]->SetBinContent(i+1,nu_per_fis[j][i] );
      //hMCSpF[j]->SetBinError(i+1,Err[j][i]);
    
    }
    cout<<isotope[j].Data()<<" : "<<hMCSpF[j]->Integral("width")<<" nu*cm^2/fission"<<endl;
  }


  if(doDrawFlux){
    c1 = new TCanvas("c1","c1");

    hMCSpF[kU8]->SetLineWidth(3);
    hMCSpF[kU5]->SetLineWidth(3);
    hMCSpF[kPu9]->SetLineWidth(3);
    hMCSpF[kPu1]->SetLineWidth(3);
    
    hMCSpF[kU8]->SetLineColor(kRed-7);
    hMCSpF[kU5]->SetLineColor(kRed);
    hMCSpF[kPu9]->SetLineColor(kBlue);
    hMCSpF[kPu1]->SetLineColor(kBlue-7);

    hMCSpF[kU8]->SetFillColor(kRed-7);
    hMCSpF[kU5]->SetFillColor(kRed);
    hMCSpF[kPu9]->SetFillColor(kBlue);
    hMCSpF[kPu1]->SetFillColor(kBlue-7);

    hMCSpF[kU8]->SetMarkerColor(kRed-7);
    hMCSpF[kU5]->SetMarkerColor(kRed);
    hMCSpF[kPu9]->SetMarkerColor(kBlue);
    hMCSpF[kPu1]->SetMarkerColor(kBlue-7);

    hMCSpF[kU8]->Draw("E3");
    hMCSpF[kU5]->Draw("same E3");
    hMCSpF[kPu9]->Draw("same E3");
    hMCSpF[kPu1]->Draw("same E3");
    hMCSpF[kU8]->GetXaxis()->SetTitle("#bar{#nu}_{e} Energy (MeV)");
    hMCSpF[kU8]->GetYaxis()->SetTitle("Flux (#bar{#nu}_{e} #times cm^{2} / fission)");
    hMCSpF[kU8]->GetYaxis()->SetTitleOffset(1.5);
    gPad->BuildLegend();
    hMCSpF[kU8]->SetTitle("");
  }

  Double_t Average_E_Fission = 0;
  for(int j=0;j<nIsotopes;j++){
    Average_E_Fission += Mean_E_Fission[j] * alpha_k[j];
  }

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

      for(int j=0;j<nIsotopes;j++){
	
	hDet[iDet][iRxt]->Add(hMCSpF[j], alpha_k[j]);

      }

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
	hDet[iDet][iRxt]->SetBinContent(ibin+1, hDet[iDet][iRxt]->GetBinContent(ibin+1) * hDet[iDet][iRxt]->GetBinWidth(ibin+1) );
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
  
  Bool_t var_baseline  = true;
  Bool_t var_MeV_fissi = true;
  Bool_t var_alpha_k   = true;
  Bool_t var_Pth       = true;
  Bool_t var_flux_pred = true;
  /*
  Bool_t var_baseline  = false;
  Bool_t var_MeV_fissi = false;
  Bool_t var_alpha_k   = false;
  Bool_t var_Pth       = true;
  Bool_t var_flux_pred = false;
  */
  TString tUncConfig = "";

  if(var_baseline && var_MeV_fissi && var_alpha_k && var_Pth && var_flux_pred) tUncConfig = "tot";
  else if(var_baseline)  tUncConfig = "baseline";
  else if(var_MeV_fissi) tUncConfig = "MeVperFiss";
  else if(var_alpha_k)   tUncConfig = "alphaK";
  else if(var_Pth)       tUncConfig = "Pth";
  else if(var_flux_pred) tUncConfig = "FluxPred";

  TGraph mVariance;
  
  TFile fMCSpF("FluxMatrixFromDCDB.root");
  TMatrixD* MCSpF_cov = (TMatrixD*)fMCSpF.Get("MCSpF_cov");
  
  TH1D* hDet_temp[nDet][nRxt];
  TH1D* hTot_temp[nDet];
  TH1D* hMCSpF_temp[nIsotopes];
  
  for(int iDet = 0; iDet<nDet; iDet++){
    hTot_temp[iDet] = new TH1D(Form("h%s_temp",tDet[iDet].Data()),Form("%s (%4.3f kt) : var",tDet[iDet].Data(),LS_density*DetVolume[iDet]*1e-3),bin_edges.size()-1,&bin_edges[0]);
    hTot_temp[iDet]->SetLineWidth(3);
    hTot_temp[iDet]->GetXaxis()->SetTitle("E_{#bar{#nu}_{e}} (MeV)");
    if(doPerYear) hTot_temp[iDet]->GetYaxis()->SetTitle("Events / year");
    else hTot_temp[iDet]->GetYaxis()->SetTitle("Events / day");
    for(int iRxt = 0; iRxt<nRxt; iRxt++){
      hDet_temp[iDet][iRxt] = new TH1D(Form("h%s%s_temp",tDet[iDet].Data(),tRxt[iRxt].Data()),Form("%s - %s : var",tDet[iDet].Data(),tRxt[iRxt].Data()),bin_edges.size()-1,&bin_edges[0]);
    }
  }
  for(int iIso = 0; iIso<nIsotopes; iIso++){
    hMCSpF_temp[iIso] = new TH1D(Form("h%s_temp",isotope[iIso].Data()),isotope[iIso].Data(),bin_edges.size()-1,&bin_edges[0]);
  }
  
  Double_t alpha_k_temp[nIsotopes][nRxt];
  Double_t P_th_temp[nRxt];
  Double_t Baseline_temp[nDet][nRxt];
  Double_t Mean_E_Fission_temp[nIsotopes];

  Double_t err_alpha_k[nIsotopes] = {3.3e-2, 6.5e-2, 4.0e-2, 11.0e-2}; // in %
  Double_t err_Pth = 0.5e-2; //in %
  //Double_t err_Pth = 0.8e-2; //in %
  //Double_t err_Pth = 2.0e-2; //in %
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

  Int_t Var2Randomize = nIsotopes*nRxt + nRxt + nDet*nRxt + nIsotopes + nIsotopes*nbin; // alpha_k + P_th + baseline + MeV per Fission + MCSpF
  //Int_t Var2Randomize = nIsotopes*nRxt + nRxt + nDet*nRxt + nIsotopes; // alpha_k + P_th + baseline + MeV per Fission
  
  TMatrixD* var_cov        = new TMatrixD(Var2Randomize,Var2Randomize);
  TMatrixD* var_cov_approx = new TMatrixD(Var2Randomize,Var2Randomize);

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
  
  //P_th
  for(int iR=0;iR<nRxt;iR++){
    for(int jR=0;jR<nRxt;jR++){
      (*var_cov)(nIsotopes*nRxt+iR, nIsotopes*nRxt+jR) = err_Pth * P_th[iR] * err_Pth * P_th[jR] * rho_Pth_r;
    }
    (*var_cov)(nIsotopes*nRxt+iR, nIsotopes*nRxt+iR) = TMath::Power(err_Pth * P_th[iR], 2);
  }
  
  //(*var_cov)(nIsotopes*nRxt+1, nIsotopes*nRxt+1) = TMath::Power(err_Pth * P_th[1], 2);
  //(*var_cov)(nIsotopes*nRxt+2, nIsotopes*nRxt+2) = TMath::Power(err_Pth * P_th[2], 2);
  //(*var_cov)(nIsotopes*nRxt+3, nIsotopes*nRxt+3) = TMath::Power(err_Pth * P_th[3], 2);
  //(*var_cov)(nIsotopes*nRxt+4, nIsotopes*nRxt+4) = TMath::Power(err_Pth * P_th[4], 2);
  //(*var_cov)(nIsotopes*nRxt+5, nIsotopes*nRxt+5) = TMath::Power(err_Pth * P_th[5], 2);
  //(*var_cov)(nIsotopes*nRxt+6, nIsotopes*nRxt+6) = TMath::Power(err_Pth * P_th[6], 2);
  //(*var_cov)(nIsotopes*nRxt+7, nIsotopes*nRxt+7) = TMath::Power(err_Pth * P_th[7], 2);
  //(*var_cov)(nIsotopes*nRxt+8, nIsotopes*nRxt+8) = TMath::Power(err_Pth * P_th[8], 2);
  //(*var_cov)(nIsotopes*nRxt+9, nIsotopes*nRxt+9) = TMath::Power(err_Pth * P_th[9], 2);

  //(*var_cov)(nIsotopes*nRxt  , nIsotopes*nRxt+1) = err_Pth * P_th[0] * err_Pth * P_th[1] * rho_Pth_r;
  //(*var_cov)(nIsotopes*nRxt+1, nIsotopes*nRxt  ) = err_Pth * P_th[0] * err_Pth * P_th[1] * rho_Pth_r;
  
  //baseline : need to add correlation btw baseline using same detector?
  for(int iR=0;iR<nRxt;iR++){
    (*var_cov)(nIsotopes*nRxt + nRxt + iR ,nIsotopes*nRxt + nRxt + iR ) = TMath::Power(err_baseline, 2);
  }
  //(*var_cov)(nIsotopes*nRxt + nRxt + 1,nIsotopes*nRxt + nRxt + 1) = TMath::Power(err_baseline, 2);

  //MeV per fission
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt    , nIsotopes*nRxt + nRxt + nDet*nRxt    ) = TMath::Power(err_MeV_fiss[0], 2);
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + 1, nIsotopes*nRxt + nRxt + nDet*nRxt + 1) = TMath::Power(err_MeV_fiss[1], 2);
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + 2, nIsotopes*nRxt + nRxt + nDet*nRxt + 2) = TMath::Power(err_MeV_fiss[2], 2);
  (*var_cov)(nIsotopes*nRxt + nRxt + nDet*nRxt + 3, nIsotopes*nRxt + nRxt + nDet*nRxt + 3) = TMath::Power(err_MeV_fiss[3], 2);
  
  for(Int_t ibin = 0; ibin < nbin*nIsotopes; ++ibin)
    for(Int_t jbin = 0; jbin < nbin*nIsotopes; ++jbin){
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
    //for(Int_t i=0;i<Var2Randomize;i++){r1[i] = new TRandom3(555+i*23);}
    
    for(int iDet = 0; iDet<nDet; iDet++){
      hTot_temp[iDet]->Reset();
      for(int iRxt = 0; iRxt<nRxt; iRxt++){
	hDet_temp[iDet][iRxt]->Reset();
      }
    }
    for(Int_t i=0;i<nIsotopes;i++){
      hMCSpF_temp[i]->Reset();
    }
    
    TVectorD* random1 = new TVectorD(Var2Randomize);
    for(Int_t i=0;i<Var2Randomize;i++){  (*random1)[i] = r1[i]->Gaus(0,1);    }

    TVectorD fluc1 = lower1 * (*random1);
    TVectorD fluc2 = lower1 * (*random1);
    TVectorD CV1(Var2Randomize);
    Double_t flucSum[nRxt] = {0.0};

    for(Int_t j=0;j<nRxt;j++){
      for(Int_t i=0;i<nIsotopes;i++){
	flucSum[j] += fluc1(i+j*nIsotopes) + alpha_k[i]; //get sum first, to normalize alpha_k later
	CV1(i+j*nIsotopes) = alpha_k[i];
      }
    }
    
    for(Int_t j=0;j<nRxt;j++){
      
      if(var_Pth){
	P_th_temp[j] = fluc1(j + nIsotopes*nRxt) + P_th[j];
	fluc2(j + nIsotopes*nRxt) = P_th_temp[j];
	CV1(j + nIsotopes*nRxt)   = P_th[j];
      }
      else{
	P_th_temp[j] = P_th[j];
      }

      for(Int_t iDet=0;iDet<nDet;iDet++){
	if(var_baseline){
	  Baseline_temp[iDet][j] = fluc1(iDet + nIsotopes*nRxt + nRxt) + Baseline[iDet][j];
	  fluc2(iDet + nIsotopes*nRxt + nRxt) = Baseline_temp[iDet][j];
	  CV1(iDet + nIsotopes*nRxt + nRxt) = Baseline[iDet][j];
	}
	else{
	  Baseline_temp[iDet][j] = Baseline[iDet][j];
	  //cout<<Baseline_temp[iDet][j]<<endl;
	}
      }
      
      for(Int_t i=0;i<nIsotopes;i++){

	if(var_MeV_fissi){
	  Mean_E_Fission_temp[i] = fluc1(i + nIsotopes*nRxt + nRxt + nRxt*nDet) + Mean_E_Fission[i];
	  fluc2(i + nIsotopes*nRxt + nRxt + nRxt*nDet) = Mean_E_Fission_temp[i];
	  CV1(i + nIsotopes*nRxt + nRxt + nRxt*nDet) = Mean_E_Fission[i];
	}
	else
	  Mean_E_Fission_temp[i] = Mean_E_Fission[i];
	
	if(var_alpha_k){
	  if( rho_alpha_k_r == 0){
	    alpha_k_temp[i][j] = (fluc1(i+nIsotopes*j)+alpha_k[i]) / flucSum[j];
	    fluc2(i+j*nIsotopes) = alpha_k_temp[i][j];
	  }
	  else{
	    alpha_k_temp[i][0] = (fluc1(i+nIsotopes*0)+alpha_k[i]) / flucSum[0];
	    alpha_k_temp[i][j] = alpha_k_temp[i][0];
	    fluc2(i+0*nIsotopes) = alpha_k_temp[i][j];
	  }
	}
	else
	  alpha_k_temp[i][j] = alpha_k[i];

	//cout<<(fluc1(i+nIsotopes*j)+alpha_k[i]) / flucSum[j] <<" ";
	//cout<<alpha_k_temp[i][0]<<" "<<alpha_k[i]<<" ";
      }
      //cout<<endl;
    }

    for(Int_t iIso=0;iIso<nIsotopes;iIso++){
      for(Int_t ibin=0;ibin<nbin;ibin++){

	Double_t binContent;
	if(var_flux_pred){
	  binContent = fluc1(nIsotopes + nIsotopes*nRxt + nRxt + nRxt*nDet + ibin + nbin*iIso)*integralCS[ibin] + hMCSpF[iIso]->GetBinContent(ibin+1);
	  //binContent = fluc1(nIsotopes + nIsotopes*nRxt + nRxt + nRxt*nDet + ibin + nbin*iIso) + hMCSpF[iIso]->GetBinContent(ibin+1);
	  fluc2(nIsotopes + nIsotopes*nRxt + nRxt + nRxt*nDet + ibin + nbin*iIso) = binContent;
	  CV1(nIsotopes + nIsotopes*nRxt + nRxt + nRxt*nDet + ibin + nbin*iIso) = hMCSpF[iIso]->GetBinContent(ibin+1);
	}
	else
	  binContent = hMCSpF[iIso]->GetBinContent(ibin+1);

	hMCSpF_temp[iIso]->SetBinContent(ibin+1, binContent);
	
      }
    }

    Double_t Average_E_Fission[nRxt] = {0.0};
    for(int i=0;i<nIsotopes;i++){
      for(Int_t j=0;j<nRxt;j++){
	Average_E_Fission[j] += Mean_E_Fission_temp[i] * alpha_k_temp[i][j];
      }
    }
    //cout<<Average_E_Fission[2]<<endl;

    Double_t offset = 0.0;
    
    for(int iDet = 0; iDet<nDet; iDet++){
      
      predVars[iDet] = new TVectorD(nBinsTot);
      
      for(int iRxt = 0; iRxt<nRxt; iRxt++){

	for(int j=0;j<nIsotopes;j++){
	  
	  hDet_temp[iDet][iRxt]->Add(hMCSpF_temp[j], alpha_k_temp[j][iRxt]);
	}

	Double_t Norm = protons_per_ton * LS_density * DetVolume[iDet]; // total number of targets
	Norm *= efficiency;
	Norm *= P_th_temp[iRxt] * Duty_cycle * GJ2MeV; // thermal power
	Norm /= (4.0 * PI * Baseline_temp[iDet][iRxt] * Baseline_temp[iDet][iRxt] * Average_E_Fission[iRxt]);
	
	if(doPerYear)
	  Norm *= year2second;
	else
	  Norm *= day2second;
		
	hDet_temp[iDet][iRxt]->Scale(Norm);
	
	for(int ibin=0;ibin<nBinsTot;ibin++){
	  hDet_temp[iDet][iRxt]->SetBinContent(ibin+1, hDet_temp[iDet][iRxt]->GetBinContent(ibin+1) * hDet_temp[iDet][iRxt]->GetBinWidth(ibin+1) );
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
    //cout<<hTot_temp[0]->Integral()<<" "<<hTot[0]->Integral()<<endl;
    
    //fill matrix
    for(Int_t ibin = 0; ibin < nBinsTot*nDet; ++ibin)
      for(Int_t jbin = 0; jbin < nBinsTot*nDet; ++jbin){

	(*simMat)(ibin,jbin)      += ((*predVar)[ibin]      - (*predSigVec)[ibin]) * ((*predVar)[jbin]      - (*predSigVec)[jbin] );

      }

    for(Int_t ibin = 0; ibin < Var2Randomize; ++ibin)
      for(Int_t jbin = 0; jbin < Var2Randomize; ++jbin){
	
	(*var_cov_approx)(ibin,jbin) += (fluc2(ibin) - CV1(ibin) ) * (fluc2(jbin) - CV1(jbin));

      }
    //cout<<fluc2(25)<<" "<<CV1(25)<<endl;
    //cout<<fluc2(10)<<" "<<CV1(10)<<endl;

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
  
  (*simMat)         *= 1./nSims;
  (*var_cov_approx) *= 1./nSims;
  
  Double_t CovSum[nDet] = {0.0};
  for(Int_t i=0;i<nDet;i++){
    for(int jbin=0;jbin<nbin;jbin++){
      for(int ibin=0;ibin<nbin;ibin++){
	CovSum[i] += (*simMat)(ibin+nbin*i,jbin+nbin*i);
      }
    }
  }

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

   TFile* outfile = new TFile(Form("FluxPredMat_%dsims_HuberHaag_10baselines_%s.root",nSims,tUncConfig.Data()),"RECREATE");

   simMat->Write("flux_cov");
   fracMat->Write("flux_frac_cov");
   corrMat->Write("flux_corr");
   mVariance.Write("variance_nsims");
   var_cov_approx->Write("input_cov_approx");

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
