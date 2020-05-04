//macro to generate the Covariance matrix of Energy scale (shape only and rate+shape)

Double_t survProb3nu(const Double_t E, const Double_t L,
		     const Double_t th13, const Double_t th12,
		     const Double_t dm31, const Double_t dm21)
{
  Double_t dm32 = dm31 - dm21;
  
  Double_t Psurv = 1.;
  
  Psurv -= TMath::Power( TMath::Sin(2*th13) ,2.0) * (
						     TMath::Power( TMath::Cos(th12)   ,2.0) * TMath::Power( TMath::Sin(1.267*dm31*L/E) ,2.0) +
						     TMath::Power( TMath::Sin(th12)   ,2.0) * TMath::Power( TMath::Sin(1.267*dm32*L/E) ,2.0)   );
  
  Psurv -= TMath::Power( TMath::Cos(th13)  , 4.0) *  TMath::Power( TMath::Sin(2*th12) ,2.0) * TMath::Power( TMath::Sin(1.267*dm21*L/E) ,2.0) ;

  return Psurv;
}

void MakeEscaleCovMat_DC(Int_t SIMS=1000, TString inputFile = "../data/events_SPMT_1e+06evt_2019_test_DYB.root", Bool_t isDYB = true, Bool_t useOsc = false, Bool_t isProspect = false){

  const Double_t PE2MeV = 52.57;
  
  std::vector<Float_t> vFarEvis;
  std::vector<Float_t> vFarWeight;
  
  TFile* fFar  = new TFile(inputFile.Data());

  TTree* tFar  = (TTree*)fFar->Get("FinalFitIBDTree");

  TGraph mVariance;
  TGraph mErrorBin12FD;

  Int_t FarEvents = tFar->GetEntries();

  Double_t EvisF,Eth,Dist;
  tFar->SetBranchAddress("myPromptEvisID",&EvisF);
  tFar->SetBranchAddress("myNeutrinoEnergy_Th",&Eth);
  tFar->SetBranchAddress("myNeutrinoDistance_Th",&Dist);
  
  std::vector<Double_t> bin_edges;
  Double_t E = 40.0;

  for(int ibin=0; ibin<100; ibin++) // 40 < E < 440 PE
    {
      bin_edges.push_back(E);
      E += 4; // 4 PE bins
    }
  bin_edges.push_back(E);
  
  Double_t Dm21 = 7.53e-5;
  //Double_t Dm21 = 7.53e-5 + 0.18e-5 * 1.0; //1sigma shift
  Double_t th13 = 0.14846;  //DYB 1809.02261
  Double_t Dm31 = 2.47e-3;
  //Double_t tan12 = 0.436;
  Double_t tan12 = 0.436 + 0.029 * 1.0; //1sigma shift
  Double_t th12 =  TMath::ATan( TMath::Sqrt(tan12) );

  if(!useOsc){
    th12 = 0.0; Dm21 = 0.0;
  }
  
  TH1D* hFar  = new TH1D("hFar ","",bin_edges.size()-1,&bin_edges[0]);

  for(Int_t evt = 0; evt<FarEvents; evt++){
    tFar->GetEntry(evt);
    Double_t weight = 1;

    //if(evt == 0) cout<<Eth<<" "<<EvisF<<endl;
    
    if(useOsc){
      weight = survProb3nu(Eth,Dist,th13,th12,Dm31,Dm21);
      //vFarEth.push_back((Float_t)Eth);
      //vFarDist.push_back((Float_t)Dist);
      vFarWeight.push_back((Float_t)weight);
    }

    hFar->Fill(EvisF, weight);
    
    vFarEvis.push_back((Float_t)EvisF);
  }

  TH1D* hFarEscl  = new TH1D("hFarEscl ","",bin_edges.size()-1,&bin_edges[0]);
  TH1D* hFarShape = (TH1D*)hFarEscl->Clone("hFarShape");

  const Int_t nPar = 3;

  TH1D* hParPrime[nPar];
    
  Double_t econstCV_Far2 = 0.0; //0.0091
  Double_t elinCV_Far2   = 1.0; //-0.0041
  Double_t equadCV_Far2  = 0.0;

  Double_t econstUnc_Far2 = 0.00314;
  Double_t elinUnc_Far2   = 0.00255;
  Double_t equadUnc_Far2  = 0.00069;

  if(!isProspect){
    econstUnc_Far2 = 0.016;
    elinUnc_Far2   = 0.0085;
    equadUnc_Far2  = 0.00069;
  }

  TMatrixD* totCov = new TMatrixD(nPar,nPar);
  
  (*totCov )(0,0)= TMath::Power(econstUnc_Far2,2);
  (*totCov )(1,1)= TMath::Power(elinUnc_Far2,2);
  (*totCov )(2,2)= TMath::Power(equadUnc_Far2,2);

  Double_t corr_AB_Far2 = -0.555785;
  Double_t corr_BC_Far2 = -0.265678;
  Double_t corr_AC_Far2 = -0.000121772;
  
  if(!isProspect){
    corr_AB_Far2 = -0.8290;
    corr_BC_Far2 = -0.0807;
    corr_AC_Far2 =  0.0014;
  }

  (*totCov )(0,1)= corr_AB_Far2* econstUnc_Far2*elinUnc_Far2;
  (*totCov )(1,0)= corr_AB_Far2* econstUnc_Far2*elinUnc_Far2;
  (*totCov )(1,2)= corr_BC_Far2* elinUnc_Far2*equadUnc_Far2;
  (*totCov )(2,1)= corr_BC_Far2* elinUnc_Far2*equadUnc_Far2;
  (*totCov )(0,2)= corr_AC_Far2* econstUnc_Far2*equadUnc_Far2;
  (*totCov )(2,0)= corr_AC_Far2* econstUnc_Far2*equadUnc_Far2;

  TRandom* r1[nPar];
  for(Int_t i=0;i<nPar;i++){ r1[i] = new TRandom3(555+i); }

  hParPrime[0] = new TH1D("hEconst","hEconst",100,econstCV_Far2-5.0* econstUnc_Far2,econstCV_Far2+5.0* econstUnc_Far2);
  hParPrime[1] = new TH1D("hElin","hElin",100,elinCV_Far2-5.0* elinUnc_Far2,elinCV_Far2+5.0* elinUnc_Far2);
  hParPrime[2] = new TH1D("hEquad","hEquad",100,equadCV_Far2-5.0* equadUnc_Far2,equadCV_Far2+5.0* equadUnc_Far2);

  TH2D* hParPrime01 = new TH2D("hParPrime01","",100,econstCV_Far2-5.0* econstUnc_Far2,econstCV_Far2+5.0* econstUnc_Far2,100,elinCV_Far2-5.0* elinUnc_Far2,elinCV_Far2+5.0* elinUnc_Far2);
  TH2D* hParPrime02 = new TH2D("hParPrime02","",100,econstCV_Far2-5.0* econstUnc_Far2,econstCV_Far2+5.0* econstUnc_Far2,100,equadCV_Far2-5.0* equadUnc_Far2,equadCV_Far2+5.0* equadUnc_Far2);
  TH2D* hParPrime12 = new TH2D("hParPrime12","",100,elinCV_Far2-5.0* elinUnc_Far2,elinCV_Far2+5.0* elinUnc_Far2,100,equadCV_Far2-5.0* equadUnc_Far2,equadCV_Far2+5.0* equadUnc_Far2);
  
  Int_t count = 0;
  
  //TVectorD simParPrime(nPar);
  
  TDecompChol cholesky1(*totCov);

  Bool_t success = cholesky1.Decompose();
  
  if(success==0){
    cout<<"Choelsky Decomposition failed "<<endl;
    exit(0);
  }
  else if(success==1)
    cout<<"Choelsky Decomposition Success "<<endl;

  TMatrixD lower1 = cholesky1.GetU(); // get upper triangular matrix
  //totCov->Print();
  lower1.Print();

  lower1.T(); // transpose to get lower triangular matrix
  
  Int_t nBinsTot = hFar->GetNbinsX();
  TMatrixD* simMat      = new TMatrixD(nBinsTot,nBinsTot);
  TMatrixD* simMatShape = new TMatrixD(nBinsTot,nBinsTot);

  Int_t nSims = SIMS;

  const Double_t FarNorm  = hFar->Integral();

  TVectorD* random1 = new TVectorD(nPar);
  
  for(Int_t isim = 0; isim<nSims; isim++){

    for(Int_t i=0;i<nPar;i++){(*random1)[i] = r1[i]->Gaus(0,1);}
    //cout<<(*random1)[0]<<" "<<(*random1)[1]<<" "<<(*random1)[2]<<endl;
    TVectorD fluc1 = lower1 * (*random1);
    
    const Double_t thisEscaleA_Far2 = econstCV_Far2 + fluc1[0];
    const Double_t thisEscaleB_Far2 = elinCV_Far2 + fluc1[1];
    const Double_t thisEscaleC_Far2 = equadCV_Far2 + fluc1[2];

    hParPrime[0]->Fill(thisEscaleA_Far2);
    hParPrime[1]->Fill(thisEscaleB_Far2);
    hParPrime[2]->Fill(thisEscaleC_Far2);

    hParPrime01->Fill(thisEscaleA_Far2,thisEscaleB_Far2);
    hParPrime02->Fill(thisEscaleA_Far2,thisEscaleC_Far2);
    hParPrime12->Fill(thisEscaleB_Far2,thisEscaleC_Far2);
    
    //simParPrime(0) = thisEscaleA_Far2;
    //simParPrime(1) = thisEscaleB_Far2;
    //simParPrime(2) = thisEscaleC_Far2;
    
    hFarEscl->Reset();
    hFarShape->Reset();

    if(isim>=count){cout <<"isim: "<<isim<<endl; count += nSims/10;}

    for(Int_t evt = 0; evt<FarEvents; evt++){
      Double_t weight = 1;
      if(useOsc)
	weight = vFarWeight[evt];

      //need to convert back to MeV Evis, since it is the same E scale of DC (approx. better modify?)
      
      Double_t EvisTemp = vFarEvis[evt] / PE2MeV;
      //Double_t EvisTemp = vFarEvis[evt];
      Double_t EvisPrime = thisEscaleA_Far2 + thisEscaleB_Far2 * EvisTemp + thisEscaleC_Far2 * EvisTemp*EvisTemp;
      EvisPrime *= PE2MeV;
      
      //if(evt == 0) cout<<EvisTemp<<" "<<vFarEvis[evt]<<" "<<EvisPrime<<endl;
      
      hFarEscl->Fill(EvisPrime, weight );
    }


    hFarShape = (TH1D*)hFarEscl->Clone("hFarShape");

    Double_t FarNormEscl  = hFarEscl->Integral();

    hFarShape->Scale(FarNorm / FarNormEscl);
    
    //fill matrix
    for(Int_t ibin = 0; ibin < nBinsTot; ++ibin)
      for(Int_t jbin = 0; jbin < nBinsTot; ++jbin){
	(*simMat)(ibin,jbin) += (hFarEscl->GetBinContent(ibin+1)-hFar->GetBinContent(ibin+1))*(hFarEscl->GetBinContent(jbin+1)-hFar->GetBinContent(jbin+1));
	(*simMatShape)(ibin,jbin) += (hFarShape->GetBinContent(ibin+1)-hFar->GetBinContent(ibin+1))*(hFarShape->GetBinContent(jbin+1)-hFar->GetBinContent(jbin+1));
      }
    
    mVariance.SetPoint(isim, (Double_t)isim, 100.0*(TMath::Sqrt( simMat->Sum() / ((Double_t)isim+1.)) /(hFar->Integral()  )));

    mErrorBin12FD.SetPoint(isim, (Double_t)isim, 100.0*(TMath::Sqrt( (*simMat)(11,11) /((Double_t)isim+1.) ) / (hFar->GetBinContent(12) )));
    
  }

  (*simMat) *= 1./nSims;
  (*simMatShape) *= 1./nSims;
  
  cout << "Matrix Norm: "<< 100.0 * (TMath::Sqrt( simMat->Sum())/(hFar->Integral() ))
       <<" percent"<<endl;
  cout << "Shape Matrix Norm: "<< 100.0 * (TMath::Sqrt( simMatShape->Sum())/(hFar->Integral() ))
       <<" percent"<<endl;

  TMatrixD* fracMat = new TMatrixD( *simMat );
  TMatrixD* corrMat = new TMatrixD( *simMat );
  TMatrixD* fracMatShape = new TMatrixD( *simMatShape );
  TMatrixD* corrMatShape = new TMatrixD( *simMatShape );

  for(Int_t ibin = 0; ibin < nBinsTot; ++ibin){
    if( hFar->GetBinContent(ibin+1) == 0. ){
      hFar->SetBinContent(ibin+1,hFar->GetBinContent(ibin));
    }
  }
  
  for(Int_t ibin = 0; ibin < nBinsTot; ++ibin)
    for(Int_t jbin = 0; jbin < nBinsTot; ++jbin)
      {
	if( hFar->GetBinContent(ibin+1) != 0. && hFar->GetBinContent(jbin+1) != 0. ){
	  if( (*simMat)(ibin,ibin)!=0.0 &&  (*simMat)(jbin,jbin)!=0.0 ){
	    (*corrMat)(ibin,jbin) /= TMath::Sqrt( (*simMat)(ibin,ibin)) * TMath::Sqrt((*simMat)(jbin,jbin));
	    (*corrMatShape)(ibin,jbin) /= TMath::Sqrt( (*simMatShape)(ibin,ibin)) * TMath::Sqrt((*simMatShape)(jbin,jbin));
	  }
	  (*fracMat)(ibin,jbin) /= (hFar->GetBinContent(ibin+1) * hFar->GetBinContent(jbin+1) );
	  (*fracMatShape)(ibin,jbin) /= (hFar->GetBinContent(ibin+1) * hFar->GetBinContent(jbin+1) );
	}
      }
   
  TFile* outfile;
  if(!isProspect)
    if(isDYB)
      outfile = new TFile(Form("../data/escaleMatrixDClikeEscale_%.0eevt_ToyMC_%dsims_th12_%4.3f_dm2_%3.2e_DYB.root",(float)FarEvents,nSims,th12,Dm21),"RECREATE");
    else
      outfile = new TFile(Form("../data/escaleMatrixDClikeEscale_%.0eevt_ToyMC_%dsims_th12_%4.3f_dm2_%3.2e_HuberHaag.root",(float)FarEvents,nSims,th12,Dm21),"RECREATE");
  else
    if(isDYB)
      outfile = new TFile(Form("../data/escaleMatrixDClikeEscaleImproved_%.0eevt_ToyMC_%dsims_th12_%4.3f_dm2_%3.2e_DYB.root",(float)FarEvents,nSims,th12,Dm21),"RECREATE");
    else
      outfile = new TFile(Form("../data/escaleMatrixDClikeEscaleImproved_%.0eevt_ToyMC_%dsims_th12_%4.3f_dm2_%3.2e_HuberHaag.root",(float)FarEvents,nSims,th12,Dm21),"RECREATE");
  
  simMat->Write("escale_cov");
  fracMat->Write("escale_frac_cov");
  corrMat->Write("escale_corr");
  simMatShape->Write("escale_cov_shape");
  fracMatShape->Write("escale_frac_cov_shape");
  corrMatShape->Write("escale_corr_shape");
  mVariance.Write("variance_nsims");
  mErrorBin12FD.Write("error_bin12");
  hFar->Write("hFar");
  hParPrime[0]->Write();
  hParPrime[1]->Write();
  hParPrime[2]->Write();
  hParPrime01->Write();
  hParPrime02->Write();
  hParPrime12->Write();
  outfile->Close();

}
