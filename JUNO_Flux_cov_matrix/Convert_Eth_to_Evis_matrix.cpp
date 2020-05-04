//macro to convert the covariance matrix from neutrino energy to visible energy
#include<stdio.h>
#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<math.h>
#include<vector>
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

using namespace std;

enum detectors{
  kFar,
  nDet
};

Double_t DetEvents[nDet];

TString stDet[nDet]={"JUNO"};

void Convert_Eth_to_Evis_matrix(Int_t SIMS=1000, Bool_t isDYB = false){


  //=================== Read Evis flux and put it in a vector : start ====================
  
  
  TFile* fDet[nDet];
  if(!isDYB)
    fDet[kFar]  = new TFile("../data/events_SPMT_1e+06evt_2020_HuberHaag.root");
  else
    fDet[kFar]  = new TFile("../data/events_SPMT_1e+06evt_2020_DYB.root");
  
  TTree* tDet[nDet];
  tDet[kFar] = (TTree*)fDet[kFar]->Get("FinalFitIBDTree");

  std::vector<Float_t> vDet[nDet];
  std::vector<Float_t> vDetEth[nDet];
  
  TGraph mVariance;
  
  Int_t FarEvents = tDet[0]->GetEntries();
  //Int_t FarEvents = 1000000;
  
  Double_t EvisDet[nDet],EthDet[nDet];

  for(int i=0;i<nDet;i++){
    tDet[i]->SetBranchAddress("myPromptEvisID"     ,&EvisDet[i]);
    tDet[i]->SetBranchAddress("myNeutrinoEnergy_Th",&EthDet[i]);
  }
  
  std::vector<Double_t> bin_edges;
  Double_t E = 40.0;
  
  for(int ibin=0; ibin<100; ibin++) // 40 < E < 440 PE
    {
      bin_edges.push_back(E);
      E += 4; // 4 PE bins
    }
  bin_edges.push_back(E);

  TH1D* hDet[nDet];
  for(int i=0;i<nDet;i++)
    hDet[i] = new TH1D(Form("h%s",stDet[i].Data()),stDet[i].Data(),bin_edges.size()-1,&bin_edges[0]);


  for(Int_t evt = 0; evt<FarEvents; evt++){
    for(int i=0;i<nDet;i++){
    
      tDet[i]->GetEntry(evt);
      
      hDet[i]->Fill(EvisDet[i]);

      vDet[i].push_back((Float_t)EvisDet[i]);
      vDetEth[i].push_back((Float_t)EthDet[i]);

      //cout<<EvisDet[i]<<" ";
    }
    //cout<<endl;
  }

  //=================== Read Evis flux and put it in a vector : end ====================


  //=================== Read Predicted flux and cov matrix : start ====================


  TFile* rxtrFile;

  if(!isDYB)
    rxtrFile = new TFile("../data/FluxPredMat_1000sims_HuberHaag_tot.root","READ");
  else
    rxtrFile = new TFile("../data/FluxPredMat_1000sims_DayaBay_tot.root","READ");
  
  TMatrixD* rxtCovMat = (TMatrixD*)rxtrFile->Get("flux_cov");

  TH1D* hDetEth[nDet];
  hDetEth[kFar] = (TH1D*)rxtrFile->Get("hJUNO");

  Int_t nBinsEtTot = hDetEth[0]->GetNbinsX() * nDet;
  Int_t nBinsEt    = hDetEth[0]->GetNbinsX();
  TVectorD *predSigVecsEt[nDet];
  TVectorD *predSigVecEt = new TVectorD(nBinsEtTot);

  for(int i=0;i<nDet;i++){    
    hDet[i]->Sumw2();
    hDet[i]->Scale( (hDetEth[i]->Integral() / FarEvents) );
  }
  
  for(Int_t iBin = 0; iBin<nBinsEt;iBin++){
    hDetEth[kFar]->SetBinError(iBin+1, TMath::Sqrt( (*rxtCovMat)(iBin,iBin) ) );
  }


  Int_t nBinsTot = hDet[0]->GetNbinsX() * nDet;

  Int_t nBins = hDet[0]->GetNbinsX();

  TVectorD *predSigVecs[nDet];

  TVectorD *predSigVec = new TVectorD(nBinsTot);

  Int_t offset = 0;
  for(Int_t det = 0; det < nDet; det++) {
    predSigVecs[det] = new TVectorD(nBins);
    for(Int_t iBin = 0; iBin<nBins;iBin++){
	(*predSigVecs[det])[iBin] = hDet[det]->GetBinContent(iBin+1);
    }
    predSigVec->SetSub(offset, *(predSigVecs[det]) );
    offset += predSigVecs[det]->GetNrows();
  }
  
  /*for(int b = 0;b<nBins;b++){
    cout<<hDet[0]->GetBinCenter(b+1);
    for(int i = 0;i<nDet;i++)
      cout<<" "<<(*predSigVec)[b+i*nBins]<<" ";
    cout<<endl;
    }*/
  
  
  //TCanvas* c1 = new TCanvas("c1","");
  //c1->Divide(2);
  //c1->cd(1);
  //rxtCovMat->Draw("colz");
  //hDetEth[0]->SetLineColor(kRed);
  //hDetEth[0]->Draw();
  /*hDetEth[1]->Draw("same");
  hDetEth[2]->Draw("same");
  hDetEth[3]->SetLineColor(kBlue); 
  hDetEth[3]->Draw("same");
  */
  //c1->cd(2);
  //hDet[0]->SetLineColor(kRed);
  //hDet[0]->Draw();
  
  offset = 0;
  for(Int_t det = 0; det < nDet; det++) {

    predSigVecsEt[det] = new TVectorD(nBinsEt);
    
    for(Int_t iBin = 0; iBin<nBinsEt;iBin++){
	(*predSigVecsEt[det])[iBin] = hDetEth[det]->GetBinContent(iBin+1);
    }
  
    predSigVecEt->SetSub(offset, *(predSigVecsEt[det]) );

    offset += predSigVecsEt[det]->GetNrows();
  }
  /*for(int b = 0;b<nBinsEt;b++){
    cout<<hDetEth[0]->GetBinCenter(b+1);
    for(int i = 0;i<nDet;i++)
      cout<<" "<<(*predSigVecsEt[i])[b]<<" ";
    cout<<endl;
    }*/
  
  //=================== Read Predicted flux and cov matrix : end ====================


  //=================== Randomization process and buile cov matrix in visible energy: start ====================
  
  TH1D* hDetRndEvis[nDet];
  TH1D* hDetRndEvisShape[nDet];
  TH1D* hDetRndEth[nDet];
  
  for(int idet=0;idet<nDet;idet++){
    hDetRndEvis[idet]      = new TH1D(Form("h%sRndEvis",stDet[idet].Data()),stDet[idet].Data(),bin_edges.size()-1,&bin_edges[0]);
    hDetRndEvisShape[idet] = new TH1D(Form("h%sRndEvisShape",stDet[idet].Data()),stDet[idet].Data(),bin_edges.size()-1,&bin_edges[0]);
    hDetRndEth[idet]       = (TH1D*)( hDetEth[idet]->Clone(Form("h%sRndEth",stDet[idet].Data()) ) );
  }

  TRandom* r1[nBinsEtTot];
  for(Int_t i=0;i<nBinsEtTot;i++){r1[i] = new TRandom3(0.0);}
  Int_t count = 0;

  TMatrixD* simMat      = new TMatrixD(nBinsTot,nBinsTot);
  TMatrixD* simMatEt    = new TMatrixD(nBinsEtTot,nBinsEtTot);
  TMatrixD* simMatShape = new TMatrixD(nBinsTot,nBinsTot);

  TVectorD *predVars[nDet], *predVarsShape[nDet];
  TVectorD *predVarsEt[nDet], *predVarsEtShape[nDet];

  TVectorD  *predVar        = new TVectorD(nBinsTot);
  TVectorD  *predVarShape   = new TVectorD(nBinsTot);

  TVectorD  *predVarEt      = new TVectorD(nBinsEtTot);
  TVectorD  *predVarEtShape = new TVectorD(nBinsEtTot);
  
  Int_t nSims = SIMS;

  TDecompChol cholesky1(*rxtCovMat);

  Bool_t success = cholesky1.Decompose();
  
  if(success==0){
    cout<<"Choelsky Decomposition failed "<<endl;
    exit(0);
  }
  else if(success==1)
    cout<<"Choelsky Decomposition Success "<<endl;

  TMatrixD lower1 = cholesky1.GetU(); // get upper triangular matrix

  //lower1.Print();

  lower1.T(); // transpose to get lower triangular matrix

  for(Int_t isim = 0; isim<nSims; isim++){

    Double_t NormTemp[nDet];
    
    for(int idet=0;idet<nDet;idet++){
      hDetRndEvis[idet]->Reset();
      hDetRndEvisShape[idet]->Reset();
      hDetRndEth[idet]->Reset();
    }

    TVectorD* random1=new TVectorD(nBinsEtTot);

    for(Int_t i=0;i<nBinsEtTot;i++){(*random1)[i] = r1[i]->Gaus(0,1);}
    
    TVectorD fluc1 = lower1 * (*random1);

    for(Int_t ibin=0;ibin<nBinsEtTot;ibin++){
      (*predVarEt)(ibin) = (*predSigVecEt)(ibin)  + fluc1[ibin];
    }

    for(Int_t idet=0;idet<nDet;idet++){
      for(Int_t b=0;b<nBinsEt;b++){
	hDetRndEth[idet]->SetBinContent(b+1, (*predVarEt)(nBinsEt*idet+b));
      }
      //cout<<hDetRndEth[idet]->GetBinContent(1)<<" ";
      //cout<<(*predVarEt)(nBinsEt*idet)<<" ";

      NormTemp[idet] = hDetRndEth[idet]->Integral();
    }
    //cout<<endl;      
    
    //c1->cd(1);
    //hDetRndEth[0]->SetLineWidth(1);
    //hDetRndEth[0]->SetLineColor(1);
    //hDetRndEth[0]->DrawClone("same");
    
    if(isim>=count){cout <<"isim: "<<isim<<endl; count += nSims/10;}

    for(Int_t evt = 0; evt<FarEvents; evt++){
      for(int idet=0;idet<nDet;idet++){

	Double_t weight = 1;

	Int_t BinDenominator = hDetEth[idet]->GetXaxis()->FindBin(vDetEth[idet][evt]);
	//Int_t BinNumerator   = hDetRndEth[idet]->GetXaxis()->FindBin(vDetEth[idet][evt]);

	Double_t Denominator = hDetEth[idet]->GetBinContent( BinDenominator );
	Double_t Numerator   = hDetRndEth[idet]->GetBinContent( BinDenominator );
	
	weight *=  Numerator / Denominator;

	//weight *=   NormTemp[idet]  / FarEvents;
	weight *=   hDetEth[idet]->Integral()  / FarEvents;

	//if(idet==0) cout<<vDetEth[idet][evt]<<" "<<vDet[idet][evt]<<" "<< Numerator / Denominator <<" "<<NormTemp[idet]  / FarEvents<<" "<<weight<<" "<<BinDenominator<<endl;
	
	hDetRndEvis[idet]->Fill(vDet[idet][evt], weight);
      }
    }
    //cout<< hDetEth[1]->Integral()<<" "<<NormTemp[1]<<" "<<hDetRndEvis[1]->Integral()<<" "<<NormTemp[1]*hDet[1]->Integral() / FarEvents << " "<< hDetEth[1]->Integral()*hDet[1]->Integral()/FarEvents<<endl;
    //cout<< hDetEth[1]->Integral()<<" "<<NormTemp[1]<<" "<<hDetRndEvis[1]->Integral()<<" "<<hDet[1]->Integral() << endl;
    // Make a histogram with this three last variables to compare distributions?? no need now, because they cannot be compared direclty (at visible energy histogram, events with Evis<1MeV are not taken into account, affecting the final rate)
    
    //c1->cd(2);
    //hDetRndEvis[0]->SetLineWidth(1);
    //hDetRndEvis[0]->SetLineColor(1);
    //hDetRndEvis[0]->DrawClone("histsame");
  
    //hDet[0]->SetLineWidth(3);

    //hDet[0]->GetXaxis()->SetRangeUser(1,7.5);

    //hDetRndEvisShape = (TH1D*)hDetRndEvis->Clone("");
    //Double_t NearNormEscl = hNearEscl->Integral();
    //hNearShape->Scale(NearNorm / NearNormEscl);    
    
    offset = 0;
    for(Int_t idet=0;idet<nDet;idet++){
      predVars[idet] = new TVectorD(nBins);
      //predVarsShape[det] = new TVectorD(nBins);
      for(Int_t iBin=0;iBin<nBins;iBin++){
	(*predVars[idet])[iBin]      = hDetRndEvis[idet]->GetBinContent(iBin+1);
	//(*predVarsShape[det])[iBin] = hFarShape->GetBinContent(iBin+1);
      }
      predVar->SetSub(offset, *(predVars[idet]) );
      //predVarShape->SetSub(offset, *(predVarsShape[idet]) );
      offset += predVars[idet]->GetNrows();
    }

    offset = 0;
    for(Int_t idet=0;idet<nDet;idet++){
      predVarsEt[idet] = new TVectorD(nBinsEt);
      for(Int_t iBin=0;iBin<nBinsEt;iBin++){
	(*predVarsEt[idet])[iBin]      = hDetRndEth[idet]->GetBinContent(iBin+1);
      }
      predVarEt->SetSub(offset, *(predVarsEt[idet]) );
      offset += predVarsEt[idet]->GetNrows();
    }
    
    //fill matrix
    for(Int_t ibin = 0; ibin < nBinsTot; ++ibin)
      for(Int_t jbin = 0; jbin < nBinsTot; ++jbin){
	(*simMat)(ibin,jbin)      += ((*predVar)[ibin]      - (*predSigVec)[ibin]) * ((*predVar)[jbin]      - (*predSigVec)[jbin] );
	//(*simMatShape)(ibin,jbin) += ((*predVarShape)[ibin] - (*predSigVec)[ibin]) * ((*predVarShape)[jbin] - (*predSigVec)[jbin] );
      }
    mVariance.SetPoint(isim, (Double_t)isim, 100.0*(TMath::Sqrt( simMat->Sum() / ((Double_t)isim+1.)) /(predSigVec->Sum() )));
    //mErrorBin12FD.SetPoint(isim, (Double_t)isim, 100.0*(TMath::Sqrt( (*simMat)(11,11) /((Double_t)isim+1.) ) / (hFar->GetBinContent(12) )));
    
    for(Int_t ibin = 0; ibin < nBinsEtTot; ++ibin)
      for(Int_t jbin = 0; jbin < nBinsEtTot; ++jbin){
	(*simMatEt)(ibin,jbin)      += ((*predVarEt)[ibin]      - (*predSigVecEt)[ibin]) * ((*predVarEt)[jbin]      - (*predSigVecEt)[jbin] );
      }
    
    delete random1;
    random1 = NULL;
  }

  //c1->cd(2);
  //hDet[0]->DrawClone("same");
  //gPad->SetLogy();
  //c1->cd(1);
  //hDetEth[0]->DrawClone("same");
  //gPad->SetLogy();
  //c1->Print("test.png");

  
  (*simMat) *= 1./nSims;
  (*simMatEt) *= 1./nSims;
  //(*simMatShape) *= 1./nSims;

  Double_t SumDet[nDet] = {0.0};
  Double_t SumDetEt[nDet] = {0.0};

  for(int idet=0;idet<nDet;idet++){
    for(Int_t ibin = 0; ibin < nBins; ++ibin)
      for(Int_t jbin = 0; jbin < nBins; ++jbin){
	SumDet[idet] += (*simMat)(ibin+idet*nBins , jbin+idet*nBins);
      }
    
    cout << "Matrix Norm("<<stDet[idet].Data() <<"): "<< 100.0 * (TMath::Sqrt( SumDet[idet])/(hDet[idet]->Integral() ))
    <<" percent"<<endl;
    
    for(Int_t ibin = 0; ibin < nBinsEt; ++ibin)
      for(Int_t jbin = 0; jbin < nBinsEt; ++jbin){
	SumDetEt[idet] += (*simMatEt)(ibin+idet*nBinsEt , jbin+idet*nBinsEt);
      }
    
    //cout << "Matrix Norm Eth ("<<stDet[idet].Data() <<"): "<< 100.0 * (TMath::Sqrt( SumDetEt[idet])/(hDetEth[idet]->Integral() ))
    //<<" percent"<<endl;

  }
  
  //cout << "Matrix Norm: "<< 100.0 * (TMath::Sqrt( simMat->Sum())/(hDet[0]->Integral() + hDet[1]->Integral()+ hDet[2]->Integral()+ hDet[3]->Integral() ))
  //<<" percent"<<endl;
  //cout << "Shape Matrix Norm: "<< 100.0 * (TMath::Sqrt( simMatShape->Sum())/(hFar->Integral() + hNear->Integral() ))
  //<<" percent"<<endl;

  TMatrixD* fracMat = new TMatrixD( *simMat );
  TMatrixD* corrMat = new TMatrixD( *simMat );
  //TMatrixD* fracMatShape = new TMatrixD( *simMatShape );
  //TMatrixD* corrMatShape = new TMatrixD( *simMatShape );

  for(Int_t ibin = 0; ibin < nBinsTot; ++ibin){
    if( (*predSigVec)[ibin] == 0. ){
      (*predSigVec)[ibin] = (*predSigVec)[ibin-1];
    }
  }
   
  for(Int_t ibin = 0; ibin < nBinsTot; ++ibin)
    for(Int_t jbin = 0; jbin < nBinsTot; ++jbin)
      {
	if( (*predSigVec)[ibin] != 0. && (*predSigVec)[jbin] != 0. ){
	  
	  if( (*simMat)(ibin,ibin)!=0.0 &&  (*simMat)(jbin,jbin)!=0.0 ){
	    (*corrMat)(ibin,jbin) /= TMath::Sqrt( (*simMat)(ibin,ibin)) * TMath::Sqrt((*simMat)(jbin,jbin));
	    //(*corrMatShape)(ibin,jbin) /= TMath::Sqrt( (*simMatShape)(ibin,ibin)) * TMath::Sqrt((*simMatShape)(jbin,jbin));
	  }
	  
	  (*fracMat)(ibin,jbin) /= (*predSigVec)[ibin]*(*predSigVec)[jbin];
	  //(*fracMatShape)(ibin,jbin) /= (*predSigVec)[ibin]*(*predSigVec)[jbin];
	  
	}
      }

  TMatrixD* fracMatEt = new TMatrixD( *simMatEt );
  TMatrixD* corrMatEt = new TMatrixD( *simMatEt );

   for(Int_t ibin = 0; ibin < nBinsEtTot; ++ibin){
     if( (*predSigVecEt)[ibin] == 0. ){
       (*predSigVecEt)[ibin] = (*predSigVecEt)[ibin-1];
     }
   }
   
   for(Int_t ibin = 0; ibin < nBinsEtTot; ++ibin)
     for(Int_t jbin = 0; jbin < nBinsEtTot; ++jbin)
       {
	 if( (*predSigVecEt)[ibin] != 0. && (*predSigVecEt)[jbin] != 0. ){

	   if( (*simMatEt)(ibin,ibin)!=0.0 &&  (*simMatEt)(jbin,jbin)!=0.0 ){
	     (*corrMatEt)(ibin,jbin) /= TMath::Sqrt( (*simMatEt)(ibin,ibin)) * TMath::Sqrt((*simMatEt)(jbin,jbin));
	   }
	   
	   (*fracMatEt)(ibin,jbin) /= (*predSigVecEt)[ibin]*(*predSigVecEt)[jbin];

	 }
       }

   TFile* outfile;

   if(!isDYB)
     outfile = new TFile(Form("../data/FluxMatrixEvis_%.0eevt_%dsims_HuberHaag.root",(float)FarEvents,nSims),"RECREATE");
   else
     outfile = new TFile(Form("../data/FluxMatrixEvis_%.0eevt_%dsims_DYB.root",(float)FarEvents,nSims),"RECREATE");


   simMatEt->Write("orig_cov");
   fracMatEt->Write("orig_frac");
   corrMatEt->Write("orig_cor");
   simMat->Write("cov_approx");
   fracMat->Write("frac_approx");
   corrMat->Write("corr_approx");
   //simMatShape->Write("escale_cov_shape");
   //fracMatShape->Write("escale_frac_cov_shape");
   //corrMatShape->Write("escale_corr_shape");
   mVariance.Write("variance_nsims");
   //mErrorBin12FD.Write("error_bin12");
   for(int i=0;i<nDet;i++)
     hDet[i]->Write();
   //t.Write();
   outfile->Close();
   
  //=================== Randomization process and buile cov matrix in visible energy: end ====================

}
