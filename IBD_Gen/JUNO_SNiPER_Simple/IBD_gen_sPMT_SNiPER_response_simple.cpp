//--------------------------------------------------------------
//--------------------------------------------------------------
//
//  Simple Reactor Anti-neutrino Event Generator
//
//  Author: T.J.C.Bezerra
//
//--------------------------------------------------------------
//--------------------------------------------------------------
//
// Upgrades:
//
// March 5th 2010: Finish of the first version
//
// August 8th 2013: Code all re-written using the kinematics of the Inverse Beta Decay as described in arvix:9903554
//
// October 17th 2018: Added cilyndrical detector option, and two reactors at different baselines
//
// December 18th 2018: Added spheric detector option.
//
// April 25th 2019: Updated E scale and resolution of JUNO SPMTs case
//                  Fixed bug on position generator of spheric detector
//
// November 13th 2019 : version using JUNO's SNiPER response && a simpler IBD kinematics
//

#include<stdio.h>
#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<math.h>
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

using namespace std;

const Double_t PI     = TMath::Pi();
const Double_t LowE   = 1.806;// 1.875;
const Double_t HighE  = 12.125;//8.125;

const Double_t m_e    = 0.51099893;//MeV
const Double_t M_p    = 938.27203;//MeV
const Double_t M_n    = 939.56536;//MeV
const Double_t dM     = M_n - M_p;//MeV
const Double_t M      = (M_p+M_n)/2;//MeV
const Double_t f      = 1.0;
const Double_t g      = 1.267;
const Double_t sigma0 = 0.952e-42/(f*f+3*g*g);
const Double_t f2     = 3.706;//mu_p - mu_n
const Double_t y2     = (dM*dM - m_e*m_e)/2.;

TH1D*  hDetHist;

enum TypeDet
  {
    kFar,
    kDet_max
  };

enum TypeRxt
  {
    kB1,  //yangjian
    kB2,  //taishan
    kRxt_max
  };

//Detector names
TString tDet[kDet_max]={"Far"};


// position in mm on a frame where the origin is a mid-point btw reactors centers

TVector3 vReactorB1Pos( -39.05e6,  0.0    , 0.0);
TVector3 vReactorB2Pos(  39.05e6,  0.0    , 0.0);
TVector3 vReactorFDPos(  0.0    ,  35.09e6, 0.0);

Double_t fluxFrac[kDet_max][kRxt_max] = {0.49 , 0.51};  //[FD][B1,B2]

Double_t isotope;

Int_t seed=0;

TF1* THespU235;
TF1* THespPu239;
TF1* THespPu241;
TF1* THespU238;

//Mean fission rates fraction (Mean of DC reactors B1 and B2)
Double_t pU235  = 0.52;
Double_t pPu239 = 0.33;
Double_t pPu241 = 0.06;
Double_t pU238  = 0.09;

Float_t baseline;

TRandom *r0, *r1, *r2, *r3, *r4, *r5, *r6, *r7, *r8, *r9, *r10;
  
Int_t GetReactor(Int_t det)
{
  Int_t wReactor;

  Double_t rRxt;

  rRxt = r0->Rndm();

  if(rRxt > fluxFrac[det][kB1])
    wReactor = kB2;
  else
    wReactor = kB1;

  return wReactor;
}

TVector3 GetPosition(Double_t Detector_Lx, Double_t Detector_Ly, Double_t Detector_Lz) // IBD interaction position at detector frame | box
{
  Double_t xv,yv,zv,rv,thetav;
  TVector3 vpos;

  xv = Detector_Lx*(2.*r1->Rndm()-1)/2;
  yv = Detector_Ly*(2.*r2->Rndm()-1)/2;
  zv = Detector_Lz*(2.*r3->Rndm()-1)/2;
  //cout<<xv<<" "<<yv<<" "<<zv<<endl;
  
  vpos.SetXYZ(xv,yv,zv);
  
  return vpos;
}
TVector3 GetPosition(Double_t Detector_R, Double_t Detector_H) // IBD interaction position at detector frame | cylinder
{
  Double_t xv,yv,zv,rv,thetav;
  TVector3 vpos;

  rv     = r1->Rndm();
  rv     = Detector_R * TMath::Sqrt(rv);
  thetav = 2.0 * PI  * r2->Rndm();
  zv     = Detector_H*(2.*r3->Rndm()-1)/2;
  xv     = rv * TMath::Cos(thetav);
  yv     = rv * TMath::Sin(thetav);
  //cout<<xv<<" "<<yv<<" "<<zv<<endl;
  
  vpos.SetXYZ(xv,yv,zv);
  
  return vpos;
}

TVector3 GetPosition(Double_t Detector_R) // IBD interaction position at detector frame | sphere
{
  Double_t xv,yv,zv,rv,thetav,phiv;
  TVector3 vpos;

  thetav = 2.0 * PI  * r2->Rndm();
  phiv = TMath::ACos(2.0 * r1->Rndm() - 1.0);
  
  Double_t rcube= r3->Rndm()*TMath::Power(Detector_R,3.0)/3.0;

  rv     = TMath::Power(rcube*3 , 1./3.);

  Double_t sinTheta = TMath::Sin(thetav);
  Double_t cosTheta = TMath::Cos(thetav);
  Double_t sinPhi   = TMath::Sin(phiv);
  Double_t cosPhi   = TMath::Cos(phiv);

  xv = rv * sinPhi * cosTheta;
  yv = rv * sinPhi * sinTheta;
  zv = rv * cosPhi;
  
  //cout<<xv<<" "<<yv<<" "<<zv<<endl;
  
  vpos.SetXYZ(xv,yv,zv);
  
  return vpos;
}

TVector3 GetBaseline(TVector3 vpos, Int_t nRxt, Int_t det)
{
  TVector3 fNeutrinoMomemtum;

  if(nRxt == kB1){
    fNeutrinoMomemtum = vpos + vReactorFDPos - vReactorB1Pos;}
  else{
    fNeutrinoMomemtum = vpos + vReactorFDPos - vReactorB2Pos;}

  baseline = fNeutrinoMomemtum.Mag();

  return fNeutrinoMomemtum;
}

Double_t GetAntiNuE(){
  Double_t Enue;

  Enue = hDetHist->GetRandom();

  return Enue;
}

Double_t GetPositronE(Double_t Enu){
  //Energy of the positron produced by an antineutrino of energy Enu
  // in the lab frame.
  //
  //If neutrino is below threshold, function returns -1.
  //

  const Double_t Mn = M_n;
  const Double_t Mp = M_p;
  const Double_t Me = m_e;
  const Double_t Delta = Mn - Mp;

  const Double_t Mdiff = -Enu + Delta + (Delta*Delta - Me*Me)/(2*Mp);

  //Checks to make sure kinematics are appropriate

  if( Mdiff < 0 )
    {
      const Double_t Epos = (-Mn + TMath::Sqrt(Mn*Mn - 4*Mp*Mdiff))/2;
      if( Epos > Me )
	return Epos;
      else
	return -1.;
    }
  else
    return -1.;

}

void IBD_Gen(int date = 2020, int nevt = 100, int initrand = 0, Int_t opt = 3, Double_t Detector_Lx=17.7e3, Double_t Detector_Ly=0.0, Double_t Detector_Lz=0.0, Int_t det=kFar, Bool_t addElecMass=false) //if Lz = 0, then cylindrical volume (Lx = radius, Ly = height)
{
  
  delete gRandom;
  gRandom = new TRandom3(initrand);

  //TFile fFluxFile("FluxPredMat_100000sims_Huber_tot_config1.root");
  TFile fFluxFile("FluxPredMat_10000sims_DayaBay_tot.root");

  hDetHist = (TH1D*)fFluxFile.Get(Form("h%s",tDet[det].Data()));

  
  
  TVector3 ps;
  TVector3 P_posi, P_neu, P_nu;
  
  Float_t E_nue, E_pos, E_neu, Ene;
  Int_t    nRxt;

  //const Double_t PE2MeV = 52.57;
  Double_t PE2MeV = 52.57;
  const Double_t Eresol = 0.18;
  
  seed = initrand;

  r0 = new TRandom3(seed);
  r1 = new TRandom3(seed);
  r2 = new TRandom3(seed);
  r3 = new TRandom3(seed);
  r4 = new TRandom3(seed);
  r5 = new TRandom3(seed);
  r6 = new TRandom3(seed);
  r7 = new TRandom3(seed);
  r8 = new TRandom3(seed);
  r9 = new TRandom3(seed);
  
  while((opt < 0) || (opt > 3))
  {
    cout<<"Choose one of the following options to print:\n 0: Print both Positron or Neutron event \n 1: Print only the Positron \n 2: Print only the Neutron\n 3: Print only Positron Energy";
    cin>>opt;
  }

  TFile *newFile;
  //newFile = new TFile(Form("events_SPMT_%.0eevt_FFIT.root",(float) nevt),"recreate");
  newFile = new TFile(Form("events_SPMT_%.0eevt_%d_test_DYB.root",(float) nevt,date),"recreate");

  TTree *newTree = new TTree("FinalFitIBDTree","FinalFitIBDTree");

  Double_t myNeutrinoEnergy_Th,myPromptEvisID,myNeutrinoDistance_Th, X, Y, Z, myPositronE;
  Int_t volumePM=1,volumeDL=1,volumeGDML=1,myNeutrinoReactor_Th=1,myPromptBackgroundType_Th=0;
  
  newTree->Branch("myNeutrinoEnergy_Th"  , &myNeutrinoEnergy_Th  , "myNeutrinoEnergy_Th/D"  );
  newTree->Branch("myPromptEvisID"       , &myPromptEvisID       , "myPromptEvisID/D"       );
  newTree->Branch("myNeutrinoDistance_Th", &myNeutrinoDistance_Th, "myNeutrinoDistance_Th/D");
  newTree->Branch("myPositronE"          , &myPositronE          , "myPositronE/D"          );

  newTree->Branch("volumePM", &volumePM, "volumePM/I");
  newTree->Branch("volumeDL", &volumeDL, "volumeDL/I");
  newTree->Branch("volumeGDML", &volumeGDML, "volumeGDML/I");
  newTree->Branch("myNeutrinoReactor_Th", &myNeutrinoReactor_Th, "myNeutrinoReactor_Th/I");
  newTree->Branch("myPromptBackgroundType_Th", &myPromptBackgroundType_Th, "myPromptBackgroundType_Th/I");
  newTree->Branch("X",&X,"X/D");
  newTree->Branch("Y",&Y,"Y/D");
  newTree->Branch("Z",&Z,"Z/D");

  Double_t lowEpos   = 0.5;
  Double_t higEpos   = 10.5;
  Int_t    nEpos     = 20;
  Double_t EposWidth = (higEpos - lowEpos) / nEpos;

  Double_t lowR3   = 0.0;
  Double_t higR3   = 5400;
  Int_t    nR3     = 20;
  Double_t R3Width = (higR3 - lowR3) / nR3;
  
  ifstream fRead("PE2MeV_vs_MeV_vs_R3.txt");
  Double_t PE2MeV_cv[20][20], PE2MeV_er[20][20], dDummy;
  Int_t iDummy;

  TRandom3 rGaus(0);

  for(int e=0;e<20;e++){
    for(int r=0;r<20;r++){
      fRead>>iDummy>>iDummy>>PE2MeV_cv[e][r]>>PE2MeV_er[e][r]>>dDummy>>dDummy>>dDummy; // Energy Counter, Volume Counter, PE2MeV_cv, PE2MeV_sigma, E_resol_%, gaus prob, nEvents on histo
      //cout<<PE2MeV_cv[e][r]<<" "<<PE2MeV_er[e][r]<<endl;
    }
  }
  
  for(int i=0;i<nevt;i++)
  {
    nRxt = GetReactor(det);
    if(Detector_Ly == 0.0 && Detector_Lz == 0.0)
      ps   = GetPosition(Detector_Lx);  // spherical detector
    else if (Detector_Lz == 0.0)
      ps   = GetPosition(Detector_Lx, Detector_Ly); // cylindrical detector
    else
      ps   = GetPosition(Detector_Lx, Detector_Ly, Detector_Lz); // box detector

    P_nu  = GetBaseline(ps, nRxt, det);
    
    E_nue = GetAntiNuE();

    E_pos = GetPositronE(E_nue); Ene = E_pos;

    if(addElecMass) E_pos += m_e;  //needed in case it is not a G4 simulation

    X = ps.X(); Y = ps.Y(); Z = ps.Z();
    Double_t R3 = pow(sqrt(X*X + Y*Y + Z*Z )/1000.0,3);

    //E_pos = 1.6;
    //R3 = 4100;

    Int_t ratioR = R3/R3Width;
    Int_t ratioE = E_pos/EposWidth -1;
    //cout<<E_pos<<" "<<E_pos/EposWidth -1<<" "<<ratioE<<" | "<<R3<<" "<<R3/R3Width<<" "<<ratioR<<" "<<endl;

    myPositronE = E_pos;

    PE2MeV = PE2MeV_cv[ratioE][ratioR] / ( (EposWidth/2.0) + ratioE*EposWidth + EposWidth);
    
    Double_t sig = (PE2MeV_er[ratioE][ratioR] / PE2MeV_cv[ratioE][ratioR]);
    Double_t nPe = PE2MeV * E_pos;

    nPe =  rGaus.Gaus(nPe,sig*nPe);      //  PE from SNiPER
    
    Ene = nPe;
    //########### Fill tree ################
    
    myNeutrinoEnergy_Th   = E_nue;
    myPromptEvisID        = Ene;
    myNeutrinoDistance_Th = baseline*1e-3;
    myNeutrinoReactor_Th  = nRxt + 1;    
    if(myNeutrinoEnergy_Th<8.125 && myPositronE>0) newTree->Fill();
    
  }

  newTree->Write();
  newFile->Close();

}

