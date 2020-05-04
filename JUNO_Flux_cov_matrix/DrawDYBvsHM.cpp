//Macro to draw the Prediction histogram and comparison, in visible energy (PE), using DYB and H.-M. as prediction

{

  TFile fHM("../data/FluxMatrixEvis_1e+06evt_1000sims_HuberHaag.root");
  TH1D* hHM = (TH1D*)fHM.Get("hJUNO");
  TH1D* hHMfrac  = new TH1D("hHMfrac","",100,40,440);
  TH1D* hHMratio = new TH1D("hHMratio","",100,40,440);
  TMatrixD* mHM = (TMatrixD*)fHM.Get("cov_approx");
  TMatrixD* mHMcor  = (TMatrixD*)fHM.Get("corr_approx");
  TMatrixD* mHMfrac = (TMatrixD*)fHM.Get("frac_approx");

  TH2D* hHMcov = new TH2D("hHMcov","",100,40,440,100,40,440);
  TH2D* hHMcor = new TH2D("hHMcor","",100,40,440,100,40,440);
  
  for(int i=0;i<100;i++){
    hHM->SetBinError(i+1, TMath::Sqrt( (*mHM)(i,i)  ) );
    hHMfrac->SetBinContent(i+1, 100.0*TMath::Sqrt( (*mHMfrac)(i,i)  ) );
    hHMratio->SetBinError(i+1, TMath::Sqrt( (*mHMfrac)(i,i)  ) );
    hHMratio->SetBinContent(i+1,1.0);
    for(int j=0;j<100;j++){
      hHMcov->SetBinContent(i+1,j+1,(*mHM)(i,j) );
      hHMcor->SetBinContent(i+1,j+1,(*mHMcor)(i,j) );
    }
  }

  TCanvas c1("c1","");
  
  hHM->SetFillColor(kBlue);
  hHM->SetLineColor(kBlue);
  hHM->SetMarkerColor(kBlue);
  hHM->GetXaxis()->SetRangeUser(40,440);
  hHM->GetXaxis()->SetTitle("Visible Energy (PE)");
  hHM->GetYaxis()->SetTitle("Events / 4 PE");
  hHM->GetYaxis()->SetTitleOffset(1.5);
  hHM->Draw("E2");

  TFile fDYB("../data/FluxMatrixEvis_1e+06evt_1000sims_DYB.root");
  TH1D* hDYB = (TH1D*)fDYB.Get("hJUNO");
  TH1D* hDYBfrac = new TH1D("hDYBfrac","",100,40,440);
  TH1D* hDYBratio = new TH1D("hDYBratio","",100,40,440);
  TMatrixD* mDYB = (TMatrixD*)fDYB.Get("cov_approx");
  TMatrixD* mDYBcor = (TMatrixD*)fDYB.Get("corr_approx");
  TMatrixD* mDYBfrac = (TMatrixD*)fDYB.Get("frac_approx");

  TH2D* hDYBcov = new TH2D("hDYBcov","",100,40,440,100,40,440);
  TH2D* hDYBcor = new TH2D("hDYBcor","",100,40,440,100,40,440);
  
  for(int i=0;i<100;i++){
    hDYB->SetBinError(i+1, TMath::Sqrt( (*mDYB)(i,i)  ) );
    hDYBfrac->SetBinContent(i+1, 100.0*TMath::Sqrt( (*mDYBfrac)(i,i)  ) );
    hDYBratio->SetBinError(i+1, TMath::Sqrt( (*mDYBfrac)(i,i)  ) );
    Double_t binRatio;
    if( hHM->GetBinContent(i+1)>0)
      binRatio = hDYB->GetBinContent(i+1) / hHM->GetBinContent(i+1);
    hDYBratio->SetBinContent(i+1,binRatio);
    for(int j=0;j<100;j++){
      hDYBcov->SetBinContent(i+1,j+1,(*mDYB)(i,j) );
      hDYBcor->SetBinContent(i+1,j+1,(*mDYBcor)(i,j) );
    }
  }

  hDYB->SetFillColor(kRed);
  hDYB->SetLineColor(kRed);
  hDYB->SetLineWidth(0);
  hDYB->SetFillStyle(3001);
  hDYB->SetMarkerColor(kRed);
  hHM->SetTitle("H.-M.");
  hDYB->SetTitle("DYB");
  //hDYB->GetXaxis()->SetRangeUser(1,8);
  hDYB->Draw("sameE2");
  gPad->BuildLegend();
  hHM->SetTitle("Prediction Spectra");

  //c1.SaveAs("DYB_vs_HM_Spectra.pdf");
  //c1.SaveAs("DYB_vs_HM_Spectra.png");
  //c1.SaveAs("DYB_vs_HM_Spectra.root");
  
  TCanvas c2("c2","");
  hHMfrac->SetFillColor(kBlue);
  hHMfrac->SetFillStyle(3001);
  hDYBfrac->SetFillColor(kRed);
  hDYBfrac->SetLineWidth(2);
  hHMfrac->SetLineWidth(2);
  hHMfrac->SetTitle("H.-M.");
  hDYBfrac->SetTitle("DYB");
  hDYBfrac->GetXaxis()->SetTitle("Visible Energy (PE)");
  hDYBfrac->GetYaxis()->SetTitle("Bin Error (%)");
  hDYBfrac->Draw();
  hHMfrac->Draw("same");
  gPad->BuildLegend(0.310172,0.691083,0.689828,0.901274);
  hDYBfrac->SetTitle("Prediction Spectra Uncertainty");

  //c2.SaveAs("DYB_vs_HM_BinFracError.pdf");
  //c2.SaveAs("DYB_vs_HM_BinFracError.png");
  //c2.SaveAs("DYB_vs_HM_BinFracError.root");

  TCanvas c3("c3","",1294,500);
  hDYBcov->SetTitle("DYB Covariance Matrix");
  hDYBcov->GetXaxis()->SetTitle("Visible Energy (PE)");
  hDYBcov->GetYaxis()->SetTitle("Visible Energy (PE)");
  hDYBcov->GetYaxis()->SetTitleOffset(1.5);
  hHMcov->SetTitle("H.-M. Covariance Matrix");
  hHMcov->GetXaxis()->SetTitle("Visible Energy (PE)");
  hHMcov->GetYaxis()->SetTitle("Visible Energy (PE)");
  hHMcov->GetYaxis()->SetTitleOffset(1.5);
  
  c3.Divide(2,1);
  c3.cd(1);
  hHMcov->Draw("colz");
  c3.cd(2);
  hDYBcov->Draw("colz");

  //c3.SaveAs("DYB_vs_HM_CovMatNoOsc.pdf");
  //c3.SaveAs("DYB_vs_HM_CovMatNoOsc.png");
  //c3.SaveAs("DYB_vs_HM_CovMatNoOsc.root");
  
  TCanvas c4("c4","",1294,500);
  const int NRGBs =7;
  int NCont=999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.000, 0.167, 0.333, 0.500, 0.667, 0.833, 1.000 };
  Double_t red[NRGBs]   = { 0.600, 1.000, 1.000, 1.000, 0.600, 0.000, 0.000 };
  Double_t green[NRGBs] = { 0.000, 0.000, 0.600, 1.000, 0.600, 0.000, 0.000 };
  Double_t blue[NRGBs]  = { 0.000, 0.000, 0.600, 1.000, 1.000, 1.000, 0.600 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

  hDYBcor->SetTitle("DYB Correlation Matrix");
  hDYBcor->GetXaxis()->SetTitle("Visible Energy (PE)");
  hDYBcor->GetYaxis()->SetTitle("Visible Energy (PE)");
  hDYBcor->GetYaxis()->SetTitleOffset(1.5);
  hHMcor->SetTitle("H.-M. Correlation Matrix");
  hHMcor->GetXaxis()->SetTitle("Visible Energy (PE)");
  hHMcor->GetYaxis()->SetTitle("Visible Energy (PE)");
  hHMcor->GetYaxis()->SetTitleOffset(1.5);
    
  c4.Divide(2,1);
  c4.cd(1);
  hHMcor->Draw("colz");
  c4.cd(2);
  hDYBcor->Draw("colz");

  hHMcor->GetZaxis()->SetRangeUser(-1,1);
  hDYBcor->GetZaxis()->SetRangeUser(-1,1);

  //c4.SaveAs("DYB_vs_HM_CorrMatNoOsc.pdf");
  //c4.SaveAs("DYB_vs_HM_CorrMatNoOsc.png");
  //c4.SaveAs("DYB_vs_HM_CorrMatNoOsc.root");


  TCanvas c5("c5","");
  hHMratio->SetFillColor(kBlue);
  hHMratio->SetLineColor(kBlue);
  hHMratio->SetMarkerColor(kBlue);
  hHMratio->GetXaxis()->SetRangeUser(40,440);
  hDYBratio->SetFillColor(kRed);
  hDYBratio->SetLineColor(kRed);
  hDYBratio->SetLineWidth(0);
  hDYBratio->SetFillStyle(3001);
  hDYBratio->SetMarkerColor(kRed);
  hHMratio->GetXaxis()->SetTitle("Visible Energy (PE)");
  hHMratio->GetYaxis()->SetTitle("Ratio wrt H.-M.");
  hHMratio->GetYaxis()->SetRangeUser(0.7,1.2);
  hHMratio->Draw("e2");
  hDYBratio->Draw("E2 same");
  hHMratio->SetTitle("H.-M.");
  hDYBratio->SetTitle("DYB");
  gPad->BuildLegend(0.318236,0.788293,0.698451,0.998049);
  hHMratio->SetTitle("Spectra Ratio");
  TLine l1(40,1,440,1);
  l1.SetLineWidth(3);
  l1.SetLineStyle(7);
  l1.Draw();
  hDYBratio->Draw("E2 same");

  //c5.SaveAs("DYB_vs_HM_ratio2.pdf");
  //c5.SaveAs("DYB_vs_HM_ratio2.png");
  //c5.SaveAs("DYB_vs_HM_ratio2.root");
}
