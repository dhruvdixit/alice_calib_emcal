#include <vector>
#include <iostream>
#include <fstream>

#include "TGeoManager.h"
#include "TH2D.h"
#include "TFile.h"
#include "TH1.h"
#include "TObjArray.h"

//#include "AliEMCALGeometry.h"
//#include "AliEMCALRecoUtils.h"
//#include "AliOADBContainer.h"
//#include "AliCalorimeterUtils.h"


using namespace std;

int getSM_ETA_PHI(int n)
{
  int sm, n0, nphi, ieta, iphi, n1;
  sm = n0 = nphi = ieta = iphi = n1 = 0;
  if(n < 11520)
    {
      sm = n / 1152;
      n0 = sm * 1152;
      nphi = 24;
    }
  else if (n < 12288)
    {
      sm = 10 + (n - 11520) / 384;
      n0 = 11520 + (sm - 10) * 384;
      nphi = 8;
    }
  else if (n < 16896)
    {
      sm = 12 + (n - 12288) / 768;
      n0 = 12288 + (sm - 12) * 768;
      nphi = 24;
    }
  else
    {
      sm = 18 + (n - 16896) / 384;
      n0 = 16896 + (sm - 18) * 384;
      nphi = 8;
    }
  n1 = n - n0;
  ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
  iphi = (n1 / 2) % nphi;
  int ret = iphi+ieta*100+sm*10000;
  return ret;
  
}

void makingCoeffComparison()//TString fBasePath)
{
  gStyle->SetOptStat(0);
  
  //......................................................
  //Get Coeff text files
  //TString fBasePath="/global/homes/d/ddixit/alice_calib_emcal";

  TH2F* hCoeffHist15o17q = new TH2F("hCoeffHist15o17q", "hCoeffHist15o17q", 100, 0.5, 1.5, 100, 0.5, 1.5);
  TH2F* hCoeffHist15o17o = new TH2F("hCoeffHist15o17o", "hCoeffHist15o17o", 100, 0.5, 1.5, 100, 0.5, 1.5);
  TH2F* hCoeffHist17q17o = new TH2F("hCoeffHist17q17o", "hCoeffHist17q17o", 100, 0.5, 1.5, 100, 0.5, 1.5);

  TH1F* hCoeffDiff15o17q = new TH1F("hCoeffDiff15o171", "hCoeffDiff15o17q", 230, -15.5, 7.5);
  TH1F* hCoeffDiff15o17o = new TH1F("hCoeffDiff15o17o", "hCoeffDiff15o17o", 230, -15.5, 7.5);
  TH1F* hCoeffDiff17q17o = new TH1F("hCoeffDiff17q17o", "hCoeffDiff17q17o", 230, -15.5, 7.5);


  TH1F* hCoeff15o = new TH1F("hCoeff15o", "hCoeff15o", 400, -0.5, 3.5);
  TH1F* hCoeff17q = new TH1F("hCoeff17q", "hCoeff17q", 400, -0.5, 3.5);
  TH1F* hCoeff17o = new TH1F("hCoeff17o", "hCoeff17o", 400, -0.5, 3.5);
  TH1F* hCoeff18m = new TH1F("hCoeff18m", "hCoeff18m", 400, -0.5, 3.5); 

  const int ncell = 17664;
  double coeff15o[ncell] = {0.0};
  double coeff17q[ncell] = {0.0};
  double coeff17o[ncell] = {0.0};
  double coeff18m[ncell] = {0.0};

  int cellID, badChn, badChnMine;
  double a, chiA, chiB;
  
  //cout << "hello" << endl;
  ifstream inFile15o("lhc15o_all_new_chiBeforeAfter.txt");
  while(inFile15o >> cellID >>  badChn >>  badChnMine >> a >>  chiA >>  chiB){
    //cout << cellID << "\t" << a << endl;
    coeff15o[cellID] = a;
  }
  //cout << "bye" << endl;

  ifstream inFile17q("lhc17q_all_coeff.txt");
  while(inFile17q >> cellID >> a){
    //cout << cellID << "\t" << a << endl;
    coeff17q[cellID] = a;
  }

  ifstream inFile17o("lhc17o_coeff.txt");
  while(inFile17o >> cellID >> a){
    //cout << cellID << "\t" << a << endl;
    coeff17o[cellID] = a;
  }

  ifstream inFile18m("lhc18m_coeff.txt");
  while(inFile18m >> cellID >> a){
    //cout << cellID << "\t" << a << endl;
    coeff18m[cellID] = a;
  }

  for(int i = 0; i < ncell; i++){
    //cout << coeff15o[i] << "\t" << coeff17q[i] << endl;
    if ((coeff15o[i] == 0) || (coeff17q[i] == 0) || (coeff17o[i] == 0)) continue;
    hCoeffHist15o17q->Fill(coeff15o[i], coeff17q[i]);
    hCoeffHist15o17o->Fill(coeff15o[i], coeff17o[i]);
    hCoeffHist17q17o->Fill(coeff17q[i], coeff17o[i]);

    hCoeffDiff15o17q->Fill(coeff15o[i]-coeff17q[i]);
    hCoeffDiff15o17o->Fill(coeff15o[i]-coeff17o[i]);
    hCoeffDiff17q17o->Fill(coeff17q[i]-coeff17o[i]);
  }
  for(int i = 0; i < ncell; i++){
    if (coeff15o[i] != 0)
      hCoeff15o->Fill(coeff15o[i]);
    if (coeff17q[i] != 0)
      hCoeff17q->Fill(coeff17q[i]);
    if (coeff17o[i] != 0)
      hCoeff17o->Fill(coeff17o[i]);
    if (coeff18m[i] != 0)
      hCoeff18m->Fill(coeff18m[i]);
  }

  hCoeffHist15o17q->SetTitle("; 15o coefficients; 17q coefficients");
  hCoeffHist15o17o->SetTitle("; 15o coefficients; 17o coefficients");
  hCoeffHist17q17o->SetTitle("; 17q coefficients; 17o coefficients");
  
  hCoeffDiff15o17q->SetTitle("; 15o coefficients - 17q coefficients;counts");hCoeffDiff15o17q->SetLineColor(kRed);hCoeffDiff15o17q->SetMarkerColor(kRed);hCoeffDiff15o17q->SetMarkerStyle(kFullCircle);
  hCoeffDiff15o17o->SetTitle("; 15o coefficients - 17o coefficients;counts");hCoeffDiff15o17o->SetLineColor(kBlue);hCoeffDiff15o17o->SetMarkerColor(kBlue);hCoeffDiff15o17o->SetMarkerStyle(kFullCircle);
  hCoeffDiff17q17o->SetTitle("; 17q coefficients - 17o coefficients;counts");hCoeffDiff17q17o->SetLineColor(kGreen);hCoeffDiff17q17o->SetMarkerColor(kGreen);hCoeffDiff17q17o->SetMarkerStyle(kFullCircle);
  
  hCoeff15o->SetTitle("; a (lhc15o);counts");hCoeff15o->SetLineColor(kRed);hCoeff15o->SetMarkerColor(kRed);hCoeff15o->SetMarkerStyle(kFullCircle);
  hCoeff17q->SetTitle("; a (lhc17q);counts");hCoeff17q->SetLineColor(kBlue);hCoeff17q->SetMarkerColor(kBlue);hCoeff17q->SetMarkerStyle(kFullCircle);
  hCoeff17o->SetTitle("; a (lhc17o);counts");hCoeff17o->SetLineColor(kGreen);hCoeff17o->SetMarkerColor(kGreen);hCoeff17o->SetMarkerStyle(kFullCircle);
  hCoeff18m->SetTitle("; a (lhc18m);counts");hCoeff18m->SetLineColor(kMagenta);hCoeff18m->SetMarkerColor(kMagenta);hCoeff18m->SetMarkerStyle(kFullCircle);
  
  TCanvas* cHist15o17q = new TCanvas("cHist15o17q", "cHist15o17q", 800, 600);
  cHist15o17q ->SetLogz();
  hCoeffHist15o17q->Draw("colz");
  TCanvas* cHist15o17o = new TCanvas("cHist15o17o", "cHist15o17o", 800, 600);
  cHist15o17o ->SetLogz();
  hCoeffHist15o17o->Draw("colz");
  TCanvas* cHist17q17o = new TCanvas("cHist17q17o", "cHist17q17o", 800, 600);
  cHist17q17o ->SetLogz();
  hCoeffHist17q17o->Draw("colz");

  TCanvas* cDiff = new TCanvas("cDiff", "cDiff", 600, 600);
  cDiff->SetLogy();
  hCoeffDiff15o17q->SetTitle(";difference between periods;counts");;
  hCoeffDiff15o17q->Draw("");
  hCoeffDiff15o17o->Draw("same");
  hCoeffDiff17q17o->Draw("same");
  TLegend* lDiff = new TLegend(0.6, 0.6, 0.8, 0.8);
  lDiff->AddEntry(hCoeffDiff15o17q, "15o - 17q");
  lDiff->AddEntry(hCoeffDiff15o17o, "15o - 17o");
  lDiff->AddEntry(hCoeffDiff17q17o, "17q - 17o");
  lDiff->Draw("same");
  
  cout << "Integral (-12, -0.5): " << hCoeffDiff15o17q->Integral(hCoeffDiff15o17q->GetXaxis()->FindBin(-12), (hCoeffDiff15o17q->GetXaxis()->FindBin(-0.5))) << endl;
  cout << "Integral (-0.5, 0.5): " << hCoeffDiff15o17q->Integral(hCoeffDiff15o17q->GetXaxis()->FindBin(-0.5), (hCoeffDiff15o17q->GetXaxis()->FindBin(0.5))) << endl;
  cout << "Integral (0.5, 6): " << hCoeffDiff15o17q->Integral(hCoeffDiff15o17q->GetXaxis()->FindBin(0.5), (hCoeffDiff15o17q->GetXaxis()->FindBin(6.0))) << endl;

  /*TCanvas* c15o = new TCanvas("c15o", "c15o", 800, 600);
  c15o->SetLogy();
  hCoeff15o->Draw("");

  TCanvas* c17q = new TCanvas("c17q", "c17q", 800, 600);
  c17q->SetLogy();
  hCoeff17q->Draw("");

  TCanvas* c17o = new TCanvas("c17o", "c17o", 800, 600);
  c17o->SetLogy();
  hCoeff17o->Draw("");*/

  TCanvas* cCoeffs = new TCanvas("cCoeffs", "cCoeffs", 800, 600);
  cCoeffs->SetLogy();
  hCoeff15o->SetTitle(";calibration coefficient (a); counts");
  hCoeff15o->Draw("");
  hCoeff17q->Draw("same");
  hCoeff17o->Draw("same");
  hCoeff18m->Draw("same");
  
  TLegend* lCoeffs = new TLegend(0.6, 0.6, 0.8, 0.8);
  lCoeffs->SetBorderSize(0);
  lCoeffs->AddEntry(hCoeff15o, "15o, 5 TeV Pb-Pb");
  lCoeffs->AddEntry(hCoeff17o, "17o, 13 TeV pp");
  lCoeffs->AddEntry(hCoeff17q, "17q, 5 TeV pp");
  lCoeffs->AddEntry(hCoeff18m, "18m, 13 TeV pp");
  lCoeffs->Draw("same");
  
}//end makeEtaPhiPlots
