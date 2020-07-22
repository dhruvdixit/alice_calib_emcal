#include <vector>

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

void drawBeforeAfter()//TString fBasePath)
{
  gStyle->SetOptStat(0);


  //......................................................
  //..Get the energy histograms
  TString fBasePath="/home/ddixit/alice_calib_emcal";
  TFile* fInputFile = new TFile(Form("%s/globallyGoodRuns//AnalysisResults_all_new.root", fBasePath.Data()), "READ");
  
  if(!fInputFile){
    cout << "Couldn't find fInputFile" << endl;
    return;
  }
  TList* fList = dynamic_cast<TList *> (dynamic_cast<TDirectoryFile *>  (fInputFile->Get("AliAnalysisTaskCalibEmcal"))->Get("histogram"));
  if(!fList){
    cout << "Couldn't find fList" << endl;
    return;
  }

  TH2D* fCellIDAmplitudeHist=(TH2D*)(fList->FindObject("_histogram_cell_id_amplitude"));
  if(!fCellIDAmplitudeHist){
    cout << "Couldn't find fCellIDAmplitudeHist" << endl;
    return;
  }

  TCanvas* c1 = new TCanvas();
  c1->SetLogy();
  c1->SetLogz();
  fCellIDAmplitudeHist->SetTitle(";Cell ID; Cell Amplitude [GeV]");
  fCellIDAmplitudeHist->Draw("colz");
  //..build two dimensional histogram with values row vs. column
  const Int_t nSM          = 20;    //Number of supermodules
  const Int_t fNMaxCols    = 48;    //eta direction
  const Int_t fNMaxRows    = 24;    //phi direction
  const Int_t nCell        = 17664; //Number of cells

  
  //energy

  //cell 1
  TH1D* fCellAmplitude_cell1B = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude_cell1B", 1+1, 1+1);
  TH1D* fCellAmplitude_cell1A = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude_cell1A", 1+1, 1+1);
  fCellAmplitude_cell1A->SetLineColor(kGreen); fCellAmplitude_cell1A->SetMarkerColor(kGreen); fCellAmplitude_cell1A->SetMarkerStyle(kFullCircle);
  fCellAmplitude_cell1B->SetLineColor(kBlack); fCellAmplitude_cell1B->SetMarkerColor(kBlack); fCellAmplitude_cell1B->SetMarkerStyle(kFullSquare);
  double scale_cell1 = 1.024263;
  for(int i = 1; i < fCellAmplitude_cell1A->GetNbinsX()+1; i++){
    fCellAmplitude_cell1A->SetBinContent(i, 0);
  }
  for(int i = 1; i < fCellAmplitude_cell1B->GetNbinsX()+1; i++){
    double energy = fCellAmplitude_cell1B->GetBinCenter(i);
    double newEnergy = energy*scale_cell1;
    double content = fCellAmplitude_cell1B->GetBinContent(i);
    int newBin = fCellAmplitude_cell1A->GetXaxis()->FindBin(newEnergy);
    if(newBin) fCellAmplitude_cell1A->SetBinContent(newBin, content);
  }

  TCanvas* c_cell1 = new TCanvas("c_cell1", "c_cell1", 800, 600);
  fCellAmplitude_cell1B->SetTitle(";Cell Amplitude [GeV];");
  c_cell1->SetLogx();
  c_cell1->SetLogy();
  fCellAmplitude_cell1B->Draw("");
  fCellAmplitude_cell1A->Draw("same");

  //cell 103
  TH1D* fCellAmplitude_cell103B = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude_cell103B", 103+1, 103+1);
  TH1D* fCellAmplitude_cell103A = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude_cell103A", 103+1, 103+1);
  fCellAmplitude_cell103A->SetLineColor(kGreen); fCellAmplitude_cell103A->SetMarkerColor(kGreen); fCellAmplitude_cell103A->SetMarkerStyle(kFullCircle);
  fCellAmplitude_cell103B->SetLineColor(kBlack); fCellAmplitude_cell103B->SetMarkerColor(kBlack); fCellAmplitude_cell103B->SetMarkerStyle(kFullSquare);
  double scale_cell103 = 3.937201;
  for(int i = 1; i < fCellAmplitude_cell103A->GetNbinsX()+1; i++){
    fCellAmplitude_cell103A->SetBinContent(i, 0);
  }
  for(int i = 1; i < fCellAmplitude_cell103B->GetNbinsX()+1; i++){
    double energy = fCellAmplitude_cell103B->GetBinCenter(i);
    double newEnergy = energy*scale_cell103;
    double content = fCellAmplitude_cell103B->GetBinContent(i);
    int newBin = fCellAmplitude_cell103A->GetXaxis()->FindBin(newEnergy);
    if(newBin) fCellAmplitude_cell103A->SetBinContent(newBin, content);
  }

  TCanvas* c_cell103 = new TCanvas("c_cell103", "c_cell103", 800, 600);
  fCellAmplitude_cell103B->SetTitle(";Cell Amplitude [GeV];");
  c_cell103->SetLogx();
  c_cell103->SetLogy();
  fCellAmplitude_cell103B->Draw("");
  fCellAmplitude_cell103A->Draw("same");

  //cell 594
  TH1D* fCellAmplitude_cell594B = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude_cell594B", 594+1, 594+1);
  TH1D* fCellAmplitude_cell594A = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude_cell594A", 594+1, 594+1);
  fCellAmplitude_cell594A->SetLineColor(kGreen); fCellAmplitude_cell594A->SetMarkerColor(kGreen); fCellAmplitude_cell594A->SetMarkerStyle(kFullCircle);
  fCellAmplitude_cell594B->SetLineColor(kBlack); fCellAmplitude_cell594B->SetMarkerColor(kBlack); fCellAmplitude_cell594B->SetMarkerStyle(kFullSquare);
  double scale_cell594 = 1.251169;
  for(int i = 1; i < fCellAmplitude_cell594A->GetNbinsX()+1; i++){
    fCellAmplitude_cell594A->SetBinContent(i, 0);
  }
  for(int i = 1; i < fCellAmplitude_cell594B->GetNbinsX()+1; i++){
    double energy = fCellAmplitude_cell594B->GetBinCenter(i);
    double newEnergy = energy*scale_cell594;
    double content = fCellAmplitude_cell594B->GetBinContent(i);
    int newBin = fCellAmplitude_cell594A->GetXaxis()->FindBin(newEnergy);
    if(newBin) fCellAmplitude_cell594A->SetBinContent(newBin, content);
  }

  TCanvas* c_cell594 = new TCanvas("c_cell594", "c_cell594", 800, 600);
  fCellAmplitude_cell594B->SetTitle(";Cell Amplitude [GeV];");
  c_cell594->SetLogx();
  c_cell594->SetLogy();
  fCellAmplitude_cell594B->Draw("");
  fCellAmplitude_cell594A->Draw("same");

    
  
}//end makeEtaPhiPlots
