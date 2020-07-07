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

void makeEtaPhiPlots()//TString fBasePath)
{
  gStyle->SetOptStat(0);


  //......................................................
  //..Get the energy histograms
  TString fBasePath="/global/homes/d/ddixit/alice_calib_emcal";
  //TFile* fInputFile = new TFile(Form("%s/PbPb_globallyGood/AnalysisResults_0031.root", fBasePath.Data()), "READ");
  //TFile* fInputFile = new TFile(Form("%s/PbPb_globallyGood/AnalysisResults_list2_new.root", fBasePath.Data()), "READ");
  //TFile* fInputFile = new TFile(Form("%s/old_runGroups/AnalysisResults_allBut2.root", fBasePath.Data()), "READ");
  //TFile* fInputFile = new TFile(Form("%s/pp_outputs/AnalysisResults_17q.root", fBasePath.Data()), "READ");
  //TFile* fInputFile = new TFile(Form("%s/outputdir_17o/AnalysisResults_17o.root", fBasePath.Data()), "READ");
  TFile* fInputFile = new TFile(Form("%s/LHC16h/AnalysisResults.root", fBasePath.Data()), "READ");
  
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

  cout << "Make histos" << endl;
  
  TString histoNameFullRange[nSM];
  TString histoNameSCRange[nSM];
  TString histoNameSmallRange[nSM];
  for(int i = 0; i < nSM; i++){
    Long_t n = i;
    TString s = "SC_CellEnergy_SM";
    TString ssc = "SC_CellEnergySCRange_SM";
    TString small = "SC_CellEnergySmallRange_SM";
    histoNameFullRange[i] = s+n;
    histoNameSCRange[i] = ssc+n;
    histoNameSmallRange[i] = small+n;
    
  }//hist names loop


  TH2F* plot2D_etaPhi_energy[nSM];
  TH2F* plot2D_etaPhi_energySCRange[nSM];
  TH2F* plot2D_etaPhi_energySmallRange[nSM];

  for(int i = 0; i < nSM; i++){
    //Full integral
    plot2D_etaPhi_energy[i] = new TH2F(histoNameFullRange[i],histoNameFullRange[i],fNMaxCols+1,-0.5,fNMaxCols+0.5, fNMaxRows+1,-0.5,fNMaxRows+0.5);
    plot2D_etaPhi_energy[i]->SetTitle(Form("SM %i Integral Range: Full;cell column (#eta direction);cell row (#phi direction); Cell Energy [GeV]", i));

    //integral between 1.5 - 5 GeV
    plot2D_etaPhi_energySCRange[i] = new TH2F(histoNameSCRange[i],histoNameSCRange[i],fNMaxCols+1,-0.5,fNMaxCols+0.5, fNMaxRows+1,-0.5,fNMaxRows+0.5);
    plot2D_etaPhi_energySCRange[i]->SetTitle(Form("SM %i Integral Range: 1.5-5 GeV;cell column (#eta direction);cell row (#phi direction); Cell Energy [GeV]", i));
	
    //integral between 0 - 0.5 GeV
    plot2D_etaPhi_energySmallRange[i] = new TH2F(histoNameSmallRange[i],histoNameSmallRange[i],fNMaxCols+1,-0.5,fNMaxCols+0.5, fNMaxRows+1,-0.5,fNMaxRows+0.5);
    plot2D_etaPhi_energySmallRange[i]->SetTitle(Form("SM %i Integral Range: 0.0-0.5 GeV;cell column (#eta direction);cell row (#phi direction); Cell Energy [GeV]", i));
  }//
  
  for(int i = 0; i < nCell; i++){      

    //location
    int locInfo = getSM_ETA_PHI(i);
    int ism = locInfo/10000;
    int ieta = (locInfo%10000)/100;
    int iphi = (locInfo%10000)%100;

    //energy
    TH1D* fCellAmplitude = fCellIDAmplitudeHist->ProjectionY("fCellAmplitude", i+1, i+1);
    double energy = fCellAmplitude->Integral();
    plot2D_etaPhi_energy[ism]->SetBinContent(ieta,iphi,energy);
    
    int energyLowBinSC = fCellIDAmplitudeHist->GetYaxis()->FindBin(1.5);//1.5 GeV Bin
    int energyHighBinSC = fCellIDAmplitudeHist->GetYaxis()->FindBin(5.0) - 1;//1.5 GeV Bin
    double energySC = fCellAmplitude->Integral(energyLowBinSC, energyHighBinSC);
    plot2D_etaPhi_energySCRange[ism]->SetBinContent(ieta,iphi,energySC);

    int energyLowBinSmall = fCellIDAmplitudeHist->GetYaxis()->FindBin(0.0);//0.0 GeV Bin
    int energyHighBinSmall = fCellIDAmplitudeHist->GetYaxis()->FindBin(0.5) - 1;//0.5 GeV Bin
    double energySmall = fCellAmplitude->Integral(energyLowBinSmall, energyHighBinSmall);
    plot2D_etaPhi_energySmallRange[ism]->SetBinContent(ieta,iphi,energySmall);
    
  }

  /*TFile* f1 = new TFile("zebraPattern_17q.root","RECREATE");
  for(int i = 0; i < nSM; i++)
    {
      plot2D_etaPhi_energy[i]->Write(histoNameFullRange[i]);
      plot2D_etaPhi_energySCRange[i]->Write(histoNameSCRange[i]);
      plot2D_etaPhi_energySmallRange[i]->Write(histoNameSmallRange[i]);
    }
  f1->Close();
  delete f1;//*/

  TCanvas* emcalSM6 = new TCanvas("emcalSM6", "emcalSM6", 1800, 600);
  emcalSM6->Divide(3, 1);
  emcalSM6->cd(1);
  plot2D_etaPhi_energy[5]->Draw("colz");
  emcalSM6->SetLogz();
  emcalSM6->cd(2);
  emcalSM6->SetLogz();
  plot2D_etaPhi_energySCRange[5]->Draw("colz");
  emcalSM6->cd(3);
  emcalSM6->SetLogz();
  plot2D_etaPhi_energySmallRange[5]->Draw("colz");
  
}//end makeEtaPhiPlots
