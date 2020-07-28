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

void makeEtaPhiPlotsCoeffs()//TString fBasePath)
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
  //TFile* fInputFile = new TFile(Form("%s/LHC16h/AnalysisResults.root", fBasePath.Data()), "READ");
  
  if(!fInputFile){
    cout << "Couldn't find fInputFile" << endl;
    return;
  }

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
  
  //..build two dimensional histogram with values row vs. column
  const Int_t nSM          = 20;    //Number of supermodules
  const Int_t fNMaxCols    = 48;    //eta direction
  const Int_t fNMaxRows    = 24;    //phi direction
  const Int_t nCell        = 17664; //Number of cells

  
  TString histoNameFullRange[nSM];
  for(int i = 0; i < nSM; i++){
    Long_t n = i;
    TString s = "SC_Coeffs_SM";
    histoNameFullRange[i] = s+n;
    
  }//hist names loop


  TH2F* plot2D_etaPhi_energy[nSM];

  for(int i = 0; i < nSM; i++){
    //Full integral
    plot2D_etaPhi_energy[i] = new TH2F(histoNameFullRange[i],histoNameFullRange[i],fNMaxCols+1,-0.5,fNMaxCols+0.5, fNMaxRows+1,-0.5,fNMaxRows+0.5);
    plot2D_etaPhi_energy[i]->SetTitle(Form("SM %i Integral Range: Full;cell column (#eta direction);cell row (#phi direction); Cell Energy [GeV]", i));

  }//
  
  for(int i = 0; i < nCell; i++){      

    //location
    int locInfo = getSM_ETA_PHI(i);
    int ism = locInfo/10000;
    int ieta = (locInfo%10000)/100;
    int iphi = (locInfo%10000)%100;

    //energy
    plot2D_etaPhi_energy[ism]->SetBinContent(ieta,iphi,energy);    
  }

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
