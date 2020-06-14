#include <vector>
#include <fstream>
#include <string>
#include "TGeoManager.h"

using namespace std;

void makingSCOADB(string inputFile, Int_t runNumber)
{

  
  const Int_t nSM          = 20;    //Number of supermodules
  const Int_t fNMaxCols    = 48;    //eta direction
  const Int_t fNMaxRows    = 24;    //phi direction
  const Int_t nCell        = 17664; //Number of cells
  //Int_t nRuns        = RunIdVec.size();
  
  //Int_t fNMaxColsAbs = 2*fNMaxCols;
  //Int_t fNMaxRowsAbs = Int_t (nSM/2)*fNMaxRows; //multiply by number of supermodules

  cout << "Make histos" << endl;
  
  TString histoName[nSM];
  for(int i = 0; i < nSM; i++)
    {
      Long_t n = i;
      TString s = "SC_CalibCoeff_SM";
      histoName[i] = s+n;
    }
  TH2F* plot2D_SC_hists[nSM];
  for(int i = 0; i < nSM; i++)
    {
      plot2D_SC_hists[i] = new TH2F(histoName[i],histoName[i],fNMaxCols+1,-0.5,fNMaxCols+0.5, fNMaxRows+1,-0.5,fNMaxRows+0.5);
      plot2D_SC_hists[i]->GetXaxis()->SetTitle("cell column (#eta direction)");
      plot2D_SC_hists[i]->GetYaxis()->SetTitle("cell row (#phi direction)");
  
    }

  cout << "read from file" << endl;

  ifstream myfile;
  myfile.open(inputFile.c_str());

  int id, bad1, bad2;
  id = bad1 = bad2 = 0;
  double coeff, chiA, chiB;
  coeff = chiA = chiB =0.0;

  if(!myfile.is_open())
    {
      cout << "File did not open" << endl;
      return;
    }
  //for(int cellID = 0; cellID < nCell; cellID++)
  while(myfile >> id >> bad1 >> bad2 >> coeff >> chiA >> chiB)
    {      
      //cout << id << "\t" << bad1 << "\t" << bad2 << "\t" << coeff << "\t" << chiA << "\t" << chiB << endl;;
      if(bad1 || bad2 || chiA > 10.0 || coeff > 1.5)
	coeff = 0;
      
      int locInfo = getSM_ETA_PHI(id);
      int ism = locInfo/10000;
      int ieta = (locInfo%10000)/100;
      int iphi = (locInfo%10000)%100;
      //cout << ism << "\t" << ieta << "\t" << iphi << "\t" << coeff << endl;
      plot2D_SC_hists[ism]->SetBinContent(ieta,iphi,coeff);
    }

  
  myfile.close();

  cout << "done" << endl;
  TFile* f1 = new TFile("EMCAL_SC_calib_coeff.root","UPDATE");
  TDirectory* runGroup = f1->mkdir(runPeriod);
  runGroup->cd();
  for(int i = 0; i < nSM; i++)
    {
      plot2D_SC_hists[i]->Write(histoName[i]);
    }
  f1->Close();
  delete f1;
  

}


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
