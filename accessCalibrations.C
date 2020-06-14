void accessCalibrations(TString runPeriod, int cellID)
{
  TFile* f1 = new TFile("EMCAL_SC_calib_coeff.root", "READ");
  f1->cd(runPeriod);
  int locInfo = getSM_ETA_PHI(cellID);
  Int_t ism = locInfo/10000;
  int ieta = (locInfo%10000)/100;
  int iphi = (locInfo%10000)%100;
  TString s = "";
  
}

double GettingCellAmplitude(int cell_id)
{
  AliVEvent* event = InputEvent();
  AliVCaloCells* emcal_cell = event->GetEMCALCells();
  double amp = emcal_cell->GetCellAmplitude(cell_id);
  return amp;
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
