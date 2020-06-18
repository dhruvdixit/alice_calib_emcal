#include <TObjArray.h>
#include <TFile.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include <AliOADBContainer.h>

#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDHeader.h>
#include <AliVCluster.h>
#include <AliVCaloCells.h>


#include "AliEmcalCorrectionSCCalibration.h"



namespace {


	void to_sm_ieta_iphi(unsigned int &sm, unsigned int &ieta,
						 unsigned int &iphi, unsigned int n)
	{
		sm = n < 11520 ? n / 1152 :
			n < 12288 ? 10 + (n - 11520) / 384 :
			n < 16896 ? 12 + (n - 12288) / 768 :
			18 + (n - 16896) / 384;

		const unsigned int n0 =
			sm < 10 ? sm * 1152 :
			sm < 12 ? 11520 + (sm - 10) * 384 :
			sm < 18 ? 12288 + (sm - 12) * 768 :
			16896 + (sm - 18) * 384;
		const unsigned int n1 = n - n0;
		const unsigned int nphi =
			sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

		ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
		iphi = (n1 / 2) % nphi;
	}


}

ClassImp(AliEmcalCorrectionCellEnergy);

//Register the class with base class
RegisterCorrectionComponent<AliEmcalCorrectionSCCalibration> AliEmcalCorrectionSCCalibration::reg("AliEmcalCorrectionSCCalibration");

#def							\
		,	\
		
AliEmcalCorrectionSCCalibration::AliEmcalCorrectionSCCalibration():
  AliEmcalCorrectionSCCalibration("AliEmcalCorrectionSCCalibration")
  ,_list(NULL)
  ,_histogram_cell_id_amplitude(NULL)
  ,_histogram_cell_id_amplitude_scale(NULL)
  ,_ncell(17664), _nsm_ieta(864), _nsm_iphi(416)
{
}

AliEmcalCorrectionSCCalibration::~AliEmcalCorrectionSCCalibration(void)
{
}


Bool_t AliEmcalCorrectionSCCalibration::Initialize()
{
  AliEmcalCorrectionComponent::Initialize();
  AliWarning("Init EMCAL cell bad channel removal");
  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
  
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);
  
  return kTRUE;
}


void AliEmcalCorrectionSCCalibration::UserCreateOutputObjects(void)
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  if(fCreateHisto)
    {
      /////////////////////////////////////////////////////////////////
      
      _histogram_cell_id_amplitude = new TH2D("_histogram_cell_id_amplitude", "", _ncell, -0.5, _ncell - 0.5, 3 * lcm_range[6], 0.1, 100);
      set_log_axis(_histogram_cell_id_amplitude->GetYaxis());
      _histogram_cell_id_amplitude->GetXaxis()->SetTitle("E(GeV)");
      fOutput->Add(_histogram_cell_id_amplitude_scale);
      
      /////////////////////////////////////////////////////////////////
      
      _histogram_cell_id_amplitude_scale = new TH2D("_histogram_cell_id_amplitude_scale", "", _ncell, -0.5, _ncell - 0.5, 3 * lcm_range[6], 0.1, 100);
      set_log_axis(_histogram_cell_id_amplitude_scale->GetYaxis());
      _histogram_cell_id_amplitude_scale->GetXaxis()->SetTitle("E(GeV)");
      fOutput->Add(_histogram_cell_id_amplitude_scale);
    }
     
}

Double_t AliEmcalCorrectionSCCalibration::emcal_scale(const Int_t cell_id)
{

}


Boot_t AliEmcalCorrectionSCCalibration::Run()
{
  AliVEvent *event = InputEvent();
  AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);
  AliAODEvent *aod_event = dynamic_cast<AliAODEvent *>(event);
  
  const Int_t run_number = event->GetRunNumber();
  const Int_t bunch_crossing = event->GetBunchCrossNumber();
  
  
  TRefArray *calo_cluster = new TRefArray();
  
  event->GetEMCALClusters(calo_cluster);
  
  AliVCaloCells *emcal_cell = event->GetEMCALCells();
  
  for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++)
    {
      AliVCluster *c = (AliVCluster *)calo_cluster->At(i);
      
      for (Int_t j = 0; j < c->GetNCells(); j++)
	{
	  const Int_t cell_id = c->GetCellsAbsId()[j];
	  const Double_t amplitude = emcal_cell->GetCellAmplitude(cell_id);
	  const Double_t scale = emcal_scale(cell_id);

	  if(fCreateHistos)
	    {
	      _histogram_cell_id_amplitude->Fill(cell_id, amplitude);
	      _histogram_cell_id_amplitude_scale->Fill(cell_id, amplitude * scale);
	    }
	}
    }
  return kTrue;
}

AliEMCALRecoUtils *AliEmcalCorrectionSCCalibration::GetEMCALRecoUtils(void)
{
	return _reco_util;
}
