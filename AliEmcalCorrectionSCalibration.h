#ifndef ALIEMCALCORRECTIONSCCALIBRATION_H
#define ALIEMCALCORRECTIONSCCALIBRATION_H

#include <vector>
#include <TList.h>
#include <TH2D.h>

#include "AliEmcalCorrectionComponent.h"

class AliEmcalCorrectionSCCalibration : public  AliEmcalCorrectionComponent
{
public:
  AliEmcalCorrectionSCCalibration();
  virtual ~AliEmcalCorrectionSCCalibration();

  //Set up and run the task
  Bool_t Initialize();
  void UserCreateOutputObjects(void);
  Bool_t Run();
  Boot_t CheckIfRunChanged();
  
protected:
  TList *_list; //!
  
  TH2D *_histogram_cell_id_amplitude; //!
  TH2D *_histogram_cell_id_amplitude_scale; //!
  
private:
  size_t _ncell; //!
  size_t _nsm_ieta; //!
  size_t _nsm_iphi; //!

  //Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEnergy> reg;

  ClassDef(AliEmcalCorrectionCellEnergy, 1); // EMCal cell energy correction component
};

#endif /* ALIEMCALCORRECTIONSCCALIBRATION*/
