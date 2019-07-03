//!  EMCal calib coefficient map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/EMCALSCCalibrations.root
*/

#include "TGeoManager.h"
#include "TH2D.h"
#include "TFile.h"

#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliOADBContainer.h"

#include <utility>
#include <fstream>

/*******************************************************************
 *  NOTE: Function to fill the 2D bad channel histograms           *
 *             using EMCAL Geo for cellID->Eta/Phi                 *
 *******************************************************************/
void Fill_histo(std::vector<std::pair<Int_t, Double_t>> vecChannels,TH2D **hSCCoeff,AliEMCALGeometry   *fEMCALGeo,AliEMCALRecoUtils  *fEMCALRecoUtils, Int_t fill=1){
    Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
    for(Int_t i=0; i<(Int_t) vecChannels.size(); i++){
      std::pair<Int_t, Double_t> chnCoeff = vecChannels.at(i);
      Int_t chn = chnCoeff.first;
      Double_t a = chnCoeff.second;
      fEMCALGeo->GetCellIndex(chn,nSupMod,nModule,nIphi,nIeta);
      fEMCALGeo->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta,iphi,ieta);
      fEMCALRecoUtils->SetEMCALChannelStatus(nSupMod, ieta, iphi);
      //cout<<endl<<"ID:"<<vecChannels.at(i).Atoi()<<" NSupMod: "<<nSupMod<<" ieta: "<<ieta<<" iphi:"<<iphi<<" fill:"<<fill<<endl; 
      hSCCoeff[nSupMod]->SetBinContent(ieta,iphi,a);
    }
}

/*******************************************************************
 *  NOTE: Sorting function to sort the final OADB file             *
 *                  by ascending runnumber                         *
 *******************************************************************/
void sortOutput(const char *fileNameOADB=""){

    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCALSingleChannelCoefficient");
    con->SetName("Old"); 

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALSingleChannelCoefficient");
    // First entry needs to be added before sorting loop
    con2->AddDefaultObject(con->GetObjectByIndex(0));
    con2->AppendObject(con->GetObjectByIndex(0),con->LowerLimit(0),con->UpperLimit(0));
    TString strTemp                             = "";

    // sorting magic happens here
    for(int i=1;i<con->GetNumberOfEntries();i++){
        largerthan                              = con2->UpperLimit(con2->GetNumberOfEntries()-1);
        currentvalue                            = -1;
        indexAdd                                = 0;
        for(int j=0;j<con->GetNumberOfEntries();j++){
            if(con->UpperLimit(j)<=largerthan) 
                continue;
            if(currentvalue < 0){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
            if(con->UpperLimit(j)<currentvalue){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
        }
        con2->AddDefaultObject(con->GetObjectByIndex(indexAdd));
        con2->AppendObject(con->GetObjectByIndex(indexAdd),con->LowerLimit(indexAdd),con->UpperLimit(indexAdd));
    }

    printf("\n\n");
    Int_t nentries2                             = con2->GetNumberOfEntries();
    for(int i=0;i<nentries2;i++){
        printf("\n Entry2 --> %d/%d -->",i,nentries2);
        printf("%d -- %d --> obj = %p , %s", con2->LowerLimit(i),con2->UpperLimit(i),con2->GetObjectByIndex(i),con2->GetObjectByIndex(i)->GetName());
    }
    printf("\n\n");

    con2->WriteToFile("EMCALSingleChannelCalibNEW.root");
}

/*******************************************************************
 *  NOTE: Function to read the bad channels from                   *
 *                  the given input file                           *
 *******************************************************************/
Bool_t readin(TString fileRuns, std::vector<pair<Int_t, Double_t>> &vec, Bool_t output)
{
  if(output) std::cout << Form("\nReading from %s...", fileRuns.Data()) << endl;
  fstream file;
  //TString fVar;
  Int_t cellID, cellBad, cellHot;
  Double_t a, chiA, chiB;
  Int_t totalN=0;
  file.open(fileRuns.Data(), ios::in);
  if(file.good())
    {
        file.seekg(0L, ios::beg);
	cellID = cellBad = cellHot = -9999;
	a = chiA = chiB = -9999.99;
	  while(!file.eof())
	  {
	    file >> cellID >> cellBad >> cellHot >> a >> chiA >> chiB;
	    if(cellBad || cellHot) a = 0;
	    if(totalN%100 == 0 && output) cout << cellID << "\t" << cellBad << "\t" << cellHot << "\t" << a << "\t" << chiA << "\t" << chiB << endl;
	    std::pair<Int_t, Double_t> chnCoeff(cellID, a);
	    vec.push_back(chnCoeff);
	    totalN++;
	    //cout << totalN << endl;
	    //if(totalN > 50) return kFALSE;
	  }
    }
  if(output) std::cout << "\"" << endl;
  file.close();
  if(totalN>0){
    return kTRUE;
  }
  else return kFALSE;
}
/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameCalibCoeff,Int_t lowRun, Int_t highRun, Int_t updateExistingMap = 0){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

    AliEMCALGeometry   *fEMCALGeo               = new AliEMCALGeometry();
    AliEMCALRecoUtils  *fEMCALRecoUtils         = new AliEMCALRecoUtils();

    fEMCALGeo                                   = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
    const int nSM=20;

    // load OADB file
    TFile *f                                    = TFile::Open(fileNameOADB);
    f->ls();
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCALSingleChannelCoefficient");
    con->SetName("Old"); 

    // make the output container
    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALSingleChannelCoefficient");

    // List of brand new arrays (including overhead arrays for updateExistingMap == 2)
    TObjArray arrayAdd(nSM);
    arrayAdd.SetName(arrName.Data());
    TObjArray arrayAddOverheadLow(nSM);
    arrayAddOverheadLow.SetName(Form("%s_OHLow",arrName.Data()));
    TObjArray arrayAddOverheadHigh(nSM);
    arrayAddOverheadHigh.SetName(Form("%s_OHHigh",arrName.Data()));

    // load OADB file to have access to existing maps
    TFile *foadb                                = TFile::Open(fileNameOADB);
    AliOADBContainer *conAP                     = (AliOADBContainer*)foadb->Get("AliEMCALSingleChannelCoefficient");
    TObjArray *arrayBClow                       = (TObjArray*)conAP->GetObject(lowRun);
    TObjArray *arrayBChigh                      = (TObjArray*)conAP->GetObject(highRun);

    TH2D *hSCCoeff[nSM];
    TH2D *hSCCoeffOverhead[nSM];
    char nameSM[50];

    // initialize histograms in array - either load existing map or create new, empty map
    for(int i=0;i<nSM;i++){
        sprintf(nameSM,"EMCALSCCalibMap_Mod%d",i);
        hSCCoeff[i]                              = NULL;
        hSCCoeffOverhead[i]                      = NULL;
        if(updateExistingMap){
            cout << "Initializing " << nameSM << " from current AliPhysics OADB file ... ";
            if(arrayBClow){
                hSCCoeff[i]                      = (TH2D*)arrayBClow->FindObject(nameSM);
                if(hSCCoeff[i] && updateExistingMap>1)
                    hSCCoeffOverhead[i]          = (TH2D*)hSCCoeff[i]->Clone(nameSM);
            } else if (arrayBChigh){
                hSCCoeff[i]                      = (TH2D*)arrayBChigh->FindObject(nameSM);
                if(hSCCoeff[i] && updateExistingMap>1)
                    hSCCoeffOverhead[i]          = (TH2D*)hSCCoeff[i]->Clone(nameSM);
            }
            if(hSCCoeff[i])
                cout << "and found it!" << endl;
        }
        if(!updateExistingMap || !hSCCoeff[i]){
            cout << "File not found! Creating new map for " << nameSM << endl;
            hSCCoeff[i]                          = new TH2D(nameSM,nameSM,48,0,48,24,0,24);
        }
    }

    // Read in the list of dead channels from the dead cell file
    //std::vector<TString> vecCalibCoefficients;
    std::vector<pair<Int_t, Double_t>> vecChnCoeffPair;
    if(!readin(fileNameCalibCoeff, vecChnCoeffPair,kTRUE)) {
        cout << "No coefficients could be found in file: " << fileNameCalibCoeff.Data() << endl;
    }


    // Fill the SC map with calibration coefficients
    Fill_histo(vecChnCoeffPair,hSCCoeff,fEMCALGeo,fEMCALRecoUtils,1);

    // add histograms to a new array
    for (Int_t mod=0;mod<nSM;mod++){
        arrayAdd.Add(hSCCoeff[mod]);
        if( updateExistingMap>1 && hSCCoeffOverhead[mod]){
            arrayAddOverheadLow.Add(hSCCoeffOverhead[mod]);
            arrayAddOverheadHigh.Add(hSCCoeffOverhead[mod]);
        }
    }

    // check old OADB file for existing BC maps in given run range and remove that old BC map
    Int_t replacedContainerLowLimit             = -1;
    Int_t replacedContainerHighLimit            = -1;
    for(int i=0;i<con->GetNumberOfEntries();i++){
        if (lowRun >= con->LowerLimit(i) && lowRun <= con->UpperLimit(i)){
            printf("\n!!! Not adding index %d for runrange %d--%d as low run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),lowRun);
            replacedContainerLowLimit           = con->LowerLimit(i);
            replacedContainerHighLimit          = con->UpperLimit(i);
        } else if (highRun >= con->LowerLimit(i) && highRun <= con->UpperLimit(i)){
            printf("\n!!! Not adding index %d for runrange %d--%d as high run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),highRun);
            replacedContainerLowLimit           = con->LowerLimit(i);
            replacedContainerHighLimit          = con->UpperLimit(i);
        }else{
            con2->AddDefaultObject(con->GetObjectByIndex(i));
            con2->AppendObject(con->GetObjectByIndex(i),con->LowerLimit(i),con->UpperLimit(i));
        }
    }
    // add new BC maps at the end of the new OADB file (will be sorted later)
    if( updateExistingMap>1 && replacedContainerLowLimit>0 && lowRun > replacedContainerLowLimit){
        con2->AddDefaultObject(&arrayAddOverheadLow);
        con2->AppendObject(&arrayAddOverheadLow, replacedContainerLowLimit,lowRun-1);
    }
    con2->AddDefaultObject(&arrayAdd);
    con2->AppendObject(&arrayAdd, lowRun,highRun);
    if(updateExistingMap>1 && replacedContainerHighLimit>0 && highRun < replacedContainerHighLimit){
        con2->AddDefaultObject(&arrayAddOverheadHigh);
        con2->AppendObject(&arrayAddOverheadHigh, highRun+1,replacedContainerHighLimit);
    }

    // temporarilty save BC map file and rename as new input file
    con2->WriteToFile("tempBC.root");
    gSystem->Exec(Form("mv tempBC.root %s",fileNameOADB));

    printf("\nMaps have been successfully added!\n\n",arrName.Data());

}
/*******************************************************************
 *  NOTE: Main function which needs to be adjusted for new BC maps *
 *******************************************************************/
void UpdateEMCAL_SCOADB(const char *fileNameOADBAli="/home/dhruv/alice_calib_emcal/EMCALSingleChannelCalibrations.root")
{
    gSystem->Load("libOADB");  
    gSystem->Load("libEMCALbase");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libEMCALrec");

    const char *fileNameOADB                ="EMCALSingleChannelCalib_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAli,fileNameOADB));

    // update OADB file with dead, bad and warm cells
    // last parameter:
    //      0: write new map for given run range
    //      1: update existing map for given run range
    //      2: same as 1 but also save over/underhead parts

    // LHC15o - March 06 2019
    updateFile(fileNameOADB,"SingleChannelCoefficients_15o","/home/dhruv/alice_calib_emcal/LHC15o/LHC15o_OADBFile_Coeffs.txt",244824,246994,1);

    // the final output will be sorted by runnumber
    sortOutput(fileNameOADB);

}
