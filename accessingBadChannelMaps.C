#include <vector>

using namespace std;

void accessingBadChannelMaps()//TString fBasePath)
{
  gStyle->SetOptStat(0);
  //......................................................
  //..Get the OADB information
  TString fBasePath="/home/ddixit/localOADB/OADB/EMCAL";
  
  AliOADBContainer *cont=new AliOADBContainer("");
  cont->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");
  
  //......................................................
  //..Initialize EMCal/DCal geometry
  vector<int> RunIdVec;
  RunIdVec.push_back(280285);
  cout << RunIdVec.at(0) << endl;
  AliCalorimeterUtils* fCaloUtils = new AliCalorimeterUtils();
  //..Create a dummy event for the CaloUtils
  AliAODEvent* aod = new AliAODEvent();
  fCaloUtils->SetRunNumber(RunIdVec.at(0));
  fCaloUtils->AccessGeometry(aod);
  AliEMCALGeometry * geom = fCaloUtils->GetEMCALGeometry();
  
  
  //.......................................................
  //..build two dimensional histogram with values row vs. column
  //..with info from OADB
  Int_t nSM          = 20;  //Number of supermodules
  Int_t fNMaxCols    = 48;  //eta direction
  Int_t fNMaxRows    = 24;  //phi direction
  Int_t nRuns        = RunIdVec.size();
  
  Int_t fNMaxColsAbs = 2*fNMaxCols;
  Int_t fNMaxRowsAbs = Int_t (nSM/2)*fNMaxRows; //multiply by number of supermodules

  TH1C* h[20];

  //Bad channel ID vector
  vector<int> deadChn;
  vector<int> warmChn;
  vector<int> badChn;
  
  TString histoName;
  histoName = Form("2DChannelMap_Flag_Bad");
  TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
  //*OADB  looks like Marcels figures*/    TH2F *plot2D_Bad_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
  plot2D_Bad_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
  plot2D_Bad_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");
  histoName = Form("2DChannelMap_Flag_Dead");
  TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs+1,-0.5,fNMaxColsAbs+0.5, fNMaxRowsAbs+1,-0.5,fNMaxRowsAbs+0.5);
  //*OADB looks like Marcels figures*/    TH2F *plot2D_Dead_OADB = new TH2F(histoName,histoName,fNMaxColsAbs,0,fNMaxColsAbs, fNMaxRowsAbs,0,fNMaxRowsAbs);
  plot2D_Dead_OADB->GetXaxis()->SetTitle("cell column (#eta direction)");
  plot2D_Dead_OADB->GetYaxis()->SetTitle("cell row (#phi direction)");


  for(Int_t iRun = 0; iRun < nRuns; iRun++)
    {
      //cout<<"------ run "<<RunIdVec.at(iRun)<<endl;
      if(iRun%5==0)cout<<"."<<flush;
      if(iRun%20==0)cout<<"Run No."<<iRun<<endl;
      TObjArray *recal=(TObjArray*)cont->GetObject(RunIdVec.at(iRun));
      if(!recal)cout<<"Error - No bad map for run Number "<<RunIdVec.at(iRun)<<" online!!"<<endl;
      
      plot2D_Bad_OADB ->Reset();
      plot2D_Dead_OADB->Reset();
      for(Int_t iSM = 0; iSM < nSM; iSM++)
	{
	  //..Bad map fore each module
	  h[iSM]=(TH1C *)recal->FindObject(Form("EMCALBadChannelMap_Mod%d",iSM));
	  //..Loop though the SM to set which cells are bad
	  for(Int_t column=0;column<48;column++)
	    {
	      for(Int_t row=0;row<24;row++)
		{
		  Int_t inRow=row;
		  Int_t inCol=column;
		  Int_t trash = 0;
		  Int_t cellColumnAbs = 0;
		  Int_t  cellRowAbs = 0;
		  int cellID=geom->GetAbsCellIdFromCellIndexes(iSM,inRow,inCol);
		  fCaloUtils->GetModuleNumberCellIndexesAbsCaloMap(cellID,0,inCol,inRow,trash,cellColumnAbs,cellRowAbs);
		  //1 = dead, 2 = bad, 3 = warm
		  if(h[iSM]->GetBinContent(column,row) == 3)//..warm
		    {
		      plot2D_Bad_OADB->SetBinContent(cellColumnAbs,cellRowAbs,3);
		      warmChn.push_back(cellID);
		    }
		  if(h[iSM]->GetBinContent(column,row) == 2)//..bad
		    {
		      plot2D_Bad_OADB->SetBinContent(cellColumnAbs,cellRowAbs,2);
		      badChn.push_back(cellID);
		    }
		  if(h[iSM]->GetBinContent(column,row)==1)//..dead
		    {
		      plot2D_Dead_OADB->SetBinContent(cellColumnAbs,cellRowAbs,1);
		      deadChn.push_back(cellID);
		    }
		}//end for over rows
	    }//end for over cols
	}//end for over SM
    }//end for over runs

  plot2D_Bad_OADB->Draw("colz");
  bool val = true;
  for(int i = 0; i < warmChn.size(); i++)
    {
      int temp = 0;
      if(i == 0)
	cout << "warm_emcal=[";

      cout << warmChn.at(i);
      
      if(val && warmChn.at(i+1) > 12288)
	{
	  cout << "]" << endl;
	  cout << "warm_dcal=[";
	  val = false;
	  temp = i;
	}
      
      
      
      if( i != warmChn.size()-1 && !temp)
	cout << ", ";
      
      if(i == warmChn.size()-1)
	cout << "]" << endl;
    }

  val = true;
  for(int i = 0; i < badChn.size(); i++)
    {
      int temp = 0;
      if(i == 0)
	cout << "bad_emcal=[";

      cout << badChn.at(i);
      
      if(val && badChn.at(i+1) > 12288)
	{
	  cout << "]" << endl;
	  cout << "bad_dcal=[";
	  val = false;
	  temp = i;
	}      
      
      if( i != badChn.size()-1 && !temp)
	cout << ", ";
      
      if(i == badChn.size()-1)
	cout << "]" << endl;
    }//bad

  val = true;
  for(int i = 0; i < deadChn.size(); i++)
    {
      int temp = 0;
      if(i == 0)
	cout << "dead_emcal=[";

      cout << deadChn.at(i);
      
      if(val && deadChn.at(i+1) > 12288)
	{
	  cout << "]" << endl;
	  cout << "dead_dcal=[";
	  val = false;
	  temp = i;
	}
      
      if( i != deadChn.size()-1 && !temp)
	cout << ", ";
      
      if(i == deadChn.size()-1)
	cout << "]" << endl;
    }//dead*/
  
}//end accesingBadChannelMaps

