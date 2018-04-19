/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//////////                                         //////////
//////////        Flow Analysis for PROOF          //////////
//////////                                         //////////
//////////    Original macro by Ante Bilandzic     //////////
//////////                                         //////////
//////////   Adapted to ROOT6 by F.A.W. Hermsen    //////////
//////////                                         //////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH($HOME/alice/workspace/v1)

#include "config.h"

#include <ctime>
#include <string>

//current file headers
#include "ProofAOTF.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom3.h"

#include "TSelector.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TObjArray.h"
#include "TFile.h"
#include "Riostream.h"

//macro loading headers
#include "TMacro.h"
#include "TROOT.h"
#include "TSystem.h"

// macro specific
#include "AliFlowEventSimpleMakerOnTheFly_mod.h"
#include "AliFlowAnalysisWithMCEventPlane_mod.h"
#include <AliFlowEventSimpleMakerOnTheFly_mod.cxx>
#include <AliFlowAnalysisWithMCEventPlane_mod.cxx>


//_____________________________________________________________________________
ProofAOTF::ProofAOTF()
{
   eventMakerOnTheFly = NULL;
   mcep = NULL;
   cutsRP = NULL;
   cutsPOI = NULL;
}

//_____________________________________________________________________________
ProofAOTF::~ProofAOTF()
{
   if (mcep) delete mcep;
   if (cutsRP) delete cutsRP;
   if (cutsPOI) delete cutsPOI;
   if (eventMakerOnTheFly) delete eventMakerOnTheFly;
}

void ProofAOTF::Begin(TTree * ) { }

void ProofAOTF::SlaveBegin(TTree * )
{
   UInt_t uiSeed = 0; // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
   if(bSameSeed){uiSeed = 44;}

   eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly_mod(uiSeed);
   eventMakerOnTheFly->SetMinMult(iMinMult);
   eventMakerOnTheFly->SetMaxMult(iMaxMult); 
   eventMakerOnTheFly->SetV1(dV1);
   eventMakerOnTheFly->SetV2(dV2);
   eventMakerOnTheFly->SetEtaRange(minEta,maxEta);
   eventMakerOnTheFly->SetPtRange(minPt,maxPt);
   eventMakerOnTheFly->SetUniformEfficiency(uniformEfficiency);
   eventMakerOnTheFly->Init();

   mcep = new AliFlowAnalysisWithMCEventPlane_mod();
   mcep->SetPtRange(minPt, maxPt);
   mcep->SetNbinsPt(ptBins);
   mcep->SetEtaRange(minEta, maxEta);
   mcep->SetNbinsEta(etaBins);
   mcep->SetHarmonic(1);
   mcep->Init();

   // e) Simple cuts for RPs: 
   cutsRP = new AliFlowTrackSimpleCuts();
   cutsRP->SetPtMax(ptMaxRP);
   cutsRP->SetPtMin(ptMinRP);
   cutsRP->SetEtaMax(etaMaxRP);
   cutsRP->SetEtaMin(etaMinRP);
   cutsRP->SetPhiMax(phiMaxRP*TMath::Pi()/180.);
   cutsRP->SetPhiMin(phiMinRP*TMath::Pi()/180.);
   if(bUseChargeRP){cutsRP->SetCharge(chargeRP);}

   // f) Simple cuts for POIs: 
   cutsPOI = new AliFlowTrackSimpleCuts();
   cutsPOI->SetPtMax(ptMaxPOI);
   cutsPOI->SetPtMin(ptMinPOI);
   cutsPOI->SetEtaMax(etaMaxPOI);
   cutsPOI->SetEtaMin(etaMinPOI);
   cutsPOI->SetPhiMax(phiMaxPOI*TMath::Pi()/180.);
   cutsPOI->SetPhiMin(phiMinPOI*TMath::Pi()/180.);
   if(bUseChargePOI){cutsPOI->SetCharge(chargePOI);}

}

Bool_t ProofAOTF::Process(Long64_t)
{

   AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(cutsRP,cutsPOI);
   mcep->Make(event);
   delete event;

   return kTRUE;
}

void ProofAOTF::SlaveTerminate()
{

   mcep->Finish();

   TList *fSlaveHistList = mcep->GetHistList();
   fSlaveHistList->SetName("cobjMCEP");
   fOutput->Add(fSlaveHistList->Clone());
}

void ProofAOTF::Terminate()
{
   TString outputFileName = "results/ProofAnalysisResults_"+to_string(time(0))+".root";  
   TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");

   TDirectoryFile *dirFileFinal = NULL;
   TString fileName = "outputMCEPanalysis"; 
   dirFileFinal = new TDirectoryFile(fileName.Data(),fileName.Data());
   dirFileFinal->Add(fOutput->FindObject("cobjMCEP")->Clone());
   dirFileFinal->Write(fOutput->GetName(), TObject::kSingleKey);
   
   outputFile->Close();
   delete outputFile;
}





