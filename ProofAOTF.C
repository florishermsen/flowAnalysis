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

// Toggle random or same seed for random generator
Bool_t bSameSeed = kFALSE;

// Set transverse momentum profile
Double_t minPt = 0.;
Double_t maxPt = 50.;
Int_t ptBins = 100; //bins for result histograms

// Set rapidity profile
Double_t minEta = -.9;
Double_t maxEta = .9;
Int_t etaBins = 60;  //bins for result histograms

// Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 60; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 61; // uniformly sampled multiplicity is < iMaxMult

// Parametrize the phi distribution, enter dVn if vn is eta-dependent:
Double_t dV1 = 0.0164; // constant harmonic v1
Double_t dV2 = 0.0; // constant harmonic v2

// g2) Configure detector's efficiency:
Bool_t uniformEfficiency = kFALSE; // if kTRUE: detectors has uniform pT efficiency
                                  // if kFALSE: you will simulate detector with non-uniform pT efficiency. 


// i) Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = minPt; // in GeV
Double_t ptMaxRP = maxPt; // in GeV
Double_t etaMinRP = minEta;
Double_t etaMaxRP = maxEta;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kTRUE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// j) Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = minPt; // in GeV
Double_t ptMaxPOI = maxPt; // in GeV
Double_t etaMinPOI = minEta;
Double_t etaMaxPOI = maxEta;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kTRUE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1


R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH($HOME/alice/workspace/v1)

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
   fProfileRP = NULL;
   fProfilePOI = NULL;
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
   fProfileRP = dynamic_cast<TProfile*>(fSlaveHistList->FindObject("FlowPro_VetaRP_MCEP"));
   fProfilePOI = dynamic_cast<TProfile*>(fSlaveHistList->FindObject("FlowPro_VetaPOI_MCEP"));
   
   fOutput->Add(fProfileRP);
   fOutput->Add(fProfilePOI);
}

void ProofAOTF::Terminate()
{
   TString outputFileName = "results/ProofAnalysisResults_"+to_string(time(0))+".root";  
   TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");

   TDirectoryFile *dirFileFinal = NULL;
   TString fileName = "outputMCEPanalysis"; 
   dirFileFinal = new TDirectoryFile(fileName.Data(),fileName.Data());

   TList *fMasterHistList = new TList();
   fProfileRP = dynamic_cast<TProfile*>(fOutput->FindObject("FlowPro_VetaRP_MCEP"));
   fProfilePOI = dynamic_cast<TProfile*>(fOutput->FindObject("FlowPro_VetaPOI_MCEP"));
   fMasterHistList->Add(fProfileRP);
   fMasterHistList->Add(fProfilePOI);
   
   fMasterHistList->SetName("cobjMCEP");
   dirFileFinal->Add(fMasterHistList);
   dirFileFinal->Write(fMasterHistList->GetName(), TObject::kSingleKey);
   
   outputFile->Close();
   delete outputFile;
}





