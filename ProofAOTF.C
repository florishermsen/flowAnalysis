

// Toggle random or same seed for random generator
Bool_t bSameSeed = kFALSE;

// Set range for rapidity
Double_t minEta = -.9;
Double_t maxEta = .9;
Double_t minEtaPlot = -.9;
Double_t maxEtaPlot = .9;
Int_t etaBins = 60;

// Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 60; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 61; // uniformly sampled multiplicity is < iMaxMult

// Parametrize the phi distribution, enter dvn if vn is eta-dependent:
Double_t dV1 = 0.0164; // constant harmonic v1
Double_t dV2 = 0.0; // constant harmonic v2 
Double_t dV3 = 0.0; // constant harmonic v3
Double_t dV4 = 0.0; // constant harmonic v4
Double_t dV5 = 0.0; // constant harmonic v5
Double_t dV6 = 0.0; // constant harmonic v6

// Parametrize the pt distribution:
//    Remark: Hardwired is Boltzmann distribution f(pt) = pt*exp[-sqrt(dMass^2+pt^2)/dT] 
Double_t dMass = 0.13957; // mass in GeV/c^2 (e.g. m_{pions} = 0.13957)
Double_t dTemperature = 0.44; // "temperature" in GeV/c (increase this parameter to get more high pt particles) 

// f) Determine how many times each sampled particle will be taken in the analysis (simulating nonflow):
Int_t nTimes = 1; // e.g. for nTimes = 2, strong 2-particle nonflow correlations are introduced 


// CURRENTLY IRRELEVANT PARAMETERS

// g1) Configure detector's acceptance:
Bool_t uniformAcceptance = kTRUE; // if kTRUE: detectors has uniform azimuthal acceptance.
                                  // if kFALSE: you will simulate detector with non-uniform acceptance in one or 
                                  // two sectors. For each sector you specify phiMin, phiMax and probability p. 
                                  // Then all particles emitted in direction phiMin < phi < phiMax will be taken 
                                  // with probability p. If p = 0, that sector is completely blocked. Set bellow 
                                  // phiMin1, phiMax1, p1 for the first sector and phiMin2, phiMax2, p2 for the second 
                                  // sector. If you set phiMin2 = phiMax2 = p2 = 0, only first non-uniform sector is 
                                  // simulated.
// 1st non-uniform sector:
Double_t phiMin1 = 60; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax1 = 120; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p1 = 0.5; // probablitity that particles emitted in [phiMin1,phiMax1] are taken
// 2nd non-uniform sector:
Double_t phiMin2 = 0.; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax2 = 0.; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p2 = 0.; // probablitity that particles emitted in [phiMin2,phiMax2] are taken

// g2) Configure detector's efficiency:
Bool_t uniformEfficiency = kTRUE; // if kTRUE: detectors has uniform pT efficiency
                                  // if kFALSE: you will simulate detector with non-uniform pT efficiency. 
                                  // Then all particles emitted in ptMin <= pt < ptMax will be taken 
                                  // with probability p, to be specified in lines just below. 
Double_t ptMin = 0.8; // non-uniform efficiency vs pT starts at pT = fPtMin
Double_t ptMax = 1.2; // non-uniform efficiency vs pT ends at pT = fPtMax
Double_t p = 0.5; // probablitity that particles emitted in [ptMin,ptMax> are taken

// i) Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = 0.0; // in GeV
Double_t ptMaxRP = 10.0; // in GeV
Double_t etaMinRP = minEta;
Double_t etaMaxRP = maxEta;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kTRUE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// j) Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = 0.0; // in GeV
Double_t ptMaxPOI = 10.0; // in GeV
Double_t etaMinPOI = -1.; // 
Double_t etaMaxPOI = 1.;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kTRUE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1

// k) Define the ranges for two subevents separated with eta gap (needed only for SP method):
Double_t etaMinA = -0.8; // minimum eta of subevent A
Double_t etaMaxA = -0.5; // maximum eta of subevent A
Double_t etaMinB = 0.5; // minimum eta of subevent B
Double_t etaMaxB = 0.8; // maximum eta of subevent B 


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
   eventMakerOnTheFly->SetMass(dMass);
   eventMakerOnTheFly->SetTemperature(dTemperature);
   eventMakerOnTheFly->SetV1(dV1);
   eventMakerOnTheFly->SetV2(dV2);
   eventMakerOnTheFly->SetV3(dV3);
   eventMakerOnTheFly->SetV3(dV4);
   eventMakerOnTheFly->SetV3(dV5);
   eventMakerOnTheFly->SetV3(dV6);
   eventMakerOnTheFly->SetEtaRange(minEta,maxEta);
   eventMakerOnTheFly->SetSubeventEtaRange(etaMinA,etaMaxA,etaMinB,etaMaxB); 
   eventMakerOnTheFly->SetNTimes(nTimes);

   if(!uniformAcceptance)
   {
      eventMakerOnTheFly->SetUniformAcceptance(kFALSE);
      eventMakerOnTheFly->SetFirstSectorPhiMin(phiMin1);
      eventMakerOnTheFly->SetFirstSectorPhiMax(phiMax1);
      eventMakerOnTheFly->SetFirstSectorProbability(p1);
      eventMakerOnTheFly->SetSecondSectorPhiMin(phiMin2);
      eventMakerOnTheFly->SetSecondSectorPhiMax(phiMax2);
      eventMakerOnTheFly->SetSecondSectorProbability(p2);
   } 
   if(!uniformEfficiency)
   {
      eventMakerOnTheFly->SetUniformEfficiency(kFALSE);
      eventMakerOnTheFly->SetPtMin(ptMin);
      eventMakerOnTheFly->SetPtMax(ptMax);
      eventMakerOnTheFly->SetPtProbability(p);
   }
   eventMakerOnTheFly->Init();

   mcep = new AliFlowAnalysisWithMCEventPlane_mod();
   mcep->SetEtaPlotRange(minEtaPlot, maxEtaPlot);
   mcep->SetEtaBins(etaBins);
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
   fOutput->Print();
}

