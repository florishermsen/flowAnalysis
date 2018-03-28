/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//////////                                         //////////
//////////       runFlowAnalysisOnTheFly.c         //////////
//////////                                         //////////
//////////    Original macro by Ante Bilandzic     //////////
//////////                                         //////////
//////////   Adapted to ROOT6 by F.A.W. Hermsen    //////////
//////////                                         //////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


// Number of events
Int_t iNevts = 100000;

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


#include <ctime>
#include <string>

#include "TStopwatch.h"
#include "TObjArray.h"
#include "Riostream.h"
#include "TFile.h"

#include "AliFlowEventSimpleMakerOnTheFly_mod.h"
#include "AliFlowAnalysisWithMCEventPlane_mod.h"
#include "AliFlowEventSimpleMakerOnTheFly_mod.cxx"
#include "AliFlowAnalysisWithMCEventPlane_mod.cxx"

void WelcomeMessage()
{
   // Welcome.
   cout<<endl;
   cout<<endl;
   cout<<"      ---- ARE YOU READY TO FLY ? ----      "<<endl;
   cout<<endl;
 
   gSystem->Sleep(500);

   cout<<endl;
   cout<<" ---- BEGIN FLOW ANALYSIS 'ON THE FLY' ---- "<<endl;
   cout<<endl;
   cout<<endl;

   gSystem->Sleep(500);

} // end of void WelcomeMessage()

int runFlowAnalysisOnTheFly()
{
   
   // Beging analysis 'on the fly'.

   // a) Formal necessities....;
   // b) Initialize the flow event maker 'on the fly';
   // c) If enabled, access particle weights from external file; 
   // d) Configure the flow analysis methods;
   // e) Simple cuts for RPs;
   // f) Simple cuts for POIs;
   // g) Create and analyse events 'on the fly'; 
   // h) Create the output file and directory structure for the final results of all methods; 
   // i) Calculate and store the final results of all methods.

   // a) Formal necessities....:
   WelcomeMessage();
   TStopwatch timer;
   timer.Start(); 
 
   // b) Initialize the flow event maker 'on the fly':
   UInt_t uiSeed = 0; // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
   if(bSameSeed){uiSeed = 44;} 
   AliFlowEventSimpleMakerOnTheFly_mod* eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly_mod(uiSeed);
   eventMakerOnTheFly->SetMinMult(iMinMult);
   eventMakerOnTheFly->SetMaxMult(iMaxMult); 
   eventMakerOnTheFly->SetV1(dV1);
   eventMakerOnTheFly->SetV2(dV2);
   eventMakerOnTheFly->SetEtaRange(minEta,maxEta);
   eventMakerOnTheFly->SetPtRange(minPt,maxPt);
   eventMakerOnTheFly->SetUniformEfficiency(uniformEfficiency);
   eventMakerOnTheFly->Init();
   
   // Configure the flow analysis method:
   AliFlowAnalysisWithMCEventPlane_mod *mcep = NULL;
   mcep = new AliFlowAnalysisWithMCEventPlane_mod();
   mcep->SetPtRange(minPt, maxPt);
   mcep->SetNbinsPt(ptBins);
   mcep->SetEtaRange(minEta, maxEta);
   mcep->SetNbinsEta(etaBins);
   mcep->SetHarmonic(1);
   mcep->Init();

   // e) Simple cuts for RPs: 
   AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
   cutsRP->SetPtMax(ptMaxRP);
   cutsRP->SetPtMin(ptMinRP);
   cutsRP->SetEtaMax(etaMaxRP);
   cutsRP->SetEtaMin(etaMinRP);
   cutsRP->SetPhiMax(phiMaxRP*TMath::Pi()/180.);
   cutsRP->SetPhiMin(phiMinRP*TMath::Pi()/180.);
   if(bUseChargeRP){cutsRP->SetCharge(chargeRP);}

   // f) Simple cuts for POIs: 
   AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
   cutsPOI->SetPtMax(ptMaxPOI);
   cutsPOI->SetPtMin(ptMinPOI);
   cutsPOI->SetEtaMax(etaMaxPOI);
   cutsPOI->SetEtaMin(etaMinPOI);
   cutsPOI->SetPhiMax(phiMaxPOI*TMath::Pi()/180.);
   cutsPOI->SetPhiMin(phiMinPOI*TMath::Pi()/180.);
   if(bUseChargePOI){cutsPOI->SetCharge(chargePOI);}
                                       
   // g) Create and analyse events 'on the fly':
   for(Int_t i=0;i<iNevts;i++) 
   {   
      // Creating the event 'on the fly':
      AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(cutsRP,cutsPOI); 
      // Passing the created event to flow analysis methods:
      mcep->Make(event);
      delete event;
   } // end of for(Int_t i=0;i<iNevts;i++)

   // h) Create the output file and directory structure for the final results of all methods: 
   TString outputFileName = "results/AnalysisResults_"+to_string(time(0))+".root";  
   TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
   const Int_t nMethods = 1;
   TString method[] = {"MCEP"};
   TDirectoryFile *dirFileFinal[nMethods] = {NULL};
   TString fileName[nMethods]; 
   for(Int_t i=0;i<nMethods;i++)
   {
      fileName[i]+="output";
      fileName[i]+=method[i].Data();
      fileName[i]+="analysis";
      dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
   }
 
   // i) Calculate and store the final results of all methods:
   mcep->Finish();
   mcep->WriteHistograms(dirFileFinal[0]);
 
   outputFile->Close();
   delete outputFile;
 
   cout<<endl;
   cout<<endl;
   cout<<" ---- LANDED SUCCESSFULLY ---- "<<endl;
   cout<<endl; 
 
   timer.Stop();
   cout << endl;
   timer.Print();
   cout << endl;
   return 0;

} // end of int runFlowAnalysisOnTheFly()

