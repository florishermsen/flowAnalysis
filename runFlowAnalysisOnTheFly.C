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


#include "config.h"

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
   AliFlowEventSimpleMakerOnTheFly_mod *eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly_mod(uiSeed);
   eventMakerOnTheFly->SetCClass(cClass);
   eventMakerOnTheFly->SetMinMult(iMinMult);
   eventMakerOnTheFly->SetMaxMult(iMaxMult); 
   eventMakerOnTheFly->SetV1(dV1);
   eventMakerOnTheFly->SetV2(dV2);
   eventMakerOnTheFly->SetEtaRange(minEta,maxEta);
   eventMakerOnTheFly->SetPtRange(minPt,maxPt);
   eventMakerOnTheFly->SetUniformEfficiency(uniformEfficiency);
   eventMakerOnTheFly->Init();
   
   // Configure the flow analysis method:
   AliFlowAnalysisWithMCEventPlane_mod *mcep = new AliFlowAnalysisWithMCEventPlane_mod();
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
   TDirectoryFile *dirFileFinal = NULL;
   TString fileName="outputMCEPanalysis";
   dirFileFinal = new TDirectoryFile(fileName.Data(),fileName.Data());
 
   // i) Calculate and store the final results of all methods:
   mcep->Finish();
   mcep->WriteHistograms(dirFileFinal);
 
   outputFile->Close();

   delete outputFile;


   if (mcep) delete mcep;
   if (cutsRP) delete cutsRP;
   if (cutsPOI) delete cutsPOI;
   if (eventMakerOnTheFly) delete eventMakerOnTheFly;

 
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

