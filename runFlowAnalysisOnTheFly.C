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
Int_t iNevts = 2000;

// Toggle random or same seed for random generator
Bool_t bSameSeed = kFALSE;

// Set range for rapidity
Double_t minEta = -1.;
Double_t maxEta = 1.;
Double_t minEtaPlot = -1.5;
Double_t maxEtaPlot = 1.5;

// Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 500; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 501; // uniformly sampled multiplicity is < iMaxMult

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

// h) Decide which flow analysis methods you will use:
Bool_t MCEP     = kTRUE; // Monte Carlo Event Plane

// i) Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = 0.0; // in GeV
Double_t ptMaxRP = 10.0; // in GeV
Double_t etaMinRP = minEta;
Double_t etaMaxRP = maxEta;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kFALSE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// j) Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = 0.0; // in GeV
Double_t ptMaxPOI = 10.0; // in GeV
Double_t etaMinPOI = -1.; // 
Double_t etaMaxPOI = 1.;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kFALSE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1

// k) Define the ranges for two subevents separated with eta gap (needed only for SP method):
Double_t etaMinA = -0.8; // minimum eta of subevent A
Double_t etaMaxA = -0.5; // maximum eta of subevent A
Double_t etaMinB = 0.5; // minimum eta of subevent B
Double_t etaMaxB = 0.8; // maximum eta of subevent B 

// l) Enable/disable usage of particle weights:
Bool_t usePhiWeights = kFALSE; // phi weights
Bool_t usePtWeights  = kFALSE; // pt weights 
Bool_t useEtaWeights = kFALSE; // eta weights
                                          
#include "TStopwatch.h"
#include "TObjArray.h"
#include "Riostream.h"
#include "TFile.h"


void CheckUserSettings()
{
 // Check if user settings make sense before taking off.

 if(iNevts <= 0)
 {
  printf("\n WARNING: nEvts <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 } 
 if(iMinMult < 0.)
 {
  printf("\n WARNING: iMinMult < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMaxMult <= 0.)
 {
  printf("\n WARNING: iMaxMult <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMinMult >= iMaxMult)
 {
  printf("\n WARNING: iMinMult >= iMaxMult !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dMass < 0.)
 {
  printf("\n WARNING: dMass < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dTemperature <= 1e-44)
 {
  printf("\n WARNING: dTemperature <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV1) > 0.5)
 {
  printf("\n WARNING: |dV1| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV2) > 0.5)
 {
  printf("\n WARNING: |dV2| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV3) > 0.5)
 {
  printf("\n WARNING: |dV3| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }

 if(!uniformAcceptance && phiMin1 > phiMax1)
 {
  cout<<" WARNING: You must have phiMin1 < phiMax1 !!!!"<<endl;
  exit(0);
 }
 if(!uniformAcceptance && !((TMath::Abs(phiMin2) < 1.e-44) && (TMath::Abs(phiMax2) < 1.e-44) && (TMath::Abs(p2) < 1.e-44)) 
    && (phiMin2 < phiMax1 || phiMin2 > phiMax2))
 {
  cout<<" WARNING: You must have phiMin2 > phiMax1 and phiMin2 < phiMax2 !!!!"<<endl;
  exit(0);
 }
 if((phiMin1 < 0 || phiMin1 > 360) || (phiMax1 < 0 || phiMax1 > 360) || 
    (phiMin2 < 0 || phiMin2 > 360) || (phiMax2 < 0 || phiMax2 > 360) )
 {
  cout<<" WARNING: You must take azimuthal angles from interval [0,360] !!!!"<<endl;
  exit(0);
 }
 if((p1 < 0 || p1 > 1) || (p2 < 0 || p2 > 1))
 {
  cout<<" WARNING: you must take p1 and p2 from interval [0,1] !!!!"<<endl;
  exit(0);
 }

} // end of void CheckUserSettings()

void WelcomeMessage()
{
 // Welcome.

 cout<<endl;
 cout<<endl;
 cout<<"      ---- ARE YOU READY TO FLY ? ----      "<<endl;
 cout<<endl;
 
 gSystem->Sleep(1544);
 
 cout<<endl;
 cout<<" ---- BEGIN FLOW ANALYSIS 'ON THE FLY' ---- "<<endl;
 cout<<endl;
 cout<<endl;

 gSystem->Sleep(1544);

} // end of void WelcomeMessage()

void SetupPar(char const *pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.
 
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_PHYSICS")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  } 
  if ( gSystem->AccessPathName(pararchivename) ) {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    
    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
} // end of void SetupPar()

void LoadLibraries() {
  

  R__ADD_INCLUDE_PATH($ALICE_ROOT)
  R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
  #include <AliFlowEventSimpleMakerOnTheFly_mod.cxx++>
  #include <AliFlowAnalysisWithMCEventPlane_mod.cxx++>

  cerr<<"Customized macros loaded...\n\n"<<endl;
  
} // end of void LoadLibraries()

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
 CheckUserSettings();
 WelcomeMessage();
 TStopwatch timer;
 timer.Start(); 
 LoadLibraries();
 
 // b) Initialize the flow event maker 'on the fly':
 UInt_t uiSeed = 0; // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
 if(bSameSeed){uiSeed = 44;} 
 AliFlowEventSimpleMakerOnTheFly_mod* eventMakerOnTheFly = new AliFlowEventSimpleMakerOnTheFly_mod(uiSeed);
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

 // c) If enabled, access particle weights from external file: 
 TFile *fileWithWeights = NULL;
 TList *listWithWeights = NULL; 
 if(usePhiWeights||usePtWeights||useEtaWeights) 
 {
  fileWithWeights = TFile::Open("weights.root","READ");
  if(fileWithWeights) 
  {
   listWithWeights = (TList*)fileWithWeights->Get("weights");
  }
  else
  {
   cout << " WARNING: the file <weights.root> with weights from the previous run was not found."<<endl;
   return 1;
  }    
 } // end of if(usePhiWeights||usePtWeights||useEtaWeights) 

 // d) Configure the flow analysis methods:
 AliFlowAnalysisWithQCumulants *qc = NULL;
 AliFlowAnalysisWithCumulants *gfc = NULL;
 AliFlowAnalysisWithFittingQDistribution *fqd = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz1sum  = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz1prod = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz2sum  = NULL;
 AliFlowAnalysisWithLeeYangZeros *lyz2prod = NULL;
 AliFlowAnalysisWithLYZEventPlane *lyzep = NULL;
 AliFlowAnalysisWithScalarProduct *sp = NULL;
 AliFlowAnalysisWithMixedHarmonics *mh = NULL;
 AliFlowAnalysisWithNestedLoops *nl = NULL;
 AliFlowAnalysisWithMCEventPlane_mod *mcep = NULL;   
 AliFlowAnalysisWithMultiparticleCorrelations *mpc = NULL;

 //fix otherwise missing decleration for C++11
 AliFlowLYZEventPlane *ep = NULL;


 // MCEP = monte carlo event plane
 if(MCEP) 
 {
  //AliFlowAnalysisWithMCEventPlane_mod *mcep = new AliFlowAnalysisWithMCEventPlane_mod();
  mcep = new AliFlowAnalysisWithMCEventPlane_mod();
  mcep->SetEtaPlotRange(minEtaPlot, maxEtaPlot);
  mcep->SetHarmonic(1); // default is v2
  mcep->Init();
 } // end of if(MCEP)

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
  if(MCEP){mcep->Make(event);}
  delete event;
 } // end of for(Int_t i=0;i<iNevts;i++)

 // h) Create the output file and directory structure for the final results of all methods: 
 TString outputFileName = "AnalysisResults.root";  
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
 const Int_t nMethods = 1;
 TString method[nMethods] = {"MCEP"};
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
 if(MCEP){mcep->Finish();mcep->WriteHistograms(dirFileFinal[0]);}
 
 outputFile->Close();
 delete outputFile;
 
 cout<<endl;
 cout<<endl;
 cout<<" ---- LANDED SUCCESSFULLY ---- "<<endl;
 cout<<endl; 
 
 timer.Stop();
 cout << endl;
 timer.Print();

 return 0;

} // end of int runFlowAnalysisOnTheFly()
