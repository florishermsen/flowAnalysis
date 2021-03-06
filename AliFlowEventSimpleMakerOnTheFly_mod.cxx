/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

/************************************ 
 * Create an event and perform full *
 * flow analysis 'on the fly'.      * 
 *                                  * 
 * author: Ante Bilandzic           * 
 *         (abilandzic@gmail.com)   *
 ************************************/ 
  
#include "Riostream.h"
#include "TMath.h"
#include "TF1.h"
#include "TH3.h"
#include "TRandom3.h"
#include "AliFlowEventSimpleMakerOnTheFly_mod.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackSimpleCuts.h"

using std::endl;
using std::cout;
ClassImp(AliFlowEventSimpleMakerOnTheFly_mod)

//========================================================================================================================================

AliFlowEventSimpleMakerOnTheFly_mod::AliFlowEventSimpleMakerOnTheFly_mod(UInt_t uiSeed):
   fCount(0),
   fCClass(0),
   fMinMult(0),
   fMaxMult(0),  
   fPtSpectra(NULL),
   fPhiDistribution(NULL),
   fEtaDistribution(NULL),
   fV1(0.),
   fV2(0.05),
   fEtaMin(-1.),
   fEtaMax(1.),
   fPtMin(0),
   fPtMax(10.),
   fPi(TMath::Pi()),
   fUniformEfficiency(kTRUE)
{
   // Constructor.
  
   // Determine seed for gRandom:
   delete gRandom;
   gRandom = new TRandom3(uiSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

} // end of AliFlowEventSimpleMakerOnTheFly_mod::AliFlowEventSimpleMakerOnTheFly_mod(UInt_t uiSeed):

//====================================================================================================================

AliFlowEventSimpleMakerOnTheFly_mod::~AliFlowEventSimpleMakerOnTheFly_mod()
{
   // Destructor.

   if(fPtSpectra){delete fPtSpectra;}
   if(fPhiDistribution){delete fPhiDistribution;}
   if(fEtaDistribution){delete fEtaDistribution;}
   if(gRandom){delete gRandom;}

} // end of AliFlowEventSimpleMakerOnTheFly_mod::~AliFlowEventSimpleMakerOnTheFly_mod() 

//====================================================================================================================

void AliFlowEventSimpleMakerOnTheFly_mod::Init()
{
   // Book all objects in this method.

   // a) Define the pt spectra;
   // b) Define the phi distribution.
   // c) Define the eta distribution.

   // a) Define the pt spectra:

   // Centrality Class - Dependent RAA
   if(fCClass==0) {
      fPtSpectra = new TF1("fPtSpectra","[0]*x*sqrt(x^2+[1]^2)*(1+([2]-1)*sqrt(x^2+[1]^2)/[3])^(-[2]/([2]-1))*(0.186179+0.730187*exp(-0.22222*(-log(250)+log(x))^2)+0.172795*exp(-0.27845*(-1.22378+log(x))^2)*TMath::Erfc(2.53292*(-1.22378+log(x))))",fPtMin,fPtMax); // realistic d-meson spectrum, WITH R_AA
   } else if(fCClass==1) {
      fPtSpectra = new TF1("fPtSpectra","[0]*x*sqrt(x^2+[1]^2)*(1+([2]-1)*sqrt(x^2+[1]^2)/[3])^(-[2]/([2]-1))*(0.264522+0.705465*exp(-0.154321*(-log(250)+log(x))^2)+0.172813*exp(-0.305176*(-1.19392+log(x))^2)*TMath::Erfc(2.20971*(-1.19392+log(x))))",fPtMin,fPtMax); // realistic d-meson spectrum, WITH R_AA
   } else {
      fPtSpectra = new TF1("fPtSpectra","[0]*x*sqrt(x^2+[1]^2)*(1+([2]-1)*sqrt(x^2+[1]^2)/[3])^(-[2]/([2]-1))*(0.596343+0.299838*exp(-0.154321*(-log(250)+log(x))^2)+0.0560588*exp(-0.78125*(-1.02962+log(x))^2)*TMath::Erfc(2.65165*(-1.02962+log(x))))",fPtMin,fPtMax); // realistic d-meson spectrum, WITH R_AA
   }
   
   fPtSpectra->SetNpx(400);
   fPtSpectra->SetParameter(0, 8.4e8); // Scaling Constant
   fPtSpectra->SetParameter(1, 1.86962); // D-Mass (GeV)
   fPtSpectra->SetParameter(2, 1.15); // q-factor
   fPtSpectra->SetParameter(3, 7.5e-2); // q-Temp


   fPtSpectra->SetParNames("Scaling Constant","D mass","Q-Factor","Q-Temperature");

   // b) Define the phi distribution:
   Double_t dPhiMin = 0.; 
   Double_t dPhiMax = TMath::TwoPi();
   fPhiDistribution = new TF1("fPhiDistribution","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))",dPhiMin,dPhiMax);
   fPhiDistribution->SetNpx(100);

   fPhiDistribution->SetParName(0,"Reaction Plane");
   fPhiDistribution->SetParameter(0,0.);
   fPhiDistribution->SetParName(1,"Directed Flow (v1)"); 
   fPhiDistribution->SetParameter(1,fV1);
   fPhiDistribution->SetParName(2,"Elliptic Flow (v2)");
   fPhiDistribution->SetParameter(2,fV2);

   // c) Define the eta distribution:
   fEtaDistribution = new TF1("fEtaDistribution","1+[0]*x^2",fEtaMin,fEtaMax); // % dip around 0

   fEtaDistribution->SetParameter(0,0.056); //10-30
   if(fCClass==1) {
      fEtaDistribution->SetParameter(0,0.061); //30-50
   } else if(fCClass==2) {
      fEtaDistribution->SetParameter(0,0.074); //60-80
   }

} // end of void AliFlowEventSimpleMakerOnTheFly_mod::Init()


//====================================================================================================================

Bool_t AliFlowEventSimpleMakerOnTheFly_mod::AcceptPt(AliFlowTrackSimple *pTrack)
{
   // For the case of non-uniform efficiency determine in this method if particle is accepted or rejected for a given pT.
   
   Double_t efficiencyBins[3][13][2] = {{ //10-30
      {2, 0},
      {3, 0.000767664},
      {4, 0.00463984},
      {5, 0.0153918},
      {6, 0.0321692},
      {7, 0.0528523},
      {8, 0.0582403},
      {10, 0.123859},
      {12, 0.145762},
      {16, 0.16983},
      {24, 0.221823},
      {36, 0.574282},
      {50, 0.553414}
   },{ //30-50
      {2, 0},
      {3, 0.00219054},
      {4, 0.012685},
      {5, 0.0190254},
      {6, 0.0351793},
      {7, 0.0726777},
      {8, 0.0877877},
      {10, 0.142911},
      {12, 0.167198},
      {16, 0.199494},
      {24, 0.250309},
      {36, 0.664265},
      {50, 1}
   },{ //60-80
      {1, 0},
      {2, 0.00134291},
      {3, 0.0119252},
      {4, 0.0249153},
      {5, 0.0467732},
      {6, 0.0523987},
      {7, 0.110127},
      {8, 0.152198},
      {10, 0.188413},
      {12, 0.302068},
      {16, 0.558807},
      {24, 0.653803},
      {50, 1},
      }
   };

   Bool_t bAccept = kTRUE;

   for ( Int_t i = 0; i < 13; i++ ) {
      if(pTrack->Pt() < efficiencyBins[fCClass][i][0]) {
         if(gRandom->Uniform(0,1) > efficiencyBins[fCClass][i][1]) {
            bAccept = kFALSE; // no mercy!
         }
         break; //very important break statement, otherwise lose ALL low energy tracks
      }
   }

   return bAccept;
 
} // end of Bool_t AliFlowEventSimpleMakerOnTheFly_mod::AcceptPt(AliFlowTrackSimple *pTrack);

//====================================================================================================================

AliFlowEventSimple* AliFlowEventSimpleMakerOnTheFly_mod::CreateEventOnTheFly(AliFlowTrackSimpleCuts const *cutsRP, AliFlowTrackSimpleCuts const *cutsPOI)
{
   // Method to create event 'on the fly'.

   // a) Determine the multiplicity of an event;
   // b) Determine the reaction plane of an event;
   // c) If v2 fluctuates uniformly event-by-event, sample its value from [fMinV2,fMaxV2];
   // d) Create event 'on the fly';
   // e) Cosmetics for the printout on the screen.

   // a) Determine the multiplicity of an event:
   //Int_t iMult = (Int_t)gRandom->Uniform(fMinMult,fMaxMult);
   Int_t iMult = fMinMult;


   // b) Determine the reaction plane of an event:
   Double_t dReactionPlane = gRandom->Uniform(0.,TMath::TwoPi());
   fPhiDistribution->SetParameter(0,dReactionPlane);

   // d) Create event 'on the fly':
   AliFlowEventSimple *pEvent = new AliFlowEventSimple(iMult);
   pEvent->SetReferenceMultiplicity(iMult);
   pEvent->SetMCReactionPlaneAngle(dReactionPlane);

   Int_t nRPs = 0; // number of particles tagged RP in this event
   Int_t nPOIs = 0; // number of particles tagged POI in this event

   for(Int_t p=0;p<iMult;p++)
   {
      AliFlowTrackSimple *pTrack = new AliFlowTrackSimple();

      pTrack->SetPt(fPtSpectra->GetRandom()); 

      // Check pT efficiency:
      if(!fUniformEfficiency && !this->AcceptPt(pTrack)) {
         delete pTrack; // very important, otherwise mem leak
         continue;
      }


      // Eta-dependent and charge-dependent v1:

      pTrack->SetEta(fEtaDistribution->GetRandom());
      pTrack->SetCharge((gRandom->Integer(2)>0.5 ? 1 : -1));

      //Double_t currentV1 = fV1*(1-1/(0.5+pTrack->Pt())); // legacy code from pt-dependent v1
      fPhiDistribution->SetParameter(1,pTrack->Eta()*pTrack->Charge()*fV1);
      pTrack->SetPhi(fPhiDistribution->GetRandom());

      // Checking the RP cuts:     
      if(cutsRP->PassesCuts(pTrack))
      {
         pTrack->TagRP(kTRUE); 
         nRPs++; 
      }
      // Checking the POI cuts:    
      if(cutsPOI->PassesCuts(pTrack))
      {
         pTrack->TagPOI(kTRUE); 
         nPOIs++;
      }
      
      pEvent->AddTrack(pTrack);
   } // end of for(Int_t p=0;p<iMult;p++)
   pEvent->SetNumberOfRPs(nRPs);
   pEvent->SetNumberOfPOIs(nPOIs);

   // introducing limited angular resolution
   // set error on event plane angle after-the-fact for use in reconstruction
   Double_t dReactionPlaneWithError;
   if(fCClass==2) {
      Double_t dReactionPlaneWithError = gRandom->Gaus(dReactionPlane, 0.942);
   } else {
      Double_t dReactionPlaneWithError = gRandom->Gaus(dReactionPlane, 0.628);
   }
   pEvent->SetMCReactionPlaneAngle(dReactionPlaneWithError);


   // e) Cosmetics for the printout on the screen:
   Int_t cycle = 100;
   if((++fCount % cycle) == 0) 
   {
      if(TMath::Abs(dReactionPlane)>1.e-44) 
      {
         cout<<" MC Reaction Plane Angle = "<<dReactionPlane<<endl;
      } else 
      {
         cout<<" MC Reaction Plane Angle is unknown :'( "<< endl;
      }     
      cout<<" # of simulated tracks  = "<<iMult<<endl;
      cout<<" # of RP tagged tracks  = "<<nRPs<<endl;
      cout<<" # of POI tagged tracks = "<<nPOIs<<endl;  
      cout <<"  .... "<<fCount<< " events processed ...."<<endl;
   } // end of if((++fCount % cycle) == 0) 

   return pEvent;
    
} // end of CreateEventOnTheFly()
 
//====================================================================================================================

