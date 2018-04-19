/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************ 
 * Create an event and perform full *
 * flow analysis 'on the fly'.      * 
 *                                  * 
 * author: Ante Bilandzic           * 
 *         (abilandzic@gmail.com)   *
 ************************************/ 

#ifndef ALIFLOWEVENTSIMPLEMAKERONTHEFLY_MOD_H
#define ALIFLOWEVENTSIMPLEMAKERONTHEFLY_MOD_H

class TF1;
class TRandom3;
class TH3F;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowTrackSimpleCuts;
        
class AliFlowEventSimpleMakerOnTheFly_mod{
   public:
      AliFlowEventSimpleMakerOnTheFly_mod(UInt_t uiSeed = 0); // constructor
      virtual ~AliFlowEventSimpleMakerOnTheFly_mod(); // destructor
      virtual void Init();   
      Bool_t AcceptPt(AliFlowTrackSimple *pTrack);  
      AliFlowEventSimple* CreateEventOnTheFly(AliFlowTrackSimpleCuts const *cutsRP, AliFlowTrackSimpleCuts const *cutsPOI); 
      // Setters and getters:
      void SetMinMult(Int_t iMinMult) {this->fMinMult = iMinMult;}
      Int_t GetMinMult() const {return this->fMinMult;} 
      void SetMaxMult(Int_t iMaxMult) {this->fMaxMult = iMaxMult;}
      Int_t GetMaxMult() const {return this->fMaxMult;}
      void SetV1(Double_t dV1) {this->fV1 = dV1;}
      Double_t GetV1() const {return this->fV1;} 
      void SetV2(Double_t dV2) {this->fV2 = dV2;}
      Double_t GetV2() const {return this->fV2;} 
      void SetEtaRange(Double_t minEta, Double_t maxEta) {this->fEtaMin = minEta;this->fEtaMax = maxEta;};
      void SetPtRange(Double_t minPt, Double_t maxPt) {this->fPtMin = minPt;this->fPtMax = maxPt;};
      void SetUniformEfficiency(Bool_t ue) {this->fUniformEfficiency = ue;}
      Bool_t GetUniformEfficiency() const {return this->fUniformEfficiency;} 

   private:
      AliFlowEventSimpleMakerOnTheFly_mod(const AliFlowEventSimpleMakerOnTheFly_mod& anAnalysis); // copy constructor
      AliFlowEventSimpleMakerOnTheFly_mod& operator=(const AliFlowEventSimpleMakerOnTheFly_mod& anAnalysis); // assignment operator
      Int_t fCount; // count number of events 
      Int_t fMinMult; // uniformly sampled multiplicity is >= iMinMult
      Int_t fMaxMult; // uniformly sampled multiplicity is < iMaxMult
      TF1 *fPtSpectra; // transverse momentum distribution (pt is sampled from hardwired Boltzmann distribution)
      TF1 *fPhiDistribution; // azimuthal distribution (phi is sampled from hardwired Fourier-like distribution)
      TF1 *fEtaDistribution; // rapidity distribution
      Double_t fV1; // harmonic v1
      Double_t fV2; // harmonic v2
      Double_t fEtaMin; // minimum eta
      Double_t fEtaMax; // maximum eta
      Double_t fPtMin; // minimum Pt
      Double_t fPtMax; // maximum Pt
      Double_t fPi; // pi
      Bool_t fUniformEfficiency; // detector has uniform efficiency vs pT, or perhaps not...

   ClassDef(AliFlowEventSimpleMakerOnTheFly_mod,1) // macro for rootcint
};
 
#endif
