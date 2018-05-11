
// Define centrality class
// 0  10-30%
// 1  30-50%
// 2  50-80%
Int_t cClass = 0;

// Number of events, only for non-PROOF
Int_t iNevts = 10;



// 10-30%
// Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 2200; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 2201; // uniformly sampled multiplicity is < iMaxMult

// Parametrize the phi distribution, enter dVn if vn is eta-dependent:
Double_t dV1 = 0.0255; // constant harmonic v1
Double_t dV2 = 0.0; // constant harmonic v2



/*
// 30-50%
	Int_t iMinMult = 700;
	Int_t iMaxMult = 701;
	Double_t dV1 = 0.0411;
	Double_t dV2 = 0.0;

// 50-80%
	Int_t iMinMult = 70;
	Int_t iMaxMult = 71;
	Double_t dV1 = 0.0525;
	Double_t dV2 = 0.0;
*/

// Toggle random or same seed for random generator
Bool_t bSameSeed = kFALSE;

// Set transverse momentum profile
Double_t minPt = 0.;
Double_t maxPt = 50.;
Int_t ptBins = 50; //bins for result histograms

// Set rapidity profile
Double_t minEta = -.9;
Double_t maxEta = .9;
Int_t etaBins = 60;  //bins for result histograms

// Configure detector's efficiency:
Bool_t uniformEfficiency = kFALSE; // if kTRUE: detectors has uniform pT efficiency
                                  // if kFALSE: simulate detector with non-uniform pT efficiency & acceptance. 


// Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = minPt; // in GeV
Double_t ptMaxRP = maxPt; // in GeV
Double_t etaMinRP = minEta;
Double_t etaMaxRP = maxEta;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kTRUE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = minPt; // in GeV
Double_t ptMaxPOI = maxPt; // in GeV
Double_t etaMinPOI = minEta;
Double_t etaMaxPOI = maxEta;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kTRUE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1


// Configure Pt cuts for extra pt-region v1 hists (not yet in macro)
Bool_t ptSubHists = kTRUE;
Double_t ptCutOffs[2] = {3,5};
