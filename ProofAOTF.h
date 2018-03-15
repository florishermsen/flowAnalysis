
#ifndef PROOFAOTF_H
#define PROOFAOTF_H

class TH1F;
class TRandom;
class TSelector;
class AliFlowEventSimpleMakerOnTheFly_mod;
class AliFlowAnalysisWithMCEventPlane_mod;
class AliFlowTrackSimpleCuts;
class ProofAOTF : public TSelector {
public :

   TH1F *fH1F;
   TRandom *fRandom;
   AliFlowEventSimpleMakerOnTheFly_mod *eventMakerOnTheFly;
   AliFlowAnalysisWithMCEventPlane_mod *mcep;
   AliFlowTrackSimpleCuts *cutsRP;
   AliFlowTrackSimpleCuts *cutsPOI;

   ProofAOTF();
   virtual ~ProofAOTF();
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(ProofAOTF,2);
};
#endif