#ifndef ALIHFEHPPBTOOL_H
#define ALIHFEHPPBTOOL_H

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#include "TString.h"
#include "TArrayD.h"
#include "TH1F.h"

class TH2F;
class TFile;

class AliHFehpPbTool : public TNamed {
    
public:
    
    // Constructors and Destructors
    AliHFehpPbTool(const Char_t *name, const Char_t *title);
    AliHFehpPbTool();
    ~AliHFehpPbTool();
    
    void SetInputFileName(TString Name) {fInputFileName = Name;}
    void SetConfigurationName(TString Name) {fConfigurationName = Name;}
    void SetpTBins(Int_t n, Double_t* array) { fpTBins.Set(n,array); };
    void SetpTBinsResults(Int_t n, Double_t* array) { fpTBinsResults.Set(n,array); };
    
    
    Bool_t ReadpTHistograms();
    Bool_t ProcesspTHistograms();
    Bool_t MergeCorrelationDistributions();
    
    
    Bool_t ConnectToInputFile(TString FileName, TString ConfigurationName);
    Bool_t CalculateHadronContamination();
    Bool_t CalculateTaggingEfficiency();
    Bool_t CompareMixedEventDistributions();
    Bool_t CalculateHFeEfficiency();
    Bool_t MergeMCCorrelationDistributions();
    Bool_t MakeMCTruthMixed();
    Bool_t ReadMCDistributions(Int_t RebinX = 2, Int_t RebinY = 2);
    Bool_t CorrectMCDistributionsByMixing();
    Bool_t CorrelationCT();
    Bool_t NormalizeHFeCorrrelationMC();
    Bool_t NormalizeHFeCorrrelation();
    Bool_t ProjectMCTo1D();
    Bool_t ProjectTo1D();
    Bool_t CorrelationCT1D();
    Bool_t ReadTaggingEfficiencyFromFile();
    Bool_t CalculateTaggingEfficiencyW();
    
    void Process();
    
    void CalculateFlow(AliHFehpPbTool* Central, AliHFehpPbTool* Peri);
    
    Bool_t ReadAndProcessCorrelationDistributions(Int_t RebinX = 2, Int_t RebinY = 2);
    TString TrainConfiguration(
                               Int_t pTBin = 0,
                               Bool_t Correlation = kTRUE,
                               Bool_t ispp = kFALSE,
                               Bool_t isMC = kFALSE,
                               Double_t ElectronDCAxy = 0.25,
                               Double_t ElectronDCAz = 1.0,
                               Double_t HadronDCAxy = 0.25,
                               Double_t HadronDCAz = 1.0,
                               Double_t TPCPIDLow = -0.5,
                               Double_t TPCPIDUp = 3.0,
                               Double_t InvariantMassCut = 0.14,
                               Double_t pTCutPartner = 0.0,
                               Double_t MultiplicityLow = 0.,
                               Double_t MultiplicityUp = 100.,
                               Double_t HadronPtCutLow = 0.3,
                               Double_t HadronPtCutUp = 2.0,
                               Double_t EtaCutLow = -0.8,
                               Double_t EtaCutUp = 0.8,
                               Double_t NonHFEangleCut = 999,
                               Int_t NHitsITS = 2,
                               Int_t SPDLayers = 0,
                               Int_t TPCNCluster = 100,
                               Int_t TPCNClusterPartner = 60,
                               Int_t TPCNClusterPID = 80);
    
    static TString TaskName(
                            Int_t pTBin = 0,
                            Bool_t Correlation = kTRUE,
                            Bool_t ispp = kFALSE,
                            Bool_t isMC = kFALSE,
                            Double_t ElectronDCAxy = 0.25,
                            Double_t ElectronDCAz = 1.0,
                            Double_t HadronDCAxy = 0.25,
                            Double_t HadronDCAz = 1.0,
                            Double_t TPCPIDLow = -0.5,
                            Double_t TPCPIDUp = 3.0,
                            Double_t InvariantMassCut = 0.14,
                            Double_t pTCutPartner = 0.0,
                            Double_t MultiplicityLow = 0.,
                            Double_t MultiplicityUp = 100.,
                            Double_t HadronPtCutLow = 0.3,
                            Double_t HadronPtCutUp = 2.0,
                            Double_t EtaCutLow = -0.8,
                            Double_t EtaCutUp = 0.8,
                            Double_t NonHFEangleCut = 999,
                            Int_t NHitsITS = 2,
                            Int_t SPDLayers = 0,
                            Int_t TPCNCluster = 100,
                            Int_t TPCNClusterPartner = 60,
                            Int_t TPCNClusterPID = 80)
    {
        TString Name = Form("%d_%d_%d_%d_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_%d_%d_%d",pTBin, Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID);
        
        return Name;
        
    };
    
    TH1F* GetHFepT() {return fHFepT;};
    TH1F* GetEff() {return fEffTagging;};
    TH1F* GetHFeEff() {return fEffHFe;};
    TH1F* GetPurity() {return fHadronContamination;};
    TH1F* GetHFeh1D(Int_t pT ) {return fHFEhNormalized1D[pT];};
    TH2F* GetHFeh2D(Int_t pT ) {return fHFEhNormalized[pT];};
    Double_t GetNEvents();
    
    void SetUseEffElectrons(){fEffCorrectionForElectrons = kTRUE;};
    
    void SetTagEff(TH1F *eff) {fEffTagging =  (TH1F*) eff->Clone("fEffTagging"); };
    void SetHFeEff(TH1F *HFeEff) {fEffHFe = (TH1F*) HFeEff->Clone("fEffHFe"); };
    
    /*
     TString GetNameFromConfiguration(Bool_t Correlation = kTRUE,
     Bool_t ispp = kFALSE,
     Bool_t isMC = kTRUE,
     Bool_t FinepTBinning = kFALSE,
     Double_t ElectronDCAxy = 0.25,
     Double_t ElectronDCAz = 1.0,
     Double_t HadronDCAxy = 0.25,
     Double_t HadronDCAz = 1.0,
     Double_t TPCPIDLow = -0.5,
     Double_t TPCPIDUp = 3.0,
     Double_t InvariantMassCut = 0.14,
     Double_t pTCutPartner = 0.0,
     Double_t MultiplicityLow = 0.,
     Double_t MultiplicityUp = 100.,
     Double_t HadronPtCutLow = 0.3,
     Double_t HadronPtCutUp = 2.0,
     Double_t EtaCutLow = -0.8,
     Double_t EtaCutUp = 0.8,
     Double_t NonHFEangleCut = 999,
     Int_t NHitsITS = 4,
     Int_t SPDLayers = 0,
     Int_t TPCNCluster = 100,
     Int_t TPCNClusterPartner = 60,
     Int_t TPCNClusterPID = 80);
     */
    
    
    
private:
    
    //void FormatTH1(TH1F*, Int_t Color, Int_t Marker, Double_t MarkerSize);
    // Data
    
    Bool_t fPrintHistograms;
    Bool_t fEffCorrectionForElectrons;
    
    TFile *fInputFile; //!
    TList *fInputList; //!
    TArrayD fpTBins;
    TArrayD fpTBinsResults;
    TString fInputFileName; //!
    TString fConfigurationName;
    
    //Spectrum
    TH1F* fInclusivepT;
    TH1F* fULSpT;
    TH1F* fLSpT;
    TH1F* fBackgroundpT;
    TH1F* fHFepT;
    
    //Corrections
    TH1F* fEffTagging;
    TH1F* fHadronContamination;
    TH1F* fEffHFe;
    
    Double_t fTOFNSigma;
    
    //Inc Same
    TH2F **fSameIncEh;
    //Inc Mixed
    TH2F **fMixedIncEh;
    //Dihadron (only same)
    TH2F **fDiHadron;
    //Background Same
    TH2F **fSameBackULSEh;
    TH2F **fSameBackLSEh;
    //Background Mixed
    TH2F **fMixedBackULSEh;
    TH2F **fMixedBackLSEh;
    
    TH2F **fSameBackNoPULSEh;
    TH2F **fSameBackNoPLSEh;
    //Total Background
    TH2F **fBackEh;
    //HFe
    TH2F **fHFEh;
    TH2F **fHFEhNormalized;
    TH2F **fBackNonIDEh;
    TH2F **fHFEhSame;
    TH2F **fHFEhMixed;
    TH2F **fBackMixedEh;
    
    TH2F **fHFEMC;
    TH2F **fHFEhMCNormalized;
    TH2F **fHFeMCMixed;
    
    TH1F **fHFEhMCNormalized1D;
    TH1F **fHFEhNormalized1D;
    
    
    
    
    ClassDef(AliHFehpPbTool, 1);
    
};

#endif
