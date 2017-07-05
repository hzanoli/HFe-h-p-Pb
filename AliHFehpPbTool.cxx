#ifndef ALIHFEHPPBTOOL_CXX
#define ALIHFEHPPBTOOL_CXX

#include "AliHFEhpPbTool.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TColor.h"
#include "TRandom3.h"
#include "TLatex.h"


ClassImp(AliHFehpPbTool)

AliHFehpPbTool::AliHFehpPbTool(const Char_t *name, const Char_t *title)
:TNamed(name,title),
fPrintHistograms(kTRUE),
fInputFile(0),
fInputList(0),
fInputFileName("merged.root"),
fConfigurationName(""),
fInclusivepT(0),
fULSpT(0),
fLSpT(0),
fBackgroundpT(0),
fHFepT(0),
fEffTagging(0),
fHadronContamination(0),
fTOFNSigma(3.),
fSameIncEh(0),
fMixedIncEh(0),
fDiHadron(0),
fSameBackULSEh(0),
fSameBackLSEh(0),
fMixedBackULSEh(0),
fMixedBackLSEh(0),
fBackEh(0),
fBackNonIDEh(0),
fHFEh(0),
fHFEhSame(0),
fHFEhMixed(0),
fBackMixedEh(0),
fHFEMC(0),
fHFEhMCNormalized(0),
fHFeMCMixed(0),
fEffHFe(0),
fYieldAS(0),
fYieldNS(0),
fEffCorrectionForElectrons(kFALSE),
fPreMerge(kFALSE),
fTPCNSigmaCenter(0),
fTPCNSigmaSTD(0),
fConfigIndex(0),
fCentralityIndex(0),
fLegendTitle(0),
fMaxDeltaEta(1.6),
fMinDeltaEta(-1.6)
{
    // TH1::AddDirectory(kFALSE);
}



AliHFehpPbTool::AliHFehpPbTool()
:TNamed("AliHFehpPbToolDefault","Default AliHFehpPbTool Analysis"),
fPrintHistograms(kTRUE),
fInputFile(0),
fInputList(0),
fInputFileName("merged.root"),
fConfigurationName(""),
fInclusivepT(0),
fULSpT(0),
fLSpT(0),
fBackgroundpT(0),
fHFepT(0),
fEffTagging(0),
fHadronContamination(0),
fTOFNSigma(3.),
fSameIncEh(0),
fMixedIncEh(0),
fDiHadron(0),
fSameBackULSEh(0),
fSameBackLSEh(0),
fMixedBackULSEh(0),
fMixedBackLSEh(0),
fBackEh(0),
fBackNonIDEh(0),
fHFEh(0),
fBackMixedEh(0),
fHFEhSame(0),
fHFEhMixed(0),
fHFEMC(0),
fHFEhMCNormalized(0),
fHFeMCMixed(0),
fEffHFe(0),
fYieldAS(0),
fYieldNS(0),
fEffCorrectionForElectrons(kFALSE),
fPreMerge(kFALSE),
fTPCNSigmaCenter(0),
fTPCNSigmaSTD(0),
fLegendTitle(0),
fMaxDeltaEta(1.6),
fMinDeltaEta(-1.6)
{
    //TH1::AddDirectory(kFALSE);
}

Bool_t AliHFehpPbTool::ConnectToInputFile(TString FileName, TString ConfigurationName)
{
    //Setting Names
    fConfigurationName = ConfigurationName;
    fInputFileName = FileName;
    
    fInputFile = new TFile(fInputFileName.Data());
    TString DirectoryName = "HFE_h_";
    DirectoryName.Append(fConfigurationName.Data());
    TString ListName = "eh_";
    ListName.Append(ConfigurationName.Data());
    
    TDirectoryFile *dir = (TDirectoryFile*) fInputFile->Get(DirectoryName);
    
    if (!dir)
        printf("No Directory to read the list: %s!\n", ConfigurationName.Data());
    
    fInputList = (TList*) dir->Get(ListName);
    
    if (!fInputList)
    {
        printf("No list to read the histograms!\n");
        return kFALSE;
    }
    else
        return kTRUE;
    
}
void AliHFehpPbTool::Process()
{
    ReadpTHistograms();
    CalculateHadronContamination();
    ProcesspTHistograms();
    ReadAndProcessCorrelationDistributions();
    MergeCorrelationDistributions();
    NormalizeHFeCorrrelation();
    ProjectTo1D();
}

void AliHFehpPbTool::ProcessAtlas()
{
    ReadpTHistograms();
    CalculateHadronContamination();
    ProcesspTHistograms();
    ReadAndProcessCorrelationDistributions();
    MergeCorrelationDistributions();
    NormalizeHFeCorrrelation();
    ProjectTo1DAtlas();
}


void AliHFehpPbTool::ProcessNHFe()
{
    ReadpTHistograms();
    CalculateHadronContamination();
    ProcesspTHistograms();
    ReadAndProcessCorrelationDistributions();
    MakeNonHFEPlots();
    MergeCorrelationDistributions();
    NormalizeHFeCorrrelation();
    ProjectTo1D();
}

void AliHFehpPbTool::ProcessInc()
{
    ReadpTHistograms();
    CalculateHadronContamination();
    ProcesspTHistograms();
    ReadAndProcessCorrelationDistributions();
    MakeIncPlots();
    MergeCorrelationDistributions();
    NormalizeHFeCorrrelation();
    ProjectTo1D();
}


TString AliHFehpPbTool::TrainConfiguration(Int_t pTBin ,
                                           Bool_t Correlation ,
                                           Bool_t ispp ,
                                           Bool_t isMC ,
                                           Double_t ElectronDCAxy ,
                                           Double_t ElectronDCAz,
                                           Double_t HadronDCAxy,
                                           Double_t HadronDCAz ,
                                           Double_t TPCPIDLow ,
                                           Double_t TPCPIDUp ,
                                           Double_t InvariantMassCut ,
                                           Double_t pTCutPartner ,
                                           Double_t MultiplicityLow,
                                           Double_t MultiplicityUp,
                                           Double_t HadronPtCutLow,
                                           Double_t HadronPtCutUp,
                                           Double_t EtaCutLow,
                                           Double_t EtaCutUp,
                                           Double_t NonHFEangleCut,
                                           Int_t NHitsITS,
                                           Int_t SPDLayers,
                                           Int_t TPCNCluster,
                                           Int_t TPCNClusterPartner ,
                                           Int_t TPCNClusterPID)
{
    TString Name = Form("%d,%d,%d,%d,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d,%d",pTBin, Correlation, ispp, isMC,   ElectronDCAxy,ElectronDCAz,HadronDCAxy,HadronDCAz,TPCPIDLow,TPCPIDUp,InvariantMassCut,pTCutPartner, MultiplicityLow, MultiplicityUp, HadronPtCutLow, HadronPtCutUp, EtaCutLow, EtaCutUp, NonHFEangleCut, NHitsITS, SPDLayers, TPCNCluster, TPCNClusterPartner, TPCNClusterPID);
    return Name;
}

Bool_t AliHFehpPbTool::ReadpTHistograms()
{
    //Reading Output from file
    
    fInclusivepT = (TH1F*) ( (TH1F*)fInputList->FindObject("fPtElec_Inc")  )->Clone("fInclusivepT");
    
    if (!fInclusivepT)
        printf("No Inclusive histogram!\n");
    
    fULSpT =  (TH1F*) ( (TH1F*)fInputList->FindObject("fPtElec_ULS")  )->Clone("fULSpT");
    if (!fULSpT)
        printf("No ULS histogram!\n");
    
    fLSpT = (TH1F*) ( (TH1F*)fInputList->FindObject("fPtElec_LS")  )->Clone("fLSpT");
    if (!fLSpT)
        printf("No LS histogram!\n");
    
    fBackgroundpT = (TH1F*) fULSpT->Clone("fBackgroundpT");
    
    fBackgroundpT->Add(fLSpT,-1);
    
    if (fBackgroundpT && fInclusivepT)
        return kTRUE;
    else
        return kFALSE;
    
}

void AliHFehpPbTool::SetEfficiencyTaggingToConstant(Double_t Value)
{
    //Assume no EfficiencyHistogram
    fEffTagging = new TH1F("fEffTagging", "Tagging Efficiency (set to constant", fpTBins.GetSize()-1, fpTBins.GetArray());
    
    for (Int_t i = 1; i <= fpTBins.GetSize()-1 ; i++)
    {
        fEffTagging->SetBinContent(i, Value);
        fEffTagging->SetBinError(i, Value/1000000000000);
    }
    
}


Bool_t AliHFehpPbTool::CompareMixedEventDistributions()
{
    ReadAndProcessCorrelationDistributions();
    TH2F **RatioULSLS;
    TH2F **RatioIncULS;
    RatioULSLS = new TH2F *[fpTBins.GetSize()-1];
    RatioIncULS = new TH2F *[fpTBins.GetSize()-1];
    TCanvas *Canvas = new TCanvas("CompareMixed","CompareMixed", 1200,900);
    Canvas->Divide(4,3);
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1 ; i++)
    {
        /*
         RatioIncULS[i] = (TH2F*) (fMixedIncEh[i]->Clone("RatioIncULS"));
         RatioIncULS[i]->Divide(fMixedIncEh[i],fMixedBackULSEh[i]);
         RatioIncULS[i]->Scale(fMixedBackULSEh[i]->GetBinContent(fMixedBackULSEh[i]->FindBin(0,0)) / fMixedIncEh[i]->GetBinContent(fMixedIncEh[i]->FindBin(0,0)) );
         Canvas->cd(i+1);
         RatioIncULS[i]->GetYaxis()->SetRangeUser(-1.6,1.6);
         RatioIncULS[i]->Draw("surf1");
         */
        
        
        RatioULSLS[i] = (TH2F*) (fMixedBackULSEh[i]->Clone("RatioULSLS"));
        RatioULSLS[i]->Add(fMixedBackLSEh[i],-1);
        RatioULSLS[i]->Scale(1./RatioULSLS[i]->GetBinContent(RatioULSLS[i]->FindBin(0,0)));
        RatioULSLS[i]->Divide(RatioULSLS[i],fMixedIncEh[i]);
        RatioULSLS[i]->Scale(fMixedIncEh[i]->GetBinContent(fMixedIncEh[i]->FindBin(0,0)));
        Canvas->cd(i+1);
        RatioULSLS[i]->GetYaxis()->SetRangeUser(-1.6,1.6);
        RatioULSLS[i]->Draw("surf1");
        
    }
    
    
}

void AliHFehpPbTool::PreMerge()
{
    Double_t pTBinsCorrelation[] = {0.5,0.75,1.0,1.25,1.5,2,2.5,3,4,6};
    fpTBins.Set((Int_t) sizeof(pTBinsCorrelation)/sizeof(pTBinsCorrelation[0]),pTBinsCorrelation);
    fPreMerge = kTRUE;
}

Bool_t AliHFehpPbTool::ReadAndProcessCorrelationDistributions(Int_t RebinX , Int_t RebinY )
{
    printf("Reading and processing correlation \n");
    //TCanvas *CorrelationHist = new TCanvas("CorrelationHist ","CorrelationHist", 1200,900 );
    
    // CorrelationHist->Divide(4,3);
    
    
    fSameIncEh =  new TH2F *[fpTBins.GetSize()-1];
    fMixedIncEh = new TH2F *[fpTBins.GetSize()-1];
    
    fDiHadron = new TH2F *[fpTBins.GetSize()-1];
    
    fSameBackULSEh = new TH2F *[fpTBins.GetSize()-1];
    fSameBackLSEh = new TH2F *[fpTBins.GetSize()-1];
    
    fMixedBackULSEh = new TH2F *[fpTBins.GetSize()-1];
    fMixedBackLSEh = new TH2F *[fpTBins.GetSize()-1];
    
    fSameBackNoPULSEh = new TH2F *[fpTBins.GetSize()-1];
    fSameBackNoPLSEh = new TH2F *[fpTBins.GetSize()-1];
    
    fBackEh = new TH2F *[fpTBins.GetSize()-1];
    fBackNonIDEh = new TH2F *[fpTBins.GetSize()-1];
    
    fBackMixedEh = new TH2F *[fpTBins.GetSize()-1];
    
    fHFEh =new TH2F *[fpTBins.GetSize()-1];
    fHFEhSame =new TH2F *[fpTBins.GetSize()-1];
    fHFEhMixed =new TH2F *[fpTBins.GetSize()-1];
    
    
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1 ; i++)
    {
        //Inclusive Distributions
        fSameIncEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_Inc%d",i));
        fMixedIncEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_Inc_EM%d",i));
        fDiHadron[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_Inc_DiHadron%d",i));
        
        
        if ( (i == (fpTBins.GetSize() - 2))  && (fPreMerge == kTRUE ))
        {
            TH2F *TempSameEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_Inc%d",i+1));
            TH2F *TempMixedEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_Inc_EM%d",i+1));
            TH2F *TempDiHadron = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_Inc_DiHadron%d",i+1));
            fSameIncEh[i]->Add(TempSameEh);
            fMixedIncEh[i]->Add(TempMixedEh);
            fDiHadron[i]->Add(TempDiHadron);
        }
        fSameIncEh[i]->SetTitle(Form("%1.2f < p^{e}_{T} < %1.2f GeV/c", fpTBins.At(i), fpTBins.At(i+1)));
        fMixedIncEh[i]->SetTitle(Form("%1.2f < p^{e}_{T} < %1.2f GeV/c", fpTBins.At(i), fpTBins.At(i+1)));
        
        
        fSameIncEh[i]->RebinY(RebinY);
        fSameIncEh[i]->RebinX(RebinX);
        fMixedIncEh[i]->RebinY(RebinY);
        fMixedIncEh[i]->RebinX(RebinX);
        
        //Hadron Contamination
        fDiHadron[i]->RebinY(RebinY);
        fDiHadron[i]->RebinX(RebinX);
        
        TH1F *DihadronpT = (TH1F*)fInputList->FindObject("fPtTrigger_Inc");
        DihadronpT = (TH1F*) DihadronpT->Rebin(fpTBins.GetSize()-1, "DihadronpTRebin", fpTBins.GetArray());
        
        
        fDiHadron[i]->Scale( fInclusivepT->GetBinContent(i+1)/DihadronpT->GetBinContent(i+1) * (1-fHadronContamination->GetBinContent(i+1)) );
        fSameIncEh[i]->Add(fDiHadron[i],-1);
        
        //Background Same
        fSameBackULSEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_ULS_Weight%d",i));
        fSameBackLSEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_LS_Weight%d",i));
        fSameBackNoPULSEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_ULS_NoP_Weight%d",i));
        fSameBackNoPLSEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_LS_NoP_Weight%d",i));
        fMixedBackULSEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_ULS_Weight_EM%d",i));
        fMixedBackLSEh[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_LS_Weight_EM%d",i));
        
        
        if ( (i == (fpTBins.GetSize() - 2))  && (fPreMerge == kTRUE ))
        {
            TH2F* TempfSameBackULSEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_ULS_Weight%d",i+1));
            TH2F* TempfSameBackLSEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_LS_Weight%d",i+1));
            TH2F* TempfSameBackNoPULSEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_ULS_NoP_Weight%d",i+1));
            TH2F* TempfSameBackNoPLSEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_LS_NoP_Weight%d",i+1));
            TH2F* TempfMixedBackULSEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_ULS_Weight_EM%d",i+1));
            TH2F* TempfMixedBackLSEh = (TH2F*) fInputList->FindObject(Form("fCEtaPhi_LS_Weight_EM%d",i+1));
            
            fSameBackULSEh[i]->Add(TempfSameBackULSEh);
            fSameBackLSEh[i]->Add(TempfSameBackLSEh);
            fSameBackNoPULSEh[i]->Add(TempfSameBackNoPULSEh);
            fSameBackNoPLSEh[i]->Add(TempfSameBackNoPLSEh);
            fMixedBackULSEh[i]->Add(TempfMixedBackULSEh);
            fMixedBackLSEh[i]->Add(TempfMixedBackLSEh);
            
        }
        
        
        
        fSameBackULSEh[i]->RebinY(RebinY);
        fSameBackULSEh[i]->RebinX(RebinX);
        fSameBackLSEh[i]->RebinY(RebinY);
        fSameBackLSEh[i]->RebinX(RebinX);
        fMixedBackULSEh[i]->RebinY(RebinY);
        fMixedBackULSEh[i]->RebinX(RebinX);
        fMixedBackLSEh[i]->RebinY(RebinY);
        fMixedBackLSEh[i]->RebinX(RebinX);
        fSameBackNoPULSEh[i]->RebinY(RebinY);
        fSameBackNoPULSEh[i]->RebinX(RebinX);
        fSameBackNoPLSEh[i]->RebinY(RebinY);
        fSameBackNoPLSEh[i]->RebinX(RebinX);
        
        printf("Processing distributions \n");
        
        
        fBackEh[i] = (TH2F*) fSameBackULSEh[i]->Clone(Form("fBackEh%d",i));
        fBackNonIDEh[i] = (TH2F*) fSameBackNoPULSEh[i]->Clone(Form("fBackNonIDEh%d",i));
        fBackMixedEh[i]  = (TH2F*) fMixedBackULSEh[i]->Clone(Form("fBackMixedEh%d",i));
        
        //ID Back
        
        fBackEh[i]->Add(fSameBackLSEh[i],-1);
        
        
        //NonID Back
        
        fBackNonIDEh[i]->Add(fSameBackNoPLSEh[i],-1);
        
        fBackNonIDEh[i]->Scale(1./(fEffTagging->GetBinContent(i+1)) - 1.);
        fBackMixedEh[i]->Add(fMixedBackLSEh[i],-1);
        fBackMixedEh[i]->Scale(1./(fEffTagging->GetBinContent(i+1)));
        
        //SemiInclusive Distribution
        fHFEhSame[i] = (TH2F*) fSameIncEh[i]->Clone(Form("fHFEhSame%d",i));
        fHFEhSame[i]->Add(fBackEh[i],-1);
        
        
        
        for (Int_t ix = 1; ix <= fHFEhSame[i]->GetNbinsX(); ix++ )
        {
            for (Int_t iy = 1; iy <= fHFEhSame[i]->GetNbinsY(); iy++ )
            {
                fHFEhSame[i]->SetBinError(ix,iy,  TMath::Sqrt( TMath::Abs(pow(fSameIncEh[i]->GetBinError(ix,iy),2) -  pow(fBackEh[i]->GetBinError(ix,iy),2) )) );
            }
        }
        
        
        
        //HFe Same Distibution
        fHFEhSame[i]->Add(fBackNonIDEh[i],-1);
        
        //HFe Mixed
        fHFEhMixed[i] = (TH2F*) fMixedIncEh[i]->Clone(Form("fHFEhMixed%d",i));
        fHFEhMixed[i]->Add(fBackMixedEh[i],-1);
        
        
        Double_t NormalizationMixed = (fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.),fHFEhMixed[i]->GetYaxis()->FindBin(0.)) + fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.)-1,fHFEhMixed[i]->GetYaxis()->FindBin(0.)) + fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.),fHFEhMixed[i]->GetYaxis()->FindBin(0.) -1) + fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.)-1,fHFEhMixed[i]->GetYaxis()->FindBin(0.)-1))/4. ;
        
        
        //Finite Bin Correction
        NormalizationMixed = NormalizationMixed / (1 - 1./1.6 * fHFEhMixed[i]->GetYaxis()->GetBinCenter( fHFEhMixed[i]->GetYaxis()->FindBin(0.) ) );
        
        
        fHFEhMixed[i]->Scale(1./NormalizationMixed);
        fHFEh[i] = (TH2F*) fHFEhSame[i]->Clone(Form("fHFEh%d",i));
        fHFEh[i]->Divide(fHFEh[i],fHFEhMixed[i]);
        
        if (fEffCorrectionForElectrons)
        {
            fHFEh[i]->Scale(1./fEffHFe->GetBinContent(i+1));
            fHFEhSame[i]->Scale(1./fEffHFe->GetBinContent(i+1));
            fHFEhMixed[i]->Scale(NormalizationMixed);
        }
        
        // CorrelationHist->cd(i+1);
        
        //fSameIncEh[i]->Draw("lego2");
        // fHFEh[i]->Draw("lego2");
        
    }
    
    //if (fHFEh[fpTBins.GetSize() - 2])
    delete fInputFile;
    return kTRUE;
    //else
    //  return kFALSE;
    
}

void AliHFehpPbTool::MakeNonHFEPlots()
{
    //change pT plots -> instead of HFE <-> Non-HFE
    fBackgroundpT->Multiply(fEffTagging); //Restore Background
    fHFepT = fBackgroundpT;
    
    //Change HFE-h -> NonHFE-h
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1 ; i++)
    {
        //Restore the efficiency scaling
        fBackNonIDEh[i]->Scale(1./(1./(fEffTagging->GetBinContent(i+1)) - 1.));
        fBackMixedEh[i]->Scale(1.*(fEffTagging->GetBinContent(i+1)));
        
        Double_t NormalizationMixed = (fBackMixedEh[i]->GetBinContent(fBackMixedEh[i]->GetXaxis()->FindBin(0.),fBackMixedEh[i]->GetYaxis()->FindBin(0.)) + fBackMixedEh[i]->GetBinContent(fBackMixedEh[i]->GetXaxis()->FindBin(0.)-1,fBackMixedEh[i]->GetYaxis()->FindBin(0.)) + fBackMixedEh[i]->GetBinContent(fBackMixedEh[i]->GetXaxis()->FindBin(0.),fBackMixedEh[i]->GetYaxis()->FindBin(0.) -1) + fBackMixedEh[i]->GetBinContent(fBackMixedEh[i]->GetXaxis()->FindBin(0.)-1,fBackMixedEh[i]->GetYaxis()->FindBin(0.)-1))/4. ;
        
        
        //Finite Bin Correction
        NormalizationMixed = NormalizationMixed / (1 - 1./1.6 * fBackMixedEh[i]->GetYaxis()->GetBinCenter( fBackMixedEh[i]->GetYaxis()->FindBin(0.) ) );
        
        fBackMixedEh[i]->Scale(1./NormalizationMixed);
        
        fBackNonIDEh[i]->Divide(fBackMixedEh[i]);
        
        fHFEh[i] = fBackNonIDEh[i];
        
    }
}

void AliHFehpPbTool::MakeIncPlots()
{
    //change pT plots -> instead of HFE <-> Non-HFE
    fHFepT = fInclusivepT;
    
    //Change HFE-h -> NonHFE-h
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1 ; i++)
    {
        //Restore the efficiency scaling
        Double_t NormalizationMixed = (fMixedIncEh[i]->GetBinContent(fMixedIncEh[i]->GetXaxis()->FindBin(0.),fMixedIncEh[i]->GetYaxis()->FindBin(0.)) + fMixedIncEh[i]->GetBinContent(fMixedIncEh[i]->GetXaxis()->FindBin(0.)-1,fMixedIncEh[i]->GetYaxis()->FindBin(0.)) + fMixedIncEh[i]->GetBinContent(fMixedIncEh[i]->GetXaxis()->FindBin(0.),fMixedIncEh[i]->GetYaxis()->FindBin(0.) -1) + fMixedIncEh[i]->GetBinContent(fMixedIncEh[i]->GetXaxis()->FindBin(0.)-1,fMixedIncEh[i]->GetYaxis()->FindBin(0.)-1))/4. ;
        
        
        //Finite Bin Correction
        NormalizationMixed = NormalizationMixed / (1 - 1./1.6 * fMixedIncEh[i]->GetYaxis()->GetBinCenter( fMixedIncEh[i]->GetYaxis()->FindBin(0.) ) );
        
        fMixedIncEh[i]->Scale(1./NormalizationMixed);
        
        fSameIncEh[i]->Divide(fMixedIncEh[i]);
        
        fHFEh[i] = fSameIncEh[i];
        
    }
}





Bool_t AliHFehpPbTool::MergeCorrelationDistributions()
{
    if (fpTBinsResults.GetSize() == fpTBins.GetSize())
    {
        printf("Not merging, since results and input have the same binning\n");
        fHFEhNormalized = new TH2F *[fpTBinsResults.GetSize() -1];
        
        for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
            fHFEhNormalized[i] = (TH2F*) fHFEh[i]->Clone(Form("fHFEhNormalized%d",i));
    }
    else {
        fHFEhNormalized = new TH2F *[fpTBinsResults.GetSize() -1];
        fHFEhSameNormalized = new TH2F *[fpTBinsResults.GetSize() -1];
        fHFEhMixedNormalized = new TH2F *[fpTBinsResults.GetSize() -1];
        //Assuming both always start at the same bin
        Int_t lastj = 0;
        for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
        {
            fHFEhNormalized[i] = (TH2F*) fHFEh[0]->Clone(Form("fHFEhNormalized%d",i));
            fHFEhSameNormalized[i] = (TH2F*) fHFEh[0]->Clone(Form("fHFEhSameNormalized%d",i));
            fHFEhMixedNormalized[i] = (TH2F*) fHFEh[0]->Clone(Form("fHFEhMixedNormalized%d",i));
            
            fHFEhNormalized[i]->Reset();
            fHFEhSameNormalized[i]->Reset();
            fHFEhMixedNormalized[i]->Reset();
            
            fHFEhNormalized[i]->GetSumw2()->Set(0);
            fHFEhSameNormalized[i]->GetSumw2()->Set(0);
            fHFEhMixedNormalized[i]->GetSumw2()->Set(0);
            
            fHFEhNormalized[i]->Sumw2();
            fHFEhSameNormalized[i]->Sumw2();
            fHFEhMixedNormalized[i]->Sumw2();
            
            for (Int_t j = lastj ; j < fpTBins.GetSize()-1 ; j++)
            {
                if (fpTBinsResults.At(i+1) >= fpTBins.At(j+1))
                {
                    fHFEhNormalized[i]->Add(fHFEh[j]);
                    fHFEhSameNormalized[i]->Add(fHFEhSame[j]);
                    fHFEhMixedNormalized[i]->Add(fHFEhMixed[j]);
                }
                else
                {
                    lastj = j;
                    break;
                }
            }
            
            fHFEhNormalized[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f",fpTBinsResults.At(i), fpTBinsResults.At(i+1)));
            fHFEhSameNormalized[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f",fpTBinsResults.At(i), fpTBinsResults.At(i+1)));
            fHFEhMixedNormalized[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f",fpTBinsResults.At(i), fpTBinsResults.At(i+1)));

        }
    }
    
    if (fHFEhNormalized[0])
        return kTRUE;
    else
        return kFALSE;
}


/*
 Bool_t AliHFehpPbTool::MergeCorrelationDistributions()
 {
 return kTRUE;
 }
 */


Bool_t AliHFehpPbTool::ReadMCDistributions(Int_t RebinX , Int_t RebinY )
{
    fHFEMC = new TH2F *[fpTBins.GetSize()-1];
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1 ; i++)
    {
        fHFEMC[i] = (TH2F*) fInputList->FindObject(Form("fCEtaPhiCutHFe%d",i));
        fHFEMC[i]->RebinY(RebinY);
        fHFEMC[i]->RebinX(RebinX);
        
    }
    fHFEMC[0]->Draw("surf1");
}

Bool_t AliHFehpPbTool::MergeMCCorrelationDistributions()
{
    fHFEhMCNormalized = new TH2F *[fpTBinsResults.GetSize() -1];
    
    //Assuming both always start at the same bin
    Int_t lastj = 0;
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        
        fHFEhMCNormalized[i] = (TH2F*) fHFEMC[0]->Clone(Form("fHFEhMCNormalized%d",i));
        fHFEhMCNormalized[i]->Reset();
        printf("================================\n");
        printf("Merging bin from %1.2f to %1.2f\n",fpTBinsResults.At(i), fpTBinsResults.At(i+1) );
        
        for (Int_t j = lastj ; j < fpTBins.GetSize()-1 ; j++)
        {
            if (fpTBinsResults.At(i+1) >= fpTBins.At(j+1))
            {
                printf("Add: Merge %1.2f to %1.2f in the bin %1.2f to %1.2f \n", fpTBins.At(j), fpTBins.At(j+1), fpTBinsResults.At(i), fpTBinsResults.At(i+1));
                fHFEhMCNormalized[i]->Add(fHFEMC[j]);
            }
            else
            {
                printf("It was merged from %1.2f to %1.2f and final value should be %1.2f\n",fpTBinsResults.At(i), fpTBins.At(j),fpTBinsResults.At(i+1));
                lastj = j;
                break;
            }
        }
        fHFEhMCNormalized[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f",fpTBinsResults.At(i), fpTBinsResults.At(i+1)));
    }
}


Bool_t AliHFehpPbTool::MakeMCTruthMixed()
{
    fHFeMCMixed = new TH2F *[fpTBinsResults.GetSize() -1];
    fHFeMCMixed[0] = (TH2F*) fHFEhMCNormalized[0]->Clone("fHFeMCMixed0");
    fHFeMCMixed[0]->Reset();
    
    for (Int_t ix = 1; ix <= fHFeMCMixed[0]->GetNbinsX() ; ix++)
    {
        for (Int_t iy = 1; iy <= fHFeMCMixed[0]->GetNbinsY() ; iy++)
        {
            fHFeMCMixed[0]->SetBinContent(ix,iy,1.-1.0/1.6* TMath::Abs(fHFeMCMixed[0]->GetYaxis()->GetBinCenter(iy)) );
            fHFeMCMixed[0]->SetBinError(ix,iy,0.0000000000000000000001);
        }
    }
    fHFeMCMixed[0]->Draw("surf1");
    //fHFeMCMixed[0]->Scale(1./fHFeMCMixed[0]->GetBinContent(fHFeMCMixed[0]->GetXaxis()->FindBin(0.),fHFeMCMixed[0]->GetYaxis()->FindBin(0.) ) );
    
    for (Int_t i = 1 ; i < fpTBinsResults.GetSize() -1; i++)
        fHFeMCMixed[i] = (TH2F*) fHFeMCMixed[0]->Clone(Form("fHFeMCMixed%d",i));
    
}

Bool_t AliHFehpPbTool::CorrectMCDistributionsByMixing()
{
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
        fHFEhMCNormalized[i]->Divide(fHFeMCMixed[i]);
    
}


Bool_t AliHFehpPbTool::CorrelationCT()
{
    fRatio2DMCCT = new TH2F *[fpTBinsResults.GetSize() -1];
    TCanvas **Canvas = new TCanvas *[fpTBinsResults.GetSize() -1];
    
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        //fHFEhNormalized[i]->RebinY(2);
        //fHFEhMCNormalized[i]->RebinY(2);
        Canvas[i] = new TCanvas(Form("MCCT%d",i),Form("MCCT%d",i),400,300);
        Canvas[i]->cd();
        fRatio2DMCCT[i] = (TH2F*) fHFEhNormalized[i]->Clone(Form("RatioDataMC%d",i));
        fRatio2DMCCT[i]->Divide(fRatio2DMCCT[i], fHFEhMCNormalized[i]);
        fRatio2DMCCT[i]->GetZaxis()->SetTitle("Rec/MC");
        fRatio2DMCCT[i]->GetYaxis()->SetRangeUser(-1.6,1.6);
        fRatio2DMCCT[i]->GetYaxis()->SetTitleOffset(1.8);
        fRatio2DMCCT[i]->GetXaxis()->SetTitleOffset(1.6);
        fRatio2DMCCT[i]->GetZaxis()->SetTitleOffset(1.4);
        fRatio2DMCCT[i]->Draw("surf1");
        Canvas[i]->SaveAs(Form("MC_test_%d.pdf",i));
    }
}



Bool_t AliHFehpPbTool::ProcesspTHistograms()
{
    //Rebin
    
    printf("========================\n");
    printf("Processing pT Distributions\n");
    printf("========================\n");
    
    fInclusivepT = (TH1F*) fInclusivepT->Rebin(fpTBins.GetSize()-1, "fInclusivepT", fpTBins.GetArray());
    
    //Remove Hadron contamination (:
    fInclusivepT->Multiply(fInclusivepT,fHadronContamination);
    
    fBackgroundpT = (TH1F*) fBackgroundpT->Rebin(fpTBins.GetSize()-1, "fBackgroundpT", fpTBins.GetArray());
    
    //TH1F *SemiInclusive = (TH1F*)fInclusivepT->Clone("SemiInclusive");
    //fInclusivepT
    fBackgroundpT->Divide(fEffTagging);
    
    //Remove sec. electrons
    fHFepT = (TH1F*)fInclusivepT->Clone("fHFEpT");
    fHFepT->Add(fBackgroundpT,-1);
    
    //Calculate TPC-TOF pT distribution
    if (fEffCorrectionForElectrons)
        fHFepT->Divide(fEffHFe);
    
    
    /*
     if (kTRUE)
     {
     for (Int_t i = 1; i <= fHFepT->GetNbinsX(); i++)
     {
     Double_t RootError = fHFepT->GetBinError(i);
     Double_t Minus = TMath::Sqrt( fBackgroundpT->GetBinError(i)*fBackgroundpT->GetBinError(i) - fInclusivepT->GetBinError(i)* fInclusivepT->GetBinError(i));
     Double_t Plus = TMath::Sqrt( TMath::Abs( fBackgroundpT->GetBinError(i)* fBackgroundpT->GetBinError(i) + fInclusivepT->GetBinError(i)* fInclusivepT->GetBinError(i)));
     printf("Def Error = %f    Minus Error = %f    Plus Error =   %f   Ratio Minus/old = %f  Ratio Plus/old = %f  \n", RootError, Minus,Plus,Minus/RootError, Plus/RootError );
     }
     }
     
     */
    if (fHFepT)
        return kTRUE;
    
    return kFALSE;
    
}


//Dont need to use MC information
Bool_t AliHFehpPbTool::CalculateHadronContamination()
{
    printf("Calculating hadron contamination \n");
    
    TH1F* ContamiantionUsingCuts[7];
    
    for (Int_t a = 0; a < 7; a++) {
        ContamiantionUsingCuts[a] = new TH1F("ContamiantionUsingCuts","Hadron Cont.; p (GeV/c) ", fpTBins.GetSize()-1, fpTBins.GetArray());
    }
    
    //TCanvas *ContaminationPlots = new TCanvas("Contamination", "Contamination", 400,300);
    
    Double_t MinimumValueToUseRightSideFit = 0.5;
    Double_t MaximumValueToUseRightSideFit = 1.0;
    
    
    Bool_t lPerformCutStudies = kTRUE;
    Bool_t DrawG = kTRUE;
    Bool_t FixElectronposition = kFALSE;
    
    TH2F **TPCTOFNsigma = new TH2F *[fpTBins.GetSize()];
    TH1F **TPCNsigma = new TH1F *[fpTBins.GetSize()];
    
    fHadronContamination = new TH1F("fHadronContamination","Hadron Contamination", fpTBins.GetSize()-1, fpTBins.GetArray());
    
    fTPCNSigmaCenter = new TH1F("fTPCNSigmaCenter","Mean TPCNSigma Values for electrons", fpTBins.GetSize()-1, fpTBins.GetArray());
    fTPCNSigmaSTD = new TH1F("fTPCNSigmaSTD","Mean TPCNSigma Values for electrons", fpTBins.GetSize()-1, fpTBins.GetArray());
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1; i++)
    {
        TPCTOFNsigma[i] = (TH2F*)fInputList->FindObject(Form("fTOFTPCnsigma_pt%d",i));
        TPCNsigma[i] = (TH1F*) TPCTOFNsigma[i]->ProjectionY(Form("TPCnsigma_pt%d",i), TPCTOFNsigma[i]->GetXaxis()->FindBin(-fTOFNSigma), TPCTOFNsigma[i]->GetXaxis()->FindBin(fTOFNSigma)) ;
        
        if  (DrawG)
        {
            TCanvas *Plot2DTPCTOF = new TCanvas(Form("Plot2DTPCTOF%d",i),Form("Plot2DTPCTOF%d",i), 900,600 );
            Plot2DTPCTOF->SetLogz();
            TPCTOFNsigma[i]->Draw("colz");
            Plot2DTPCTOF->SaveAs(Form("Plot2DTPCTOF_Cent_%d_pT%d_Config%d.png",fCentralityIndex, i,fConfigIndex));
        }
        
        TPCNsigma[i]->GetYaxis()->SetTitle("Entries");
        
        TPCNsigma[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f ",fpTBins.At(i),fpTBins.At(i+1)));
        
        //TPCNsigma[i]->Rebin();
        if (fpTBins.At(i) >= 4.0)
            TPCNsigma[i]->Rebin(2);
        
        
        //Use Pion as 2 gaussians
        
        TF1 *pion = new TF1("Pions", "gausn(0)+gausn(3)",-12,8);
        
        pion->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetMaximumBin()), TPCNsigma[i]->GetXaxis()->GetBinLowEdge(TPCNsigma[i]->GetMaximumBin()), 0.8, 0.1 * TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetMaximumBin()), TPCNsigma[i]->GetXaxis()->GetBinLowEdge(TPCNsigma[i]->GetMaximumBin()), 0.8 );
        
        Double_t PionFitMin,PionFitMax(-2);
        
        PionFitMin = TPCNsigma[i]->GetXaxis()->GetBinLowEdge(TPCNsigma[i]->GetMaximumBin()) - 0.7;
        
        TPCNsigma[i]->Fit(pion,"QN","", PionFitMin, PionFitMax);
        
        TF1 *elec = new TF1("Electrons", "gausn(0)", -10, 10);
        elec->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(0.)), 0, 1);
        elec->SetParLimits(1, -1,1);
        TPCNsigma[i]->Fit(elec,"QN","",-0.5, 3);
        
        TF1 *Kaon = new TF1("Kaon", "gausn(0)",-10,10);
        
        if ( (fpTBins.At(i) >= MinimumValueToUseRightSideFit) && (fpTBins.At(i) < MaximumValueToUseRightSideFit) )
            Kaon->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(7.)), 7, 2);
        else
            Kaon->FixParameter(0,0);
        
        TPCNsigma[i]->Fit(Kaon,"QN","",5.5, 10);
        
        
        TF1 *TotalFit = new TF1("Total Fit", "gausn(0)+gausn(3)+gausn(6) +gaus(9)",PionFitMin,10); //Elec + pion (2 gaus) + [kaon (1 gaus) -- only until 2 GeV/c
        
        Double_t Parameters[12] = {elec->GetParameter(0), elec->GetParameter(1), elec->GetParameter(2), pion->GetParameter(0), pion->GetParameter(1), pion->GetParameter(2), pion->GetParameter(3), pion->GetParameter(4), pion->GetParameter(5), Kaon->GetParameter(0), Kaon->GetParameter(1), Kaon->GetParameter(2)};
        
        TotalFit->SetParameters(Parameters);
        Double_t TotalFitMax = 10;
        
        TotalFit->SetParLimits(1,-1,1);
        TotalFit->SetParLimits(2,0.5,1.5);
        
        TotalFit->SetParLimits(10,0,10);
        TotalFit->SetParLimits(11,0,10);
        
        if (! ((fpTBins.At(i) >= MinimumValueToUseRightSideFit) && (fpTBins.At(i) < MaximumValueToUseRightSideFit)) )
        {
            TotalFit->FixParameter(11,0);
            TotalFit->FixParameter(10,0);
            TotalFit->FixParameter(9,0);
            TotalFitMax = 2.;
        }
        
        
        
        if (i == 0 || i == 1)
        {
            TotalFit->FixParameter(11,Kaon->GetParameter(2));
            TotalFit->FixParameter(10,Kaon->GetParameter(1));
            TotalFit->FixParameter(9,Kaon->GetParameter(0));
            
        }
        
        if (fpTBins.At(i) == 4)
        {
            TotalFit->FixParameter(6,0);
            TotalFit->FixParameter(7,0);
            TotalFit->FixParameter(8,0);
            
        }
        
        
        
        
        TPCNsigma[i]->Fit(TotalFit,"LQI0","", PionFitMin, TotalFitMax);
        
        elec->SetParameters(TotalFit->GetParameter(0), TotalFit->GetParameter(1), TotalFit->GetParameter(2));
        
        pion->SetParameters(TotalFit->GetParameter(3), TotalFit->GetParameter(4), TotalFit->GetParameter(5),TotalFit->GetParameter(6), TotalFit->GetParameter(7), TotalFit->GetParameter(8));
        
        Kaon->SetParameters(TotalFit->GetParameter(9), TotalFit->GetParameter(10), TotalFit->GetParameter(11));
        
        if(DrawG)
        {
            TCanvas *HadronContamination = new TCanvas(Form("HadronContamination%d",i),Form("HadronContamination%d",i), 1.5*400,1.5*300);
            HadronContamination->SetLogy();
           // HadronContamination->SetGridx();
            //HadronContamination->SetGridy();
            TPCNsigma[i]->SetTitle("");
            TPCNsigma[i]->Draw();
            TPCNsigma[i]->SetLineColor(1);
            TPCNsigma[i]->SetMarkerColor(1);
            TPCNsigma[i]->SetMarkerSize(0.3);
            TPCNsigma[i]->SetMarkerStyle(24);
            
            TPCNsigma[i]->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} - <TPC d#it{E}/d#it{x}>|_{e} (#sigma)");
            
            TPCNsigma[i]->GetYaxis()->SetTitleSize(0.045);
            TPCNsigma[i]->GetXaxis()->SetTitleSize(0.045);
            
            TPCNsigma[i]->GetYaxis()->SetRangeUser(0.5,10E7);
            
            
            TotalFit->SetLineColor(3);
            TotalFit->Draw("same");
            
            TH1F* Ratio = (TH1F*) TPCNsigma[i]->Clone("Data/Fit");
            Ratio->SetMarkerColor(kAzure);
            Ratio->SetLineColor(kAzure);
            Ratio->Divide(TotalFit);
            Ratio->GetXaxis()->SetRangeUser(-5.,3.);
            Ratio->Draw("same");
            
            Ratio->GetXaxis()->SetRangeUser(PionFitMin,3);
            pion->SetLineColor(9);
            pion->Draw("same");
            pion->SetLineStyle(2);
            
            elec->Draw("same");
            elec->SetLineStyle(2);
            elec->SetLineColor(2);
            elec->SetLineWidth(3);
            
            Kaon->SetLineStyle(2);
            Kaon->SetLineColor(kBlue);
            if (((fpTBins.At(i) >= MinimumValueToUseRightSideFit) && (fpTBins.At(i) < MaximumValueToUseRightSideFit)))
                Kaon->Draw("same");
            
            TLatex *   tex = new TLatex(0.14,0.86,"ALICE Preliminary");
            tex->SetNDC();
            tex->SetTextAlign(12);
            tex->SetTextFont(42);
            tex->SetTextSize(0.03782506);
            tex->SetLineWidth(2);
            tex->Draw();
            tex = new TLatex(0.14,0.81,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
            tex->SetNDC();
            tex->SetTextAlign(12);
            tex->SetTextFont(42);
            tex->SetTextSize(0.03782506);
            tex->SetLineWidth(2);
            tex->Draw();
            tex = new TLatex(0.14,0.75,"0-20% ZNA class");
            tex->SetNDC();
            tex->SetTextAlign(12);
            tex->SetTextFont(42);
            tex->SetTextSize(0.03782506);
            tex->SetLineWidth(2);
            tex->Draw();
            tex = new TLatex(0.14,0.69,"|#eta| < 0.8");
            tex->SetNDC();
            tex->SetTextAlign(12);
            tex->SetTextFont(42);
            tex->SetTextSize(0.03782506);
            tex->SetLineWidth(2);
            tex->Draw();
            tex = new TLatex(0.47,0.83,Form("%1.1f < #it{p}_{T} < %1.1f GeV/#it{c}",fpTBins.At(i),fpTBins.At(i+1)));
            tex->SetNDC();
            tex->SetTextAlign(12);
            tex->SetTextFont(42);
            tex->SetTextSize(0.03782506);
            tex->SetLineWidth(2);
            tex->Draw();

            
            
            TLegend leg(0.6677852,0.5614657,0.8892617,0.7813239);
            leg.SetTextSize(0.03782506);
            leg.SetBorderSize(0);
            leg.AddEntry(TPCNsigma[i],"Data", "lp");
            
            
            //leg.SetHeader(fLegendTitle.Data());
            leg.AddEntry(pion,"Pion", "l");
            leg.AddEntry(elec,"Electron", "l");
            leg.AddEntry(TotalFit, "Total Fit", "l");
            leg.AddEntry(Ratio, "Data/Fit", "lp");
            if (((fpTBins.At(i) >= MinimumValueToUseRightSideFit) && (fpTBins.At(i) < MaximumValueToUseRightSideFit)))
                leg.AddEntry(Kaon, "Kaon", "l");
            //leg.AddEntry((TObject*)0,Form("#chi^{2}/NDF = %1.2f", TotalFit->GetChisquare()/TotalFit->GetNDF()), "l");
            leg.Draw("same");
            //HadronContamination->BuildLegend();
            
            HadronContamination->SetTickx();
            HadronContamination->SetTicky();
            
            HadronContamination->Print(Form("Contamination_Cent_%d_pT%d_Config%d.pdf",fCentralityIndex, i,fConfigIndex));
            //HadronContamination->SaveAs(Form("Contamination_Cent_%d_pT%d_Config%d.root",fCentralityIndex, i,fConfigIndex));
        }
        
        
        
        Float_t LowerCut[7] = {-1};
        Float_t UpperCut = 3;
        Double_t ContaminationInSample = 0;
        
        for (Int_t iLowCut = 0; iLowCut < 1; iLowCut++)
        {
            Double_t ElectronIntegral = elec->Integral(LowerCut[iLowCut],UpperCut);
            Double_t PionIntegral = pion->Integral(LowerCut[iLowCut],UpperCut);
            Double_t KaonIntegral= Kaon->Integral(LowerCut[iLowCut],UpperCut);
            
            //KaonIntegral = 0;
            Double_t Contamination = (PionIntegral+ KaonIntegral)/(PionIntegral + ElectronIntegral + KaonIntegral);
            Double_t ErrorPionKaon = TMath::Sqrt(PionIntegral + KaonIntegral);
            Double_t DemError = TMath::Sqrt(PionIntegral + KaonIntegral + ElectronIntegral);
            Double_t ContaminationError = Contamination * TMath::Sqrt(pow(ErrorPionKaon/(PionIntegral + KaonIntegral),2) + pow(DemError/(PionIntegral + KaonIntegral+ElectronIntegral),2));
            ContaminationInSample = Contamination;
            
            ContamiantionUsingCuts[iLowCut]->SetBinContent(i+1, Contamination);
            ContamiantionUsingCuts[iLowCut]->SetBinError(i+1, ContaminationError);
            
            /// if (iLowCut == -0.5)
            
        }
        
        //calculate error associated, for now assume that
        
        
        //REMOVE THIS FOR PRODUCTION
        fHadronContamination->SetBinContent(i+1,1-ContaminationInSample);
        
        
        //fHadronContamination->SetBinContent(i+1,1);
        fHadronContamination->SetBinError(i+1,(1)/100000000000000000000.);
        
        
        fTPCNSigmaCenter->SetBinContent(i+1, TotalFit->GetParameter(1));
        fTPCNSigmaCenter->SetBinError(i+1, TotalFit->GetParError(1));
        
        fTPCNSigmaSTD->SetBinContent(i+1, TotalFit->GetParameter(2));
        fTPCNSigmaSTD->SetBinError(i+1, TotalFit->GetParError(2));
        
        
        /*
         //Double_t BinCounting = TPCNsigma[i]->Integrate(-0.5,3.0);
         printf("Bin: %d, Elec: %1.3f Pion = %1.3f; Total = %1.3f;Cont: = %1.3f \n", i+1, ElectronIntegral, PionIntegral, ElectronIntegral+PionIntegral, Contamination);
         */
    }
    
    /*
    ContaminationPlots->cd();
    Int_t iLowCut = 0;
    
    ContamiantionUsingCuts[iLowCut]->SetMarkerColor(kRed);
    ContamiantionUsingCuts[iLowCut]->GetYaxis()->SetRangeUser(10E-3,10);
    ContamiantionUsingCuts[iLowCut]->SetLineColor(kRed);
    ContamiantionUsingCuts[iLowCut]->SetMarkerSize(0.5);
    ContamiantionUsingCuts[iLowCut]->SetMarkerStyle(kFullCircle);
    ContamiantionUsingCuts[iLowCut]->Draw("same");
    
    
    
     Int_t Colors[7] = {TColor::GetColor("#e41a1c"), TColor::GetColor("#377eb8"), TColor::GetColor("#4daf4a"), TColor::GetColor("#984ea3"), TColor::GetColor("#ff7f00"), TColor::GetColor("#ffff33"), TColor::GetColor("#a65628")};
     Float_t LowerCut[7] = {-3,-2.5,-2,-1.5,-1.0,-0.5,0.};
     
     for (Int_t iLowCut = 0; iLowCut < 7; iLowCut++) {
     ContamiantionUsingCuts[iLowCut]->SetMarkerColor(Colors[iLowCut]);
     ContamiantionUsingCuts[iLowCut]->GetYaxis()->SetRangeUser(10E-3,10);
     ContamiantionUsingCuts[iLowCut]->SetLineColor(Colors[iLowCut]);
     ContamiantionUsingCuts[iLowCut]->SetTitle(Form("%1.1f < n#sigma TPC < 3.0 ",LowerCut[iLowCut]));
     ContamiantionUsingCuts[iLowCut]->SetMarkerSize(0.5);
     ContamiantionUsingCuts[iLowCut]->SetMarkerStyle(kFullCircle);
     ContamiantionUsingCuts[iLowCut]->Draw("same");
     }
     
     
     ContaminationPlots->SetLogy();
     ContaminationPlots->BuildLegend(0.6748166,0.1442308,0.8777506,0.5048077);
     */
    
    //HadronContamination->cd(fpTBins.GetSize());
    fHadronContamination->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // fHadronContamination->Draw();
    return kTRUE;
    
}
//NEEDS MC

Bool_t AliHFehpPbTool::CalculateTaggingEfficiency()
{
    printf("========================\n");
    printf("Processing pT Distributions\n");
    printf("========================\n");
    
    TH1F* ULS = (TH1F*) fInputList->FindObject("fEtaCutElectronBKULSMainSources");
    TH1F* LS = (TH1F*) fInputList->FindObject("fEtaCutElectronBKLSMainSources");
    
    //TH1F* ULS = (TH1F*) fInputList->FindObject("fPtElec_ULS");
    //TH1F* LS = (TH1F*) fInputList->FindObject("fPtElec_LS");
    
    
    TH1F *Total = (TH1F*) fInputList->FindObject("fEtaCutElectronBKNoTag");
    
    ULS->Add(LS,-1);
    
    ULS->Rebin(2);Total->Rebin(2);
    ULS = (TH1F*) ULS->Rebin(fpTBins.GetSize()-1, "BackgroundForTaggingEff", fpTBins.GetArray());
    Total = (TH1F*) Total->Rebin(fpTBins.GetSize()-1, "TotalBackgroundForTaggingEff", fpTBins.GetArray());
    
    ULS->Divide(ULS,Total,1,1,"B");
    
    
    fEffTagging = (TH1F*) ULS->Clone("fEffTagging");
    TFile *FileEff = new TFile("Eff_new_method_both_4_compare_me.root", "RECREATE");
    FileEff->cd();
    fEffTagging->Write("Eff");
    //fEffTagging->GetSumw2()->Set(0);
    //fEffTagging->Draw();
    
    return kTRUE;
    
}

void AliHFehpPbTool::Print2DCorrelation()
{
    for (int i = 0; i < fpTBins.GetSize() - 1 ; i++)
    {
        TCanvas print2Dplot("localprint","localprint",1200,300);
        
        print2Dplot.Divide(3,1);
        print2Dplot.cd(1);
        FormatTH2Titles(fHFEhSame[i]);
        fHFEhSame[i]->GetZaxis()->SetTitle("#it{N}^{ assoc}_{same}(#Delta#eta,#Delta#varphi)");
        fHFEhSame[i]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        fHFEhSame[i]->GetXaxis()->CenterTitle(true);
        fHFEhSame[i]->GetYaxis()->CenterTitle(true);

        
        fHFEhSame[i]->Draw("surf1");
        print2Dplot.cd(2);
        FormatTH2Titles(fHFEhMixed[i]);
        
        Double_t NormalizationMixed = (fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.),fHFEhMixed[i]->GetYaxis()->FindBin(0.)) + fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.)-1,fHFEhMixed[i]->GetYaxis()->FindBin(0.)) + fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.),fHFEhMixed[i]->GetYaxis()->FindBin(0.) -1) + fHFEhMixed[i]->GetBinContent(fHFEhMixed[i]->GetXaxis()->FindBin(0.)-1,fHFEhMixed[i]->GetYaxis()->FindBin(0.)-1))/4. ;
        
        
        //Finite Bin Correction
        NormalizationMixed = NormalizationMixed / (1 - 1./1.6 * fHFEhMixed[i]->GetYaxis()->GetBinCenter( fHFEhMixed[i]->GetYaxis()->FindBin(0.) ) );
        
        
        fHFEhMixed[i]->Scale(1./NormalizationMixed);

        fHFEhMixed[i]->Draw("surf1");
        fHFEhMixed[i]->GetZaxis()->SetTitle("#it{N}^{assoc}_{mixed}(#Delta#eta,#Delta#varphi)/#it{N}(0,0)");
        fHFEhMixed[i]->GetXaxis()->SetTitle("#Delta#varphi(rad)");
        fHFEhMixed[i]->GetXaxis()->CenterTitle(true);
        fHFEhMixed[i]->GetYaxis()->CenterTitle(true);

        print2Dplot.cd(3);
        FormatTH2Titles(fHFEh[i]);
        fHFEh[i]->Draw("surf1");
        fHFEh[i]->GetZaxis()->SetTitle("N_{eh}/#varepsilon_{tracking}");
        fHFEh[i]->GetXaxis()->SetTitle("#Delta#varphi(rad)");
        
        print2Dplot.Print(Form("HfeEffCorrected_%d_%d_%d.pdf",fConfigIndex,fCentralityIndex,i));
        
        
        TCanvas print2Dsame("localprintsame","localprint",600,450);
        gStyle->SetOptTitle(0);
        fHFEhSame[i]->Draw("surf1");
        
        TLatex* tex = new TLatex(0.02,0.08,"ALICE Preliminary");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.02,0.03,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.74,0.90, Form("%d < #it{p}_{T}^{ e} < %d GeV/#it{c}", (Int_t)fpTBins.At(i),(Int_t)fpTBins.At(i+1)));
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.7,0.84," 0.3 < #it{p}_{T}^{ assoc} < 2 GeV/#it{c}");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.68,0.96," 0-20% ZNA class");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.05,0.94,"(c,b) #rightarrow e - charged particles correlation");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.04728132);
        tex->SetLineColor(2);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.05,0.87,"Same Event");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineColor(2);
        tex->SetLineWidth(2);
        tex->Draw();
        
        print2Dsame.Print(Form("Same_%d_%d_%d.pdf",fConfigIndex,fCentralityIndex,i));
        print2Dsame.SaveAs(Form("Print2DSame_%d_%d_%d.root",fConfigIndex,fCentralityIndex,i));
        
        
        
        
        
        TCanvas print2Dmixed("localprintmixed","localprint",600,450);
        gStyle->SetOptTitle(0);
        fHFEhMixed[i]->Draw("surf1");
        
        tex = new TLatex(0.02,0.08,"ALICE Preliminary");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.02,0.03,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.74,0.90, Form("%d < #it{p}_{T}^{ e} < %d GeV/#it{c}", (Int_t)fpTBins.At(i),(Int_t)fpTBins.At(i+1)));
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.7,0.84," 0.3 < #it{p}_{T}^{ assoc} < 2 GeV/#it{c}");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.68,0.96," 0-20% ZNA class");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.05,0.94,"(c,b) #rightarrow e - charged particles correlation");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.04728132);
        tex->SetLineColor(2);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.05,0.87,"Mixed Event");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineColor(2);
        tex->SetLineWidth(2);
        tex->Draw();
        
        print2Dmixed.Print(Form("Mixed_%d_%d_%d.pdf",fConfigIndex,fCentralityIndex,i));
        print2Dsame.SaveAs(Form("Mixed_%d_%d_%d.root",fConfigIndex,fCentralityIndex,i));

        
        
        
        fHFEhMixed[i]->Scale(NormalizationMixed);
    }
    
    for (int i = 0; i < fpTBinsResults.GetSize() - 1 ; i++)
    {
        TCanvas *merged = new TCanvas(Form("merged%d",i) ,"merged",400,300);
        merged->cd();
        FormatTH2Titles(fHFEhNormalized[i]);
        fHFEhNormalized[i]->Draw("surf1");
        
        fHFEhNormalized[i]->GetZaxis()->SetTitle("#it{N}^{ assoc}(#Delta#eta,#Delta#varphi)/#it{N}^{(c,b) #rightarrow e}");
        fHFEhNormalized[i]->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        fHFEhNormalized[i]->GetXaxis()->CenterTitle(true);
        fHFEhNormalized[i]->GetYaxis()->CenterTitle(true);
        
        
        TLatex* tex = new TLatex(0.02,0.08,"ALICE Preliminary");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.02,0.03,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.74,0.90, Form("%d < #it{p}_{T}^{ e} < %d GeV/#it{c}", (Int_t)fpTBinsResults.At(i),(Int_t)fpTBinsResults.At(i+1)));
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.7,0.84," 0.3 < #it{p}_{T}^{ assoc} < 2 GeV/#it{c}");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.68,0.96," 0-20% ZNA class");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.05,0.94,"(c,b) #rightarrow e - charged particles correlation");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.04728132);
        tex->SetLineColor(2);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.05,0.87,"#frac{1}{#it{N}^{(c,b) #rightarrow e}} #frac{Same Event}{Mixed Event}");
        tex->SetNDC();
        tex->SetTextAlign(12);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03782506);
        tex->SetLineColor(2);
        tex->SetLineWidth(2);
        tex->Draw();

        merged->Print(Form("Hfe_merged_%d_%d_%d.pdf",fConfigIndex,fCentralityIndex,i));
    }
    
    //fHFEhNormalized[pT]; //
}

void AliHFehpPbTool::FormatTH2Titles(TH2F* histo)
{
    histo->GetYaxis()->SetRangeUser(-1.6,1.6);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitleOffset(1.3);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetYaxis()->SetTitleOffset(1.3);
    histo->GetZaxis()->SetTitleSize(0.035);
    histo->GetZaxis()->SetTitleOffset(1.4);
}

Bool_t AliHFehpPbTool::CalculateTaggingEfficiencyW()
{
    
    printf("========================\n");
    printf("eff with w\n");
    printf("========================\n");
    
    TH1F* ULS = (TH1F*) fInputList->FindObject("fEtaCutElectronBKULSMainSources_WithMotherW");
    TH1F* LS = (TH1F*) fInputList->FindObject("fEtaCutElectronBKLSMainSources_WithMotherW");
    TH1F *Total = (TH1F*) fInputList->FindObject("fEtaCutElectronBKNoTag_WithMotherW");
    
    ULS->Add(LS,-1);
    
    ULS = (TH1F*) ULS->Rebin(fpTBins.GetSize()-1, "BackgroundForTaggingEff", fpTBins.GetArray());
    Total = (TH1F*) Total->Rebin(fpTBins.GetSize()-1, "TotalBackgroundForTaggingEff", fpTBins.GetArray());
    
    ULS->Divide(ULS,Total,1,1,"B");
    
    fEffTagging = (TH1F*) ULS->Clone("fEffTagging");
    
    //delete ULS;
    //delete Total;
    
    //TFile *FileCompare = new TFile("CompareYvonne.root","RECREATE");
    //FileCompare->cd();
    //fEffTagging->Write("MyEff");
    
    
    if (fEffTagging)
        return kTRUE;
    
    
    return kFALSE;
    
}

Bool_t AliHFehpPbTool::CalculateTaggingEfficiencyPureHijing()
{
    
    printf("========================\n");
    printf("eff Using only Hijing \n");
    printf("========================\n");
    
    TH1F* ULS = (TH1F*) fInputList->FindObject("fElectronBKGNoEnhULS");
    TH1F* LS = (TH1F*) fInputList->FindObject("fElectronBKGNoEnhLS");
    TH1F *Total = (TH1F*) fInputList->FindObject("fElectronBKGNoEnhTotalNumber");
    
    ULS->Add(LS,-1);
    
    ULS = (TH1F*) ULS->Rebin(fpTBins.GetSize()-1, "BackgroundForTaggingEff", fpTBins.GetArray());
    Total = (TH1F*) Total->Rebin(fpTBins.GetSize()-1, "TotalBackgroundForTaggingEff", fpTBins.GetArray());
    
    ULS->Divide(ULS,Total,1,1,"B");
    
    fEffTagging = (TH1F*) ULS->Clone("fEffTagging");
    
    if (fEffTagging)
        return kTRUE;
    
    
    return kFALSE;
    
}

Bool_t AliHFehpPbTool::CalculateTaggingEfficiencyWToData()
{
    
    printf("========================\n");
    printf("Eff using W to Data \n");
    printf("========================\n");
    
    if (!fInputList)
        printf("No input list");
    else
        ("list ok");
    
    TH1F* ULS = (TH1F*) fInputList->FindObject("fElectronBKGWToDataULS");
    TH1F* LS = (TH1F*) fInputList->FindObject("fElectronBKGWToDataLS");
    TH1F *Total = (TH1F*) fInputList->FindObject("fElectronBKGWToDataTotal");
    
    if (!ULS || !LS || !Total)
        printf("Error reading histograms for efficiency calculation\n");
    else
        printf("reading ok");
    
    
    ULS->Add(LS,-1);
    
    ULS = (TH1F*) ULS->Rebin(fpTBins.GetSize()-1, "BackgroundForTaggingEff", fpTBins.GetArray());
    Total = (TH1F*) Total->Rebin(fpTBins.GetSize()-1, "TotalBackgroundForTaggingEff", fpTBins.GetArray());
    
    ULS->Divide(ULS,Total,1,1,"B");
    
    fEffTagging = (TH1F*) ULS->Clone("fEffTagging");
    
    if (fEffTagging)
        return kTRUE;
    
    
    printf("========================\n");
    printf("Finished Eff using W to Data \n");
    printf("========================\n");
    
    
    return kFALSE;
    
    
    
}


Bool_t AliHFehpPbTool::CalculateTaggingEfficiencyWToDataNoEnh()
{
    
    printf("========================\n");
    printf("Eff using W to Data \n");
    printf("========================\n");
    
    if (!fInputList)
        printf("No input list");
    else
        ("list ok");
    
    TH1F* ULS = (TH1F*) fInputList->FindObject("fElectronBKGNoEnhULS_WithW");
    TH1F* LS = (TH1F*) fInputList->FindObject("fElectronBKGNoEnhLS_WithW");
    TH1F *Total = (TH1F*) fInputList->FindObject("fElectronBKGNoEnhTotalNumber_WithW");
    
    if (!ULS || !LS || !Total)
        printf("Error reading histograms for efficiency calculation\n");
    else
        printf("reading ok");
    
    
    ULS->Add(LS,-1);
    
    ULS = (TH1F*) ULS->Rebin(fpTBins.GetSize()-1, "BackgroundForTaggingEff", fpTBins.GetArray());
    Total = (TH1F*) Total->Rebin(fpTBins.GetSize()-1, "TotalBackgroundForTaggingEff", fpTBins.GetArray());
    
    ULS->Divide(ULS,Total,1,1,"B");
    
    fEffTagging = (TH1F*) ULS->Clone("fEffTagging");
    
    if (fEffTagging)
        return kTRUE;
    
    
    printf("========================\n");
    printf("Finished Eff using W to Data \n");
    printf("========================\n");
    
    
    return kFALSE;
    
    
    
}




Bool_t AliHFehpPbTool::ReadTaggingEfficiencyFromFile()
{
    TFile *FileEff = new TFile("Eff_new_method.root");
    fEffTagging = (TH1F*) FileEff->Get("Eff");
    //fEffTagging->Draw();
    
}

Double_t AliHFehpPbTool::GetNEvents()
{
    TH1F *Events = (TH1F*) fInputList->FindObject("fNevent");
    return Events->GetBinContent(1);
    
    
}

Bool_t AliHFehpPbTool::CalculateHFeEfficiency()
{
    //TCanvas *HFeEfficiency = new TCanvas("HFeEff", "HFeEff",400,300);
    
    TH1F *Reco = (TH1F*) fInputList->FindObject("fEtaCutElectronRecoHFEMC");
    Reco = (TH1F*) Reco->Rebin(fpTBins.GetSize()-1, "HFFindinfEff", fpTBins.GetArray());
    //Reco->Rebin(5);
    TH1F *MC = (TH1F*) fInputList->FindObject("fEtaCutElectronGeneratedSignalPtEtaZvtx");
    MC = (TH1F*) MC->Rebin(fpTBins.GetSize()-1, "HFFindinfEff", fpTBins.GetArray());
    Reco->Divide(Reco,MC,1,1,"B");
    
    //Reco->SetTitle("HeavyFlavor electron efficiency");
    //Reco->Draw();
    fEffHFe = (TH1F*) Reco->Clone("fEffHFe");
    
    if (fEffHFe)
        return kTRUE;
    
    return kFALSE;
}

//
void AliHFehpPbTool::CalculateFlow1D(AliHFehpPbTool* Reference, Bool_t FitConstant)
{
    TH1F **JetReference = new TH1F *[fpTBinsResults.GetSize()];
    fFlowHistograms = new TH1F *[fpTBinsResults.GetSize()];
    TF1 **FlowFunction = new TF1 *[fpTBinsResults.GetSize()];
    
    for (Int_t i = 0; i < fpTBinsResults.GetSize() - 1 ; i++)
    {
        if (FitConstant)
            FlowFunction[i] = new TF1(Form("fit%d",i),"[0]",-0.5*TMath::Pi(),1.5*TMath::Pi());
        else
        {
            FlowFunction[i] = new TF1(Form("fit%d",i),"[0]*(1 + 2 * [1] * TMath::Cos(x) + 2 * [2] * TMath::Cos(2*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            FlowFunction[i]->SetParameters(10.,-0.05, 0.1);
            //FlowFunction[i]->FixParameter(1,0);
        }
        
        JetReference[i] = Reference->GetHFeh1DSub(i);
        
        
        fFlowHistograms[i] = (TH1F*) fHFEhNormalized1D[i]->Clone(Form("FlowHistograms%d",i));
        fFlowHistograms[i]->Add(JetReference[i],-1);
       // fFlowHistograms[i]->Fit(FlowFunction[i],"0");
        
        
    }
    
    
}

/*
 void AliHFehpPbTool::CalculateFlow1FromFitTemplate(AliHFehpPbTool* Reference)
{
    TH1F **JetReference = new TH1F *[fpTBinsResults.GetSize()];
    fFlowHistograms = new TH1F *[fpTBinsResults.GetSize()];
    TF1 **FlowFunction = new TF1 *[fpTBinsResults.GetSize()];
    
    for (Int_t i = 0; i < fpTBinsResults.GetSize() - 1 ; i++)
    {
        if (FitConstant)
            FlowFunction[i] = new TF1(Form("fit%d",i),"[0]",-0.5*TMath::Pi(),1.5*TMath::Pi());
        else
        {
            FlowFunction[i] = new TF1(Form("fit%d",i),"[0]*(1 + 2 * [1] * TMath::Cos(x) + 2 * [2] * TMath::Cos(2*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            FlowFunction[i]->SetParameters(10.,-0.05, 0.1);
            //FlowFunction[i]->FixParameter(1,0);
        }
        
        JetReference[i] = Reference->GetHFeh1DSub(i);
        
        
        fFlowHistograms[i] = (TH1F*) fHFEhNormalized1D[i]->Clone(Form("FlowHistograms%d",i));
        fFlowHistograms[i]->Add(JetReference[i],-1);
        // fFlowHistograms[i]->Fit(FlowFunction[i],"0");
        
        
    }
    
}
*/

TF1* AliHFehpPbTool::DoTemplateJetFit(TH1F* Correlation)
{
    TF1 *FitGaussians = new  TF1("TemplateFit","gaus(0)+gaus(3)+[6]",-2,5);
    FitGaussians->SetParameters(1.1,0,0.3,0.6,TMath::Pi(),3,Correlation->GetBinContent(10));
    FitGaussians->FixParameter(1,0);
    FitGaussians->FixParameter(4,TMath::Pi());
    Correlation->Fit(FitGaussians,"I0EM");
    
    return FitGaussians;
    
}

TF1* AliHFehpPbTool::DoTemplatev2Fit(TH1F* Correlation, TF1* JetFit)
{
    TF1 *FitGaussians = new  TF1("TemplateFit","gaus(0)+gaus(3)+[6]*(1+2*[7]*cos(x)+2*[8]*cos(2*x))",-2,5);
    FitGaussians->SetParameters(1,0,1,1,TMath::Pi(),2,Correlation->GetBinContent(10),-0.01,0.027);
    FitGaussians->FixParameter(0,JetFit->GetParameter(0));
    FitGaussians->FixParameter(1,JetFit->GetParameter(1));
    FitGaussians->FixParameter(2,JetFit->GetParameter(2));
    FitGaussians->FixParameter(3,JetFit->GetParameter(3));
    FitGaussians->FixParameter(4,JetFit->GetParameter(4));
    FitGaussians->FixParameter(5,JetFit->GetParameter(5));
    
    Correlation->Fit(FitGaussians,"QI0EM");
    
    return FitGaussians;
}


void AliHFehpPbTool::SystematicFromv2Fit(TH1F* Correlation,TF1* JetFit)
{
   //Obtain values that will be varied: Sigma and gaussian normalization
    
    Double_t p0NS(JetFit->GetParameter(0)), p2NS(JetFit->GetParameter(2)),p0AS(JetFit->GetParameter(3)),p2AS(JetFit->GetParameter(5));
    Double_t p0NSError(JetFit->GetParError(0)), p2NSError(JetFit->GetParError(2)),p0ASError(JetFit->GetParError(3)),p2ASError(JetFit->GetParError(5));
    
    TF1* FirstResult = DoTemplatev2Fit(Correlation, JetFit);
    Double_t v2FromFit(FirstResult->GetParameter(8));
    printf("v2delta = %1.6f +- %1.6f\n",v2FromFit, FirstResult->GetParError(8));
    
    //For each of them, varry randomly with gaussian distribution
    
    TH1F *v2Results = new TH1F("v2FitVariations", "v2FitVariations;v2;Counts;", 100,0,0.01);
    TRandom3 Generator;
    for (int i = 0; i < 10000; i++)
    {
        Double_t p0NSRand(Generator.Gaus(p0NS,p0NSError)),p2NSRand(Generator.Gaus(p2NS,p2NSError)),p0ASRand(Generator.Gaus(p0AS,p0ASError)),p2ASRand(Generator.Gaus(p2AS,p2ASError));
        JetFit->SetParameters(p0NSRand, 0, p2NSRand, p0ASRand, TMath::Pi(), p2ASRand);
        
        JetFit->SetParameters(p0NSRand, 0, p2NSRand, p0ASRand, TMath::Pi(), p2ASRand);
        
        TF1* RandonResults = DoTemplatev2Fit(Correlation, JetFit);
        //printf("v2delta = %1.6f, Chi2Red = %1.2f \n",RandonResults->GetParameter(8),RandonResults->GetChisquare()/RandonResults->GetNDF());
        
        if (RandonResults->GetChisquare()/RandonResults->GetNDF() < 2.0)
            v2Results->Fill(RandonResults->GetParameter(8));
    }
    
    v2Results->Draw();
}


void AliHFehpPbTool::DrawTemplatev2Fit(TH1F* Correlation)
{
    //TCanvas *Canvas = new TCanvas
    TF1 *Fit = Correlation->GetFunction("TemplateFit");
    Fit->SetParName(6,"P");
    Fit->SetParName(7,"V_{1#Delta}");
    Fit->SetParName(8,"V_{2#Delta}");


    Fit->SetLineColor(TColor::GetColor("#2196F3"));
    TF1 *Gaussians = new  TF1("Jet Component","gaus(0)+gaus(3)+[6]",-2,5);
    Gaussians->SetParameters(Fit->GetParameter(0),Fit->GetParameter(1), Fit->GetParameter(2), Fit->GetParameter(3), Fit->GetParameter(4), Fit->GetParameter(5),Fit->GetParameter(6));
    TF1 *V1 = new  TF1("TemplateFit","[0]*(1 + 2*[1]*cos(x))",-2,5);
    Gaussians->SetLineColor(TColor::GetColor("#673AB7"));
    Gaussians->SetLineStyle(4);
    V1->SetParameters(Fit->GetParameter(6),Fit->GetParameter(7));
    V1->SetLineColor(TColor::GetColor("#6A1B9A"));
    V1->SetLineStyle(3);
    
    TF1 *V2 = new  TF1("TemplateFit","[0]*(1 + 2*[1]*cos(2*x))",-2,5);
    V2->SetParameters(Fit->GetParameter(6),Fit->GetParameter(8));
    V2->SetLineColor(TColor::GetColor("#B71C1C"));
    V2->SetLineStyle(2);
    
    Correlation->SetLineColor(kBlack);
    Correlation->SetMarkerSize(0.5);
    Correlation->SetMarkerStyle(21);
    Correlation->SetMarkerColor(kBlack);
    Correlation->Draw();
    Gaussians->Draw("same");
    V1->Draw("same");
    V2->Draw("same");
    Fit->Draw("same");
    
}

void AliHFehpPbTool::DrawTemplateJetFit(TH1F* Correlation)
{
    //TCanvas *Canvas = new TCanvas
    TF1 *Fit = Correlation->GetFunction("TemplateFit");
    Fit->SetParName(0,"NS norm.");
    Fit->SetParName(2,"NS #sigma");
    Fit->SetParName(3,"AS norm.");
    Fit->SetParName(5,"AS #sigma");
    Fit->SetParName(6,"Constant");

    Fit->SetLineColor(TColor::GetColor("#2196F3"));
    TF1 *Gaussians = new  TF1("Jet Component","gaus(0)+gaus(3)+[6]",-2,5);
    Gaussians->SetParameters(Fit->GetParameter(0),Fit->GetParameter(1), Fit->GetParameter(2), Fit->GetParameter(3), Fit->GetParameter(4), Fit->GetParameter(5),Fit->GetParameter(6));
    Gaussians->SetLineColor(TColor::GetColor("#673AB7"));
    
    
    //Gaussians->SetLineStyle(4);
    /*
    V1->SetParameters(Fit->GetParameter(6),Fit->GetParameter(7));
    V1->SetLineColor(TColor::GetColor("#6A1B9A"));
    V1->SetLineStyle(3);
    
    TF1 *V2 = new  TF1("TemplateFit","[0]*(1 + 2*[1]*cos(2*x))",-2,5);
    V2->SetParameters(Fit->GetParameter(6),Fit->GetParameter(8));
    V2->SetLineColor(TColor::GetColor("#B71C1C"));
    V2->SetLineStyle(2);
     */
    
    Correlation->GetYaxis()->SetTitle("C_{LM}");
    
    Correlation->SetLineColor(kBlack);
    Correlation->SetMarkerSize(0.5);
    Correlation->SetMarkerStyle(21);
    Correlation->SetMarkerColor(kBlack);
    Correlation->Draw();
    Gaussians->Draw("same");
    //V1->Draw("same");
    //V2->Draw("same");
    ///Fit->Draw("same");
    
}






Bool_t AliHFehpPbTool::NormalizeHFeCorrrelation()
{
    printf("Normalizing by the spectra \n");
    TH1F *HFepTRebin = (TH1F*)fHFepT->Rebin(fpTBinsResults.GetSize()-1, "HFepTRebin", fpTBinsResults.GetArray());
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        fHFEhNormalized[i]->Scale(1./HFepTRebin->GetBinContent(i+1));
        fHFEhSameNormalized[i]->Scale(1./HFepTRebin->GetBinContent(i+1));
        
        //Mixed requires a bit different workaround...
        Double_t NormalizationMixed = (fHFEhMixedNormalized[i]->GetBinContent(fHFEhMixedNormalized[i]->GetXaxis()->FindBin(0.),fHFEhMixedNormalized[i]->GetYaxis()->FindBin(0.)) + fHFEhMixedNormalized[i]->GetBinContent(fHFEhMixedNormalized[i]->GetXaxis()->FindBin(0.)-1,fHFEhMixedNormalized[i]->GetYaxis()->FindBin(0.)) + fHFEhMixedNormalized[i]->GetBinContent(fHFEhMixedNormalized[i]->GetXaxis()->FindBin(0.),fHFEhMixedNormalized[i]->GetYaxis()->FindBin(0.) -1) + fHFEhMixedNormalized[i]->GetBinContent(fHFEhMixedNormalized[i]->GetXaxis()->FindBin(0.)-1,fHFEhMixedNormalized[i]->GetYaxis()->FindBin(0.)-1))/4. ;
        
        //Finite Bin Correction
        NormalizationMixed = NormalizationMixed / (1 - 1./1.6 * fHFEhMixedNormalized[i]->GetYaxis()->GetBinCenter( fHFEhMixedNormalized[i]->GetYaxis()->FindBin(0.) ) );
        
        fHFEhMixedNormalized[i]->Scale(1./NormalizationMixed);
        
    }
}

Bool_t AliHFehpPbTool::NormalizeHFeCorrrelationMC()
{
    TH1F *pTMC = (TH1F*) fInputList->FindObject("fEtaCutElectronGeneratedSignalPtEtaZvtx");
    if (!pTMC) printf("error reading histogram for normalizaton");
    pTMC = (TH1F*) pTMC->Rebin(fpTBinsResults.GetSize()-1, "ptMC", fpTBinsResults.GetArray());
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
        fHFEhMCNormalized[i]->Scale(1./pTMC->GetBinContent(i+1));
}

Bool_t AliHFehpPbTool::ProjectMCTo1D()
{
    fHFEhMCNormalized1D = new TH1F *[fpTBinsResults.GetSize() -1];
    
    
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        fHFEhMCNormalized1D[i] = (TH1F*) fHFEhMCNormalized[i]->ProjectionX(Form("fHFEhMCNormalized1D%d",i),1, fHFEhMCNormalized[i]->GetNbinsY());
        fHFEhMCNormalized1D[i]->Scale(1.,"width");
        
        
        /*
         
         Int_t MinBin[3] = {999,999,999};
         Double_t MinValue[3] = {999,999,999};
         
         for (Int_t j = 1; j <=  fHFEhMCNormalized1D[i]->GetNbinsX() ; j++)
         {
         if (fHFEhMCNormalized1D[i]->GetBinContent(j) <= MinValue[0])
         {
         MinBin[0] = j;
         MinValue[0] = fHFEhMCNormalized1D[i]->GetBinContent(j);
         }
         
         
         }
         
         for (Int_t j = 1; j <=  fHFEhMCNormalized1D[i]->GetNbinsX() ; j++)
         {
         if (j == MinBin[0])
         continue;
         
         if (fHFEhMCNormalized1D[i]->GetBinContent(j) <= MinValue[1])
         {
         MinBin[1] = j;
         MinValue[1] =fHFEhMCNormalized1D[i]->GetBinContent(j);
         }
         }
         
         for (Int_t j = 1; j <=  fHFEhMCNormalized1D[i]->GetNbinsX() ; j++)
         {
         if (j == MinBin[0] || j == MinBin[1] )
         continue;
         
         if (fHFEhMCNormalized1D[i]->GetBinContent(j) <= MinValue[2])
         {
         MinBin[2] = j;
         MinValue[2] =fHFEhMCNormalized1D[i]->GetBinContent(j);
         }
         }
         
         Double_t Pedestal = (MinValue[0] + MinValue[1] + MinValue[2])/3.;
         
         TF1 *PedestalFunction = new TF1("PedestalFunction",Form("%1.4f",Pedestal), -3,7);
         
         fHFEhMCNormalized1D[i]->Add(PedestalFunction,-1);
         
         */
        
    }
    
    
    
}

Bool_t AliHFehpPbTool::ProjectTo1D()
{
    Int_t FirstBin = fHFEhNormalized[0]->GetYaxis()->FindBin(fMinDeltaEta+0.01);
    Int_t LastBin = fHFEhNormalized[0]->GetYaxis()->FindBin(fMaxDeltaEta-0.01);

    printf("Projecting to 1D\n");
    printf("From %1.2f to %1.2f \n",fHFEhNormalized[0]->GetYaxis()->GetBinLowEdge(FirstBin),fHFEhNormalized[0]->GetYaxis()->GetBinLowEdge(LastBin)+ fHFEhNormalized[0]->GetYaxis()->GetBinWidth(LastBin));
    fHFEhNormalized1D = new TH1F *[fpTBinsResults.GetSize() -1];
    fHFEhNormSub1D = new TH1F *[fpTBinsResults.GetSize() -1];
    
    fBaseline = new TH1F("fBaseline", "Baseline", fpTBinsResults.GetSize() -1, fpTBinsResults.GetArray());
    
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        fHFEhNormalized1D[i] = (TH1F*) fHFEhNormalized[i]->ProjectionX(Form("fHFEhNormalized1D%d",i),FirstBin, LastBin);
        fHFEhNormalized1D[i]->Scale(1.,"width");
        fHFEhNormalized1D[i]->Scale(1.0/(fMaxDeltaEta-fMinDeltaEta)); // Normalize by eta size
        fHFEhNormSub1D[i] = (TH1F*)  fHFEhNormalized1D[i]->Clone(Form("fHFEhNormSub1D%d",i));
        
        Int_t MinBin[3] = {1,2,3};
        Double_t MinValue[3] = {fHFEhNormalized1D[i]->GetBinContent(1),fHFEhNormalized1D[i]->GetBinContent(2),fHFEhNormalized1D[i]->GetBinContent(3)};
        
        for (Int_t j = 4; j <=  fHFEhNormalized1D[i]->GetNbinsX() ; j++)
        {
            if (fHFEhNormalized1D[i]->GetBinContent(j) <= MinValue[0])
            {
                MinBin[0] = j;
                MinValue[0] =fHFEhNormalized1D[i]->GetBinContent(j);
            }
            
        }
        
        for (Int_t j = 4; j <=  fHFEhNormalized1D[i]->GetNbinsX() ; j++)
        {
            if (j == MinBin[0])
            continue;
            
            if (fHFEhNormalized1D[i]->GetBinContent(j) <= MinValue[1])
            {
                MinBin[1] = j;
                MinValue[1] =fHFEhNormalized1D[i]->GetBinContent(j);
            }
        }
        
        for (Int_t j = 4; j <=  fHFEhNormalized1D[i]->GetNbinsX() ; j++)
        {
            if (j == MinBin[0] || j == MinBin[1] )
            continue;
            
            if (fHFEhNormalized1D[i]->GetBinContent(j) <= MinValue[2])
            {
                MinBin[2] = j;
                MinValue[2] =fHFEhNormalized1D[i]->GetBinContent(j);
            }
        }
        
        Double_t MinValueError[3] = {fHFEhNormalized1D[i]->GetBinError(MinBin[0]),fHFEhNormalized1D[i]->GetBinError(MinBin[1]),fHFEhNormalized1D[i]->GetBinError(MinBin[2])};
        
        
        Double_t Pedestal = (MinValue[0]*1./(MinValueError[0]*MinValueError[0]) + MinValue[1]*1./(MinValueError[1]*MinValueError[1]) + MinValue[2]*1./(MinValueError[2]*MinValueError[2]));
        
        Pedestal = Pedestal/(1./(MinValueError[0]*MinValueError[0]) + 1./(MinValueError[1]*MinValueError[1]) + 1./(MinValueError[2]*MinValueError[2]) );
        
        Double_t PedestaError = 1./TMath::Sqrt(1./(MinValueError[0]*MinValueError[0]) + 1./(MinValueError[1]*MinValueError[1]) + 1./(MinValueError[2]*MinValueError[2]));
        printf("\n\n ================== \n pedestal error = %1.7f\n", PedestaError);
        
        TF1 *PedestalFunction = new TF1("PedestalFunction",Form("%1.4f",Pedestal), -3,7);
        
        fHFEhNormSub1D[i]->Add(PedestalFunction,-1);
        
        fBaseline->SetBinContent(i+1,Pedestal);
        fBaseline->SetBinError(i+1,PedestaError);
        
        
    }
    
}


void AliHFehpPbTool::ProjectTo1DAtlas()
{
    Int_t FirstBin = fHFEhNormalized[0]->GetYaxis()->FindBin(fMinDeltaEta+0.01);
    Int_t LastBin = fHFEhNormalized[0]->GetYaxis()->FindBin(fMaxDeltaEta-0.01);
    
    printf("Projecting to 1D\n");
    printf("From %1.2f to %1.2f \n",fHFEhNormalized[0]->GetYaxis()->GetBinLowEdge(FirstBin),fHFEhNormalized[0]->GetYaxis()->GetBinLowEdge(LastBin)+ fHFEhNormalized[0]->GetYaxis()->GetBinWidth(LastBin));
    fHFEhNormalized1D = new TH1F *[fpTBinsResults.GetSize() -1];
    fHFEhNormSub1D = new TH1F *[fpTBinsResults.GetSize() -1];
    
    fBaseline = new TH1F("fBaseline", "Baseline", fpTBinsResults.GetSize() -1, fpTBinsResults.GetArray());
    
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        TH1F *Same = (TH1F*) fHFEhSameNormalized[i]->ProjectionX(Form("fHFEhSameNormalized1D%d",i),FirstBin, LastBin);
        TH1F *Mixed = (TH1F*) fHFEhMixedNormalized[i]->ProjectionX(Form("fHFEhMixedNormalized1D%d",i),FirstBin, LastBin);
        Mixed->Scale(1./6.);
        Same->Divide(Mixed);
        
        fHFEhNormalized1D[i] = Same;
        fHFEhNormalized1D[i]->Scale(1.,"width");
        fHFEhNormalized1D[i]->Scale(1.0/(fMaxDeltaEta-fMinDeltaEta)); // Normalize by eta size
        fHFEhNormSub1D[i] = (TH1F*)  fHFEhNormalized1D[i]->Clone(Form("fHFEhNormSub1D%d",i));
        
        Int_t MinBin[5] = {1,2,3};
        Double_t MinValue[3] = {fHFEhNormalized1D[i]->GetBinContent(1),fHFEhNormalized1D[i]->GetBinContent(2),fHFEhNormalized1D[i]->GetBinContent(3)};
        
        for (Int_t j = 4; j <=  fHFEhNormalized1D[i]->GetNbinsX() ; j++)
        {
            if (fHFEhNormalized1D[i]->GetBinContent(j) <= MinValue[0])
            {
                MinBin[0] = j;
                MinValue[0] =fHFEhNormalized1D[i]->GetBinContent(j);
            }
            
        }
        
        for (Int_t j = 4; j <=  fHFEhNormalized1D[i]->GetNbinsX() ; j++)
        {
            if (j == MinBin[0])
                continue;
            
            if (fHFEhNormalized1D[i]->GetBinContent(j) <= MinValue[1])
            {
                MinBin[1] = j;
                MinValue[1] =fHFEhNormalized1D[i]->GetBinContent(j);
            }
        }
        
        for (Int_t j = 4; j <=  fHFEhNormalized1D[i]->GetNbinsX() ; j++)
        {
            if (j == MinBin[0] || j == MinBin[1] )
                continue;
            
            if (fHFEhNormalized1D[i]->GetBinContent(j) <= MinValue[2])
            {
                MinBin[2] = j;
                MinValue[2] =fHFEhNormalized1D[i]->GetBinContent(j);
            }
        }
        
        Double_t MinValueError[3] = {fHFEhNormalized1D[i]->GetBinError(MinBin[0]),fHFEhNormalized1D[i]->GetBinError(MinBin[1]),fHFEhNormalized1D[i]->GetBinError(MinBin[2])};
        
        
        Double_t Pedestal = (MinValue[0]*1./(MinValueError[0]*MinValueError[0]) + MinValue[1]*1./(MinValueError[1]*MinValueError[1]) + MinValue[2]*1./(MinValueError[2]*MinValueError[2]));
        
        Pedestal = Pedestal/(1./(MinValueError[0]*MinValueError[0]) + 1./(MinValueError[1]*MinValueError[1]) + 1./(MinValueError[2]*MinValueError[2]) );
        
        Double_t PedestaError = 1./TMath::Sqrt(1./(MinValueError[0]*MinValueError[0]) + 1./(MinValueError[1]*MinValueError[1]) + 1./(MinValueError[2]*MinValueError[2]));
        printf("\n\n ================== \n pedestal error = %1.7f\n", PedestaError);
        
        TF1 *PedestalFunction = new TF1("PedestalFunction",Form("%1.4f",Pedestal), -3,7);
        
        fHFEhNormSub1D[i]->Add(PedestalFunction,-1);
        
        fBaseline->SetBinContent(i+1,Pedestal);
        fBaseline->SetBinError(i+1,PedestaError);
        
    }
    
}


Bool_t AliHFehpPbTool::SubtractPedestal(Double_t Pedestal, Int_t pT)
{
    TF1 *PedestalFunction = new TF1("PedestalFunction",Form("%1.4f",Pedestal), -3,7);
    delete fHFEhNormSub1D[pT];
    fHFEhNormSub1D[pT] = (TH1F*)  fHFEhNormalized1D[pT]->Clone(Form("fHFEhNormSub1D%d",pT));
    fHFEhNormSub1D[pT]->Add(PedestalFunction,-1);
    
}

Bool_t AliHFehpPbTool::CalculateYield(Bool_t Flow)
{
    TF1 **FitYield = new TF1 *[fpTBinsResults.GetSize() -1];
    
    fYieldAS = new TH1F("fYieldAS", "Yields AS",fpTBinsResults.GetSize() -1, fpTBinsResults.GetArray());
    fYieldNS = new TH1F("fYieldNS", "Yields NS",fpTBinsResults.GetSize() -1, fpTBinsResults.GetArray());
    
    for (Int_t pT = 0 ; pT < fpTBinsResults.GetSize() -1 ; pT++)
    {
        //TCanvas *Yield = new TCanvas(Form(Form(""),));
        //NS gaussian -> [0] = area, [1] = mean , [2] = sigma
        //AS gaussian -> [3] = area, [4] = mean , [5] = sigma
        //[6] = pedestal
        //FitYield[pT] = new TF1 (Form("FitYield%d",pT), "[0]*(1 + 2*[1]*TMath::Cos(2*x)) + gaus(2) + gaus (5)", -0.5*TMath::Pi(),1.5*TMath::Pi());
        FitYield[pT] = new TF1 (Form("FitYield%d",pT), "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2])) + [3]/TMath::Sqrt(2.*TMath::Pi())/[5]*TMath::Exp(-(x-[4])*(x-[4])/2./([5]*[5])) +[6] ",-0.5*TMath::Pi(),1.5*TMath::Pi());
        
        FitYield[pT]->SetParameters(2., 0. , 0.7 , 1., TMath::Pi(), 0.7, fHFEhNormalized1D[pT]->GetBinContent(fHFEhNormalized1D[pT]->GetNbinsX()/2));
        // if (Flow)
        //{
        //  TF1 *FlowFunction = fFlowHistograms[pT]->GetFunction(Form("fit%d",pT));
        //FitYield[pT]->FixParameter(0, FlowFunction->GetParameter(0));
        // FitYield[pT]->FixParameter(1, FlowFunction->GetParameter(1));
        //}
        //else
        //{
        //FitYield[pT]->FixParameter(0, Pedestal);
        //FitYield[pT]->FixParameter(1, 0);
        // }
        
        FitYield[pT]->FixParameter(1, 0); //Near side peak
        FitYield[pT]->FixParameter(4, TMath::Pi()); // Away side peak
        //FitYield[pT]->FixParameter(6, 0);
        
        FitYield[pT]->SetParLimits(5, 0.5,3);
        FitYield[pT]->SetParLimits(2, 0.3,3);
        
        
        //fHFEhNormSub1D[pT]->Fit(FitYield[pT]);
        fHFEhNormalized1D[pT]->Fit(FitYield[pT], "0+");
        
        double yieldNS = FitYield[pT]->GetParameter(0);
        double errorNS = FitYield[pT]->GetParError(0);
        
        double yieldAS =  FitYield[pT]->GetParameter(3);
        double errorAS =FitYield[pT]->GetParError(3);
        
        printf("Yield NS: %1.4f +-%1.4f, AN: %1.4f+-%1.4f \n", yieldNS,errorNS,yieldAS, errorAS);
        fYieldNS->SetBinContent(pT+1, yieldNS);
        fYieldNS->SetBinError(pT+1, errorNS);
        
        fYieldAS->SetBinContent(pT+1, yieldAS);
        fYieldAS->SetBinError(pT+1, errorAS);
        
        SubtractPedestal(FitYield[pT]->GetParameter(6),pT);
        
    }
    
    
    
}



Bool_t AliHFehpPbTool::CorrelationCT1D()
{
    TCanvas **Canvas = new TCanvas *[fpTBinsResults.GetSize() -1];
    TH1F **Ratio1D = new TH1F *[fpTBinsResults.GetSize() -1];
    
    gStyle->SetOptFit(1);
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        Canvas[i] = new TCanvas(Form("CorrelationCT1D%d",i),Form("CorrelationCT1D%d",i),400,300);
        Canvas[i]->cd();
        Ratio1D[i] = (TH1F*) fHFEhNormalized1D[i]->Clone(Form("RatioDataMC1D%d",i));
        Ratio1D[i]->Divide( Ratio1D[i], fHFEhMCNormalized1D[i]);
        Ratio1D[i]->Fit("pol0");
        Canvas[i]->SaveAs(Form("CorrelationCT1D%d",i));
        
    }
}

Bool_t AliHFehpPbTool::CalculateV22PC()
{
    fHFehProjectionForV2NonSub = new TH1F *[fpTBinsResults.GetSize() -1];
    
    for (int pT = 0; pT < fpTBinsResults.GetSize() -1; pT++)
    {
        fHFehProjectionForV2NonSub[pT] = (TH1F*)fHFEhNormalized1D[pT]->Clone(Form("fHFehProjectionForV2NonSub%d",pT));
        fHFehProjectionForV2NonSub[pT]->Reset();
        
        for (Int_t ix =  1 ; ix <= fHFEhNormalized[pT]->GetNbinsX(); ix++)
        {
            Double_t Mean = 0;
            Double_t Normalization =0;
            
            if (ix <= fHFEhNormalized[pT]->GetNbinsX()/2)
            {
                
                for (Int_t iy = fHFEhNormalized[pT]->GetYaxis()->FindBin(-1.55) ; iy <= fHFEhNormalized[pT]->GetYaxis()->FindBin(-0.85) ; iy++)
                {
                    Mean += fHFEhNormalized[pT]->GetBinContent(ix,iy)/( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
                    Normalization += 1./( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
                    //printf("iy = %d normalization = %f mean = %f \n",iy,Normalization,Mean);
                }
                
                for (Int_t iy = fHFEhNormalized[pT]->GetYaxis()->FindBin(0.85) ; iy <= fHFEhNormalized[pT]->GetYaxis()->FindBin(1.55) ; iy++)
                {
                    Mean += fHFEhNormalized[pT]->GetBinContent(ix,iy)/( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
                    Normalization += 1./( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
                    //printf("iy = %d normalization = %f mean = %f \n",iy,Normalization,Mean);
                }
                
                fHFehProjectionForV2NonSub[pT]->SetBinContent(ix, Mean/Normalization);
                fHFehProjectionForV2NonSub[pT]->SetBinError(ix, TMath::Sqrt(1./Normalization));
            }
            else
            {
                for (Int_t iy = fHFEhNormalized[pT]->GetYaxis()->FindBin(-1.55) ; iy <= fHFEhNormalized[pT]->GetYaxis()->FindBin(1.55) ; iy++)
                {
                    Mean += fHFEhNormalized[pT]->GetBinContent(ix,iy)/( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
                    Normalization += 1./( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
                    //printf("iy = %d normalization = %f mean = %f \n",iy,Normalization,Mean);
                }
                fHFehProjectionForV2NonSub[pT]->SetBinContent(ix, Mean/Normalization);
                fHFehProjectionForV2NonSub[pT]->SetBinError(ix, TMath::Sqrt(1./Normalization));
                
                
            }
            
            
            
            
        }
    }
    
    
}


void AliHFehpPbTool::CalculateMCWeight()
{
    Double_t bins[83] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.5, 9.0, 9.5, 10.0, 11, 12, 13, 14, 15, 16 , 18, 20, 25, 30, 35, 40, 50};
    Int_t Nbins =82;
    
    TFile *FileW = new TFile("BackgroundW.root", "RECREATE");
    
    TH1F *temp = (TH1F*) fInputList->FindObject("fPtMCpi0_NoMother");
    TH1F* PionEnh = (TH1F*) temp->Rebin(Nbins,"PionEnh",bins);
    PionEnh->SetTitle("#pi^{0} from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    PionEnh->Scale(1.,"width");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_NoMother");
    TH1F *EtaEnh = (TH1F*) temp->Rebin(Nbins,"EtaEnh",bins);
    EtaEnh->SetTitle("#eta from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    EtaEnh->Scale(1.,"width");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_PureHijing");
    TH1F* PionHijing = (TH1F*) temp->Rebin(Nbins,"PionHijing",bins);
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    PionHijing->Scale(1.,"width");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_PureHijing");
    TH1F *EtaHijing = (TH1F*) temp->Rebin(Nbins,"EtaHijing",bins);
    EtaHijing->Scale(1.,"width");
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    TCanvas *PionCanvas = new TCanvas("Pion","Pion",400,300);
    TCanvas *EtaCanvas = new TCanvas("Eta","Eta",400,300);
    
    TCanvas *PionCanvasW = new TCanvas("PionW","PionW",400,300);
    TCanvas *EtaCanvasW = new TCanvas("EtaW","EtaW",400,300);
    
    PionCanvas->cd();
    PionCanvas->SetLogy();
    PionHijing->SetMarkerSize(0.3);
    PionHijing->SetMarkerStyle(kOpenCircle);
    PionHijing->SetLineColor(kBlack);
    PionHijing->SetMarkerColor(kBlack);
    PionHijing->DrawCopy();
    
    PionEnh->SetLineColor(kRed);
    PionEnh->SetMarkerColor(kRed);
    PionEnh->SetMarkerSize(0.3);
    PionEnh->SetMarkerStyle(kOpenCircle);
    PionEnh->DrawCopy("same");
    PionCanvas->BuildLegend();
    PionCanvas->SaveAs("PionDistributions.pdf");
    
    EtaCanvas->cd();
    EtaCanvas->SetLogy();
    EtaHijing->SetMarkerStyle(kOpenCircle);
    EtaHijing->SetMarkerSize(0.3);
    EtaHijing->SetMarkerColor(kBlack);
    EtaHijing->SetLineColor(kBlack);
    EtaHijing->GetYaxis()->SetRangeUser(0.01,10E9);
    EtaHijing->DrawCopy();
    
    EtaEnh->SetLineColor(kRed);
    EtaEnh->SetMarkerStyle(kOpenCircle);
    EtaEnh->SetMarkerColor(kRed);
    EtaEnh->SetMarkerSize(0.3);
    EtaEnh->DrawCopy("same");
    EtaCanvas->BuildLegend();
    EtaCanvas->SaveAs("EtaDistributions.pdf");
    
    PionCanvasW->cd();
    EtaCanvasW->SetLogy();
    PionHijing->Divide(PionEnh);
    PionHijing->Draw();
    PionHijing->SetTitle("Weight for #pi^{0}");
    FileW->cd();
    PionHijing->Write("Pi0W");
    
    EtaCanvasW->cd();
    EtaCanvasW->SetLogy();
    EtaHijing->Divide(EtaEnh);
    EtaHijing->SetTitle("Weight for #eta");
    EtaHijing->Draw();
    FileW->cd();
    EtaHijing->Write("EtaW");
}

void AliHFehpPbTool::CompareMCMotherDistributions(TString FilemTScalling)
{
    TFile *mTScalling = new TFile(FilemTScalling.Data());
    TF1* HagFunctionPi0 = (TF1*) mTScalling->Get("HagFunctionPi0");
    TF1* HagFunctionEta = (TF1*) mTScalling->Get("HagFunctionEta");
    
    TH1F* PionHijing = (TH1F*) fInputList->FindObject("fPtMCpi0_PureHijing");
    TH1F* PionEnh = (TH1F*) fInputList->FindObject("fPtMCpi0_NoMother");
    
    TH1F* EtaHijing = (TH1F*) fInputList->FindObject("fPtMCEta_PureHijing");
    TH1F* EtaEnh = (TH1F*) fInputList->FindObject("fPtMCEta_NoMother");
    
    
    TH1F* RatioPurHijingPion = (TH1F*) PionHijing->Clone("RatioPurHijingPion");
    TH1F* RatioPurHijingEta = (TH1F*) EtaHijing->Clone("RatioPurHijingEta");
    
    for (Int_t i = 1; i <= RatioPurHijingPion->GetNbinsX(); i++)
    {
        RatioPurHijingPion->SetBinContent(i,RatioPurHijingPion->GetBinContent(i)/RatioPurHijingPion->GetBinCenter(i));
        RatioPurHijingPion->SetBinError(i,RatioPurHijingPion->GetBinError(i)/RatioPurHijingPion->GetBinCenter(i));
        
        RatioPurHijingEta->SetBinContent(i,RatioPurHijingEta->GetBinContent(i)/RatioPurHijingEta->GetBinCenter(i));
        RatioPurHijingEta->SetBinError(i,RatioPurHijingEta->GetBinError(i)/RatioPurHijingEta->GetBinCenter(i));
    }
    
    RatioPurHijingPion->Scale(1.,"width");
    RatioPurHijingPion->Divide(HagFunctionPi0);
    
    RatioPurHijingEta->Scale(1.,"width");
    RatioPurHijingEta->SetLineColor(kRed);
    RatioPurHijingEta->SetMarkerColor(kRed);
    RatioPurHijingEta->Divide(HagFunctionEta);
    
    TCanvas *RatiosToMtScale = new TCanvas("Ratio","Ratio",400,300);
    RatiosToMtScale->SetLogy();
    
    RatioPurHijingPion->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingEta->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingPion->Draw("same");
    RatioPurHijingPion->GetYaxis()->SetTitle("Ratio Hijing/m_{T} scaled data");
    RatioPurHijingPion->SetTitle("#pi^{0}");
    RatioPurHijingEta->SetTitle("#eta");
    RatioPurHijingEta->Draw("same");
    
    RatiosToMtScale->BuildLegend();
    
    RatiosToMtScale->SaveAs("RatioHijinToMtScale.pdf");
    
    PionHijing->Add(PionEnh);
    EtaHijing->Add(EtaEnh);
    
    //PionHijing->Scale(1./PionHijing->GetEntries());
    PionHijing->Scale(1.,"width");
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,PionHijing->GetBinError(i)/PionHijing->GetBinCenter(i));
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,EtaHijing->GetBinError(i)/EtaHijing->GetBinCenter(i));
        
    }
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    PionHijing->GetXaxis()->SetRangeUser(0,50);
    PionHijing->Draw();
    
    EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,50);
    EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoData.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
}


void AliHFehpPbTool::CalculateWToData()
{
    Double_t Nevents = GetNEvents();
    
    TF1 *HagFunctionPi0 = new TF1("levy","1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.135*(7.331-2.)))*pow(1.+(sqrt(0.135*0.135+x*x)-0.135)/(7.331*0.1718),-7.331)",0,50);//p-Pb pi0 AllCent
    
    TF1 *HagFunctionEta =new TF1("levy1","0.48*((((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.13498*0.13498+25)-0.13498)/(7.331*0.1718),-7.331)) / (((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+25)-0.13498)/(7.331*0.1718),-7.331)))*(x/sqrt(x*x + 0.54751*0.54751 - 0.13498*0.13498))*1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+x*x)-0.13498)/(7.331*0.1718),-7.331)",0,50);
    
    
    TH1F* PionHijing = (TH1F*) fInputList->FindObject("fPtMCpi0_PureHijing");
    TH1F* PionEnh = (TH1F*) fInputList->FindObject("fPtMCpi0_NoMother");
    
    TH1F* EtaHijing = (TH1F*) fInputList->FindObject("fPtMCEta_PureHijing");
    TH1F* EtaEnh = (TH1F*) fInputList->FindObject("fPtMCEta_NoMother");
    
    
    TH1F* RatioPurHijingPion = (TH1F*) PionHijing->Clone("RatioPurHijingPion");
    TH1F* RatioPurHijingEta = (TH1F*) EtaHijing->Clone("RatioPurHijingEta");
    
    for (Int_t i = 1; i <= RatioPurHijingPion->GetNbinsX(); i++)
    {
        RatioPurHijingPion->SetBinContent(i,RatioPurHijingPion->GetBinContent(i)/RatioPurHijingPion->GetBinCenter(i));
        RatioPurHijingPion->SetBinError(i,RatioPurHijingPion->GetBinError(i)/RatioPurHijingPion->GetBinCenter(i));
        
        RatioPurHijingEta->SetBinContent(i,RatioPurHijingEta->GetBinContent(i)/RatioPurHijingEta->GetBinCenter(i));
        RatioPurHijingEta->SetBinError(i,RatioPurHijingEta->GetBinError(i)/RatioPurHijingEta->GetBinCenter(i));
    }
    
    RatioPurHijingPion->Scale(1.,"width");
    RatioPurHijingPion->Divide(HagFunctionPi0);
    
    RatioPurHijingEta->Scale(1.,"width");
    RatioPurHijingEta->SetLineColor(kRed);
    RatioPurHijingEta->SetMarkerColor(kRed);
    RatioPurHijingEta->Divide(HagFunctionEta);
    
    TCanvas *RatiosToMtScale = new TCanvas("Ratio","Ratio",400,300);
    RatiosToMtScale->SetLogy();
    
    RatioPurHijingPion->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingEta->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingPion->Draw("same");
    RatioPurHijingPion->GetYaxis()->SetTitle("Ratio Hijing/m_{T} scaled data");
    RatioPurHijingPion->SetTitle("#pi^{0}");
    RatioPurHijingEta->SetTitle("#eta");
    RatioPurHijingEta->Draw("same");
    
    RatiosToMtScale->BuildLegend();
    
    RatiosToMtScale->SaveAs("RatioHijinToWeight.pdf");
    
    //PionHijing->Add(PionEnh);
    //EtaHijing->Add(EtaEnh);
    
    //PionHijing->Scale(1./PionHijing->GetEntries());
    PionHijing->Scale(1.,"width");
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        
        //Double_t PionIntegral = HagFunctionPi0->Integral(PionHijing->GetBinLowEdge(i), PionHijing->GetBinLowEdge(i)+ PionHijing->GetBinWidth(i));
        //Double_t EtaIntegral = HagFunctionEta->Integral(EtaHijing->GetBinLowEdge(i), EtaHijing->GetBinLowEdge(i)+ EtaHijing->GetBinWidth(i));
        
        
        Double_t PionIntegral = 1.;
        Double_t EtaIntegral = 1.;
        
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/(PionHijing->GetBinCenter(i)*PionIntegral));
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/(EtaHijing->GetBinCenter(i)*EtaIntegral));
    }
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    PionHijing->GetXaxis()->SetRangeUser(0,50);
    
    
    EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,50);
    EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    /*
     for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
     {
     PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
     EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
     
     }
     */
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Draw();
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoDataNoIntegral.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
}

void AliHFehpPbTool::CalculateWToDataNoEnh()
{
    Double_t bins[83] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.5, 9.0, 9.5, 10.0, 11, 12, 13, 14, 15, 16,18,20,100};
    Int_t Nbins = 82-4;
    
    //Double_t bins[59] = {0.100,0.120,0.140,0.160,0.180,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.200,2.400,2.600,2.800,3.000,3.200,3.400,3.600,3.800,4.000,4.500,5.000,5.500,6.000,6.500,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,18.000,20.00};
    // Int_t Nbins = 58;
    
    //Double_t bins[45] =  {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};
    
    //Int_t Nbins = 44;
    
    Double_t Nevents = GetNEvents();
    
    TF1 *HagFunctionPi0 = new TF1("levy","0.5*1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.135*(7.331-2.)))*pow(1.+(sqrt(0.135*0.135+x*x)-0.135)/(7.331*0.1718),-7.331)",0,100);//p-Pb pi0 AllCent
    
    TF1 *HagFunctionEta =new TF1("levy1","0.48*((((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.13498*0.13498+25)-0.13498)/(7.331*0.1718),-7.331)) / (((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+25)-0.13498)/(7.331*0.1718),-7.331)))*(x/sqrt(x*x + 0.54751*0.54751 - 0.13498*0.13498))*1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+x*x)-0.13498)/(7.331*0.1718),-7.331)",0,100);
    
    
    
    TH1F *temp = (TH1F*) fInputList->FindObject("fPtMCpi0_NoMother");
    TH1F* PionEnh = (TH1F*) temp->Rebin(Nbins,"PionEnh",bins);
    PionEnh->SetTitle("#pi^{0} from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_NoMother");
    TH1F *EtaEnh = (TH1F*) temp->Rebin(Nbins,"EtaEnh",bins);
    EtaEnh->SetTitle("#eta from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_PureHijing");
    TH1F* PionHijing = (TH1F*) temp->Rebin(Nbins,"PionHijing",bins);
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_PureHijing");
    TH1F *EtaHijing = (TH1F*) temp->Rebin(Nbins,"EtaHijing",bins);
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    
    
    TH1F* RatioPurHijingPion = (TH1F*) PionHijing->Clone("RatioPurHijingPion");
    TH1F* RatioPurHijingEta = (TH1F*) EtaHijing->Clone("RatioPurHijingEta");
    
    for (Int_t i = 1; i <= RatioPurHijingPion->GetNbinsX(); i++)
    {
        RatioPurHijingPion->SetBinContent(i,RatioPurHijingPion->GetBinContent(i)/RatioPurHijingPion->GetBinCenter(i));
        RatioPurHijingPion->SetBinError(i,RatioPurHijingPion->GetBinError(i)/RatioPurHijingPion->GetBinCenter(i));
        
        RatioPurHijingEta->SetBinContent(i,RatioPurHijingEta->GetBinContent(i)/RatioPurHijingEta->GetBinCenter(i));
        RatioPurHijingEta->SetBinError(i,RatioPurHijingEta->GetBinError(i)/RatioPurHijingEta->GetBinCenter(i));
    }
    
    RatioPurHijingPion->Scale(1.,"width");
    RatioPurHijingPion->Divide(HagFunctionPi0);
    
    RatioPurHijingEta->Scale(1.,"width");
    RatioPurHijingEta->SetLineColor(kRed);
    RatioPurHijingEta->SetMarkerColor(kRed);
    RatioPurHijingEta->Divide(HagFunctionEta);
    
    TCanvas *RatiosToMtScale = new TCanvas("Ratio","Ratio",400,300);
    RatiosToMtScale->SetLogy();
    
    RatioPurHijingPion->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingEta->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingPion->Draw("same");
    RatioPurHijingPion->GetYaxis()->SetTitle("Ratio Hijing/m_{T} scaled data");
    RatioPurHijingPion->SetTitle("#pi^{0}");
    RatioPurHijingEta->SetTitle("#eta");
    RatioPurHijingEta->Draw("same");
    
    RatiosToMtScale->BuildLegend();
    
    RatiosToMtScale->SaveAs("RatioHijinToWeight.pdf");
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Scale(1.,"width");
    EtaHijing->Scale(1.,"width");
    
    PionHijing->Scale(1./(2*TMath::Pi()));
    EtaHijing->Scale(1./(2*TMath::Pi()));
    
    PionHijing->Scale(1./2.4);
    EtaHijing->Scale(1./2.4);
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,PionHijing->GetBinError(i)/PionHijing->GetBinCenter(i));
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,EtaHijing->GetBinError(i)/EtaHijing->GetBinCenter(i));
        
    }
    
    PionHijing->SaveAs("Pionspectra.root");
    EtaHijing->SaveAs("EtaSpectra.root");
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    PionHijing->GetXaxis()->SetRangeUser(0,50);
    
    
    EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,50);
    EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
        EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
        
    }
    
    
    
    PionHijing->Draw();
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoData_Inverse_Fit.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
    
}

void AliHFehpPbTool::CalculateWToDataUseEnh()
{
    Double_t bins[111] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.5, 9.0, 9.5, 10.0, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,100};
    
    Int_t Nbins = 110;
    
    Double_t Nevents = GetNEvents();
    
    //[0]  = A, [1] = nTsallis, [2] = T, [3] = mass
    
    /*
     TF1 *HagFunctionPi0 = new TF1("pi0", "[0]/2*TMath::Pi() * (([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2))) * pow ((1 + (TMath::Sqrt(x*x + [3]*[3])- [3])/([1]*[2])),-[1])" ,0,50);//p-Pb pi0 AllCent
     
     HagFunctionPi0->SetParameters(8.15719, 7.19545, 0.168252, 134.9766/1000.);
     
     TF1 *HagFunctionEta = new TF1("eta", "[0]/2*TMath::Pi() * (([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2))) * pow ((1 + (TMath::Sqrt(x*x + [3]*[3])- [3])/([1]*[2])),-[1])" ,0,50);
     
     HagFunctionEta->SetParameters(0.909368, 7.50074, 0.269817, 547.862/1000.);
     */
    
    TF1 *HagFunctionPi0 = new TF1("levy","1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.135*(7.331-2.)))*pow(1.+(sqrt(0.135*0.135+x*x)-0.135)/(7.331*0.1718),-7.331)",0,100);//p-Pb pi0 AllCent
    
    TF1 *HagFunctionEta =new TF1("levy1","0.48*((((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.13498*0.13498+25)-0.13498)/(7.331*0.1718),-7.331)) / (((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+25)-0.13498)/(7.331*0.1718),-7.331)))*(x/sqrt(x*x + 0.54751*0.54751 - 0.13498*0.13498))*1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+x*x)-0.13498)/(7.331*0.1718),-7.331)",0,100);
    
    
    
    TH1F *temp = (TH1F*) fInputList->FindObject("fPtMCpi0_NoMother");
    TH1F* PionEnh = (TH1F*) temp->Rebin(Nbins,"PionEnh",bins);
    PionEnh->SetTitle("#pi^{0} from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_NoMother");
    TH1F *EtaEnh = (TH1F*) temp->Rebin(Nbins,"EtaEnh",bins);
    EtaEnh->SetTitle("#eta from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_PureHijing");
    TH1F* PionHijing = (TH1F*) temp->Rebin(Nbins,"PionHijing",bins);
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_PureHijing");
    TH1F *EtaHijing = (TH1F*) temp->Rebin(Nbins,"EtaHijing",bins);
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    
    
    TH1F* RatioPurHijingPion = (TH1F*) PionHijing->Clone("RatioPurHijingPion");
    TH1F* RatioPurHijingEta = (TH1F*) EtaHijing->Clone("RatioPurHijingEta");
    
    for (Int_t i = 1; i <= RatioPurHijingPion->GetNbinsX(); i++)
    {
        RatioPurHijingPion->SetBinContent(i,RatioPurHijingPion->GetBinContent(i)/RatioPurHijingPion->GetBinCenter(i));
        RatioPurHijingPion->SetBinError(i,RatioPurHijingPion->GetBinError(i)/RatioPurHijingPion->GetBinCenter(i));
        
        RatioPurHijingEta->SetBinContent(i,RatioPurHijingEta->GetBinContent(i)/RatioPurHijingEta->GetBinCenter(i));
        RatioPurHijingEta->SetBinError(i,RatioPurHijingEta->GetBinError(i)/RatioPurHijingEta->GetBinCenter(i));
    }
    
    RatioPurHijingPion->Scale(1.,"width");
    RatioPurHijingPion->Divide(HagFunctionPi0);
    
    RatioPurHijingEta->Scale(1.,"width");
    RatioPurHijingEta->SetLineColor(kRed);
    RatioPurHijingEta->SetMarkerColor(kRed);
    RatioPurHijingEta->Divide(HagFunctionEta);
    
    TCanvas *RatiosToMtScale = new TCanvas("Ratio","Ratio",400,300);
    RatiosToMtScale->SetLogy();
    
    RatioPurHijingPion->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingEta->GetXaxis()->SetRangeUser(0,20);
    RatioPurHijingPion->Draw("same");
    RatioPurHijingPion->GetYaxis()->SetTitle("Ratio Hijing/m_{T} scaled data");
    RatioPurHijingPion->SetTitle("#pi^{0}");
    RatioPurHijingEta->SetTitle("#eta");
    RatioPurHijingEta->Draw("same");
    
    RatiosToMtScale->BuildLegend();
    
    RatiosToMtScale->SaveAs("RatioHijinToWeight.pdf");
    
    PionHijing->Add(PionEnh);
    EtaHijing->Add(EtaEnh);
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Scale(1.,"width");
    EtaHijing->Scale(1.,"width");
    
    PionHijing->Scale(1./(2*TMath::Pi()));
    EtaHijing->Scale(1./(2*TMath::Pi()));
    
    PionHijing->Scale(1./1.6);
    EtaHijing->Scale(1./1.6);
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,PionHijing->GetBinError(i)/PionHijing->GetBinCenter(i));
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,EtaHijing->GetBinError(i)/EtaHijing->GetBinCenter(i));
        
    }
    
    PionHijing->SaveAs("Pionspectra.root");
    EtaHijing->SaveAs("EtaSpectra.root");
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    PionHijing->GetXaxis()->SetRangeUser(0,50);
    
    
    EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,50);
    EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    /*
     for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
     {
     PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
     EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
     
     }
     */
    
    PionHijing->Draw();
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoData_Inverse.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
    
}

void AliHFehpPbTool::CalculateWToDataNoFit()
{
    Double_t bins[59] = {0.100,0.120,0.140,0.160,0.180,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.200,2.400,2.600,2.800,3.000,3.200,3.400,3.600,3.800,4.000,4.500,5.000,5.500,6.000,6.500,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,18.000,20.00};
    
    Double_t ybinsDATA[58] = {39.55, 34.89, 30.87, 27.06,23.6, 18.92, 13.96, 10.47, 7.994, 6.194, 4.852, 3.837, 3.081, 2.513, 2.053, 1.701, 1.406, 1.172, 0.9803, 0.8255, 0.6971, 0.5445, 0.3991, 0.2944, 0.2208, 0.1662, 0.1271, 0.09805, 0.07637, 0.05982, 0.04697, 0.03344, 0.02158, 0.01417, 0.009503, 0.006552, 0.004579, 0.003216, 0.002324, 0.001701, 0.00124, 0.000756, 0.0003868, 0.0002116, 0.0001213, 7.17E-05, 4.44E-05, 2.36E-05, 1.10E-05, 5.45E-06, 2.91E-06, 1.64E-06, 9.96E-07, 6.40E-07, 3.99E-07, 2.63E-07, 1.51E-07, 7.32E-08};
    
    Int_t Nbins = 58;
    
    Double_t Nevents = GetNEvents();
    
    //[0]  = A, [1] = nTsallis, [2] = T, [3] = mass
    
    /*
     TF1 *HagFunctionPi0 = new TF1("pi0", "[0]/2*TMath::Pi() * (([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2))) * pow ((1 + (TMath::Sqrt(x*x + [3]*[3])- [3])/([1]*[2])),-[1])" ,0,50);//p-Pb pi0 AllCent
     
     HagFunctionPi0->SetParameters(8.15719, 7.19545, 0.168252, 134.9766/1000.);
     
     TF1 *HagFunctionEta = new TF1("eta", "[0]/2*TMath::Pi() * (([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2))) * pow ((1 + (TMath::Sqrt(x*x + [3]*[3])- [3])/([1]*[2])),-[1])" ,0,50);
     
     HagFunctionEta->SetParameters(0.909368, 7.50074, 0.269817, 547.862/1000.);
     */
    
    
    TH1F *temp = (TH1F*) fInputList->FindObject("fPtMCpi0_NoMother");
    TH1F* PionEnh = (TH1F*) temp->Rebin(Nbins,"PionEnh",bins);
    PionEnh->SetTitle("#pi^{0} from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_NoMother");
    TH1F *EtaEnh = (TH1F*) temp->Rebin(Nbins,"EtaEnh",bins);
    EtaEnh->SetTitle("#eta from enhanced (with no mother);p_{T} (GeV/c);Counts/bin size ");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_PureHijing");
    TH1F* PionHijing = (TH1F*) temp->Rebin(Nbins,"PionHijing",bins);
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_PureHijing");
    TH1F *EtaHijing = (TH1F*) temp->Rebin(Nbins,"EtaHijing",bins);
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    
    //PionHijing->Add(PionEnh);
    //EtaHijing->Add(EtaEnh);
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Scale(1.,"width");
    EtaHijing->Scale(1.,"width");
    
    PionHijing->Scale(1./(2*TMath::Pi()));
    EtaHijing->Scale(1./(2*TMath::Pi()));
    
    PionHijing->Scale(1./2.4);
    EtaHijing->Scale(1./2.4);
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,PionHijing->GetBinError(i)/PionHijing->GetBinCenter(i));
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,EtaHijing->GetBinError(i)/EtaHijing->GetBinCenter(i));
        
    }
    
    
    PionHijing->SaveAs("Pionspectra.root");
    EtaHijing->SaveAs("EtaSpectra.root");
    
    //Divide by data!
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        //PionHijing->SetBinContent(i,(0.5*ybinsDATA[i-1])/PionHijing->GetBinContent(i));
        PionHijing->SetBinError(i,0);
        
        //EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        //EtaHijing->SetBinError(i,0);
        
    }
    
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    //PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    //PionHijing->GetXaxis()->SetRangeUser(0,20);
    
    
    //EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,20);
    //EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    /*
     for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
     {
     PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
     EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
     
     }
     */
    
    PionHijing->Draw();
    //EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoData_Inverse.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
    
}

void AliHFehpPbTool::CalculateWToDataPionFromAndreaNoEnh()
{
    Double_t bins[83] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.5, 9.0, 9.5, 10.0, 11, 12, 13, 14, 15, 16,18,20,100};
    Int_t Nbins = 82-4-1;
    
    //Double_t bins[59] = {0.100,0.120,0.140,0.160,0.180,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.200,2.400,2.600,2.800,3.000,3.200,3.400,3.600,3.800,4.000,4.500,5.000,5.500,6.000,6.500,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,18.000,20.00};
    // Int_t Nbins = 58;
    
    //Double_t bins[45] =  {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};
    
    //Int_t Nbins = 44;
    
    Double_t Nevents = GetNEvents();
    
    TF1 *HagFunctionPi0 = new TF1("levy","1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.135*(7.331-2.)))*pow(1.+(sqrt(0.135*0.135+x*x)-0.135)/(7.331*0.1718),-7.331)",0,100);//p-Pb pi0 AllCent
    
    TF1 *HagFunctionEta =new TF1("levy1","0.48*((((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.13498*0.13498+25)-0.13498)/(7.331*0.1718),-7.331)) / (((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+25)-0.13498)/(7.331*0.1718),-7.331)))*(x/sqrt(x*x + 0.54751*0.54751 - 0.13498*0.13498))*1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+x*x)-0.13498)/(7.331*0.1718),-7.331)",0,100);
    
    
    TH1F *EtaHijing = (TH1F*) fInputList->FindObject("fPtMCEta_kNoMother");
    
    TH1F* temp = (TH1F*) fInputList->FindObject("fPtMCEta_kNoFeedDown");
    EtaHijing->Add(temp);
    
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_kBeauty");
    EtaHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_kCharm");
    EtaHijing->Add(temp);
    
    
    TH1F *PionHijing = (TH1F*) fInputList->FindObject("fPtMCpi0_kNoMother");
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kNoFeedDown");
    PionHijing->Add(temp);
    
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kBeauty");
    PionHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kCharm");
    PionHijing->Add(temp);
    
    
    PionHijing =  (TH1F*) PionHijing->Rebin(Nbins,"PionHijing",bins);
    EtaHijing =  (TH1F*) EtaHijing->Rebin(Nbins,"PionHijing",bins);
    
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Scale(1.,"width");
    EtaHijing->Scale(1.,"width");
    
    PionHijing->Scale(1./(2*TMath::Pi()));
    EtaHijing->Scale(1./(2*TMath::Pi()));
    
    PionHijing->Scale(1./2.4);
    EtaHijing->Scale(1./2.4);
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,PionHijing->GetBinError(i)/PionHijing->GetBinCenter(i));
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,EtaHijing->GetBinError(i)/EtaHijing->GetBinCenter(i));
        
    }
    
    PionHijing->SaveAs("Pionspectra.root");
    EtaHijing->SaveAs("EtaSpectra.root");
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    PionHijing->GetXaxis()->SetRangeUser(0,50);
    
    
    EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,50);
    EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
        EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
        
    }
    
    
    
    PionHijing->Draw();
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoData_Inverse_Fit_NoEnh.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
}

void AliHFehpPbTool::CalculateWToDataNoFitPion()
{
    Double_t bins[61+4+6] = {0.0,0.02,0.04,0.06,0.08,0.100,0.120,0.140,0.160,0.180,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.200,2.400,2.600,2.800,3.000,3.200,3.400,3.600,3.800,4.000,4.500,5.000,5.500,6.000,6.500,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,18.000,20.00,25,30,35,40,45,50,100};
    
    Int_t Nbins = 58+1+1+4+6;
    
    Double_t Nevents = GetNEvents();
    
    TH1F *EtaHijing = (TH1F*) fInputList->FindObject("fPtMCEta_kNoMother");
    
    TH1F* temp = (TH1F*) fInputList->FindObject("fPtMCEta_kLightMesons");
    EtaHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_kNoFeedDown");
    EtaHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_kBeauty");
    EtaHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCEta_kCharm");
    EtaHijing->Add(temp);
    
    
    TH1F *PionHijing = (TH1F*) fInputList->FindObject("fPtMCpi0_kNoMother");

    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kNoFeedDown");
    PionHijing->Add(temp);
    
    //temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kNoIsPrimary");
    //PionHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kLightMesons");
    PionHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kBeauty");
    PionHijing->Add(temp);
    
    temp = (TH1F*) fInputList->FindObject("fPtMCpi0_kCharm");
    PionHijing->Add(temp);
    
    
    PionHijing =  (TH1F*) PionHijing->Rebin(Nbins,"PionHijing",bins);
    EtaHijing =  (TH1F*) EtaHijing->Rebin(Nbins,"EtaHijing",bins);
    
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Scale(1.,"width");
    EtaHijing->Scale(1.,"width");
    
    PionHijing->Scale(1./(2*TMath::Pi()));
    EtaHijing->Scale(1./(2*TMath::Pi()));
    
    PionHijing->Scale(1./1.6);
    EtaHijing->Scale(1./1.6);
    
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,0);
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,0);
        
    }
    
    TF1 *PionFitFunction = new TF1("#pi^{#pm}", "[0]*(([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*pow((1 + (TMath::Sqrt([4]*[4]+x*x)-[3])/([1]*[2])),-[1])",0,100);
    //C, n, T, m(pion), m(meson)
    PionFitFunction->SetParameters(0.5*2.70114,6.79896,1.50142e-01,139.57018/1000.,139.57018/1000.);
    
    TF1 *EtaFitFunction = new TF1("#eta", "[0]*(([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*pow((1 + (TMath::Sqrt([4]*[4]+x*x)-[3])/([1]*[2])),-[1])",0,100);
    //C, n, T, m(pion), m(meson)
    EtaFitFunction->SetParameters(0.48*0.5*2.70114,6.79896,1.50142e-01,139.57018/1000.,547.862/1000.);
    
    PionHijing->SaveAs("Pionspectra.root");
    EtaHijing->SaveAs("EtaSpectra.root");
    PionHijing->Divide(PionFitFunction);
    EtaHijing->Divide(EtaFitFunction);
    
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    
    
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,20);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    /*
     for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
     {
     PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
     EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
     
     }
     */
    
    /*
    
     TFile *ResultsMinjung = new TFile("$ALICE_PHYSICS/PWGHF/hfe/macros/nonHFEcorrect.root");
     TH1F* WeightsRun1 = (TH1F*) ResultsMinjung->Get("hRatio_pPb_5.023TeV_HIJING_pion");
     WeightsRun1->Draw();
     
     TH1F* WeightsRun12 = (TH1F*) ResultsMinjung->Get("hRatio_pPb_5.023TeV_HIJING_eta");
     WeightsRun12->Draw("same");
    */
    
    
    PionHijing->SetLineColor(kGreen);
    EtaHijing->SetLineColor(kMagenta);
    PionHijing->Draw("same");
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundW_Primary.root", "RECREATE");
    PionHijing->Write("Pi0");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
    
}

void AliHFehpPbTool::CalculateWToDataNoFitPionEnh()
{
    Double_t bins[61+4+4] = {0.0,0.02,0.04,0.06,0.08,0.100,0.120,0.140,0.160,0.180,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.200,2.400,2.600,2.800,3.000,3.200,3.400,3.600,3.800,4.000,4.500,5.000,5.500,6.000,6.500,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,18.000,20.00,25,30,35,40,50};
    
    Double_t ybinsDATA[58] = {39.55, 34.89, 30.87, 27.06,23.6, 18.92, 13.96, 10.47, 7.994, 6.194, 4.852, 3.837, 3.081, 2.513, 2.053, 1.701, 1.406, 1.172, 0.9803, 0.8255, 0.6971, 0.5445, 0.3991, 0.2944, 0.2208, 0.1662, 0.1271, 0.09805, 0.07637, 0.05982, 0.04697, 0.03344, 0.02158, 0.01417, 0.009503, 0.006552, 0.004579, 0.003216, 0.002324, 0.001701, 0.00124, 0.000756, 0.0003868, 0.0002116, 0.0001213, 7.17E-05, 4.44E-05, 2.36E-05, 1.10E-05, 5.45E-06, 2.91E-06, 1.64E-06, 9.96E-07, 6.40E-07, 3.99E-07, 2.63E-07, 1.51E-07, 7.32E-08};
    
    Int_t Nbins = 58+1+1+4+4;
    
    Double_t Nevents = GetNEvents();
    
    //[0]  = A, [1] = nTsallis, [2] = T, [3] = mass
    
    /*
     TF1 *HagFunctionPi0 = new TF1("pi0", "[0]/2*TMath::Pi() * (([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2))) * pow ((1 + (TMath::Sqrt(x*x + [3]*[3])- [3])/([1]*[2])),-[1])" ,0,50);//p-Pb pi0 AllCent
     
     HagFunctionPi0->SetParameters(8.15719, 7.19545, 0.168252, 134.9766/1000.);
     
     TF1 *HagFunctionEta = new TF1("eta", "[0]/2*TMath::Pi() * (([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2))) * pow ((1 + (TMath::Sqrt(x*x + [3]*[3])- [3])/([1]*[2])),-[1])" ,0,50);
     
     HagFunctionEta->SetParameters(0.909368, 7.50074, 0.269817, 547.862/1000.);
     */
    
    
    TH1F *EtaHijing = (TH1F*) fInputList->FindObject("fPtMCEta_kNoMother_Enh");
    
    
    TH1F *PionHijing = (TH1F*) fInputList->FindObject("fPtMCpi0_kNoMother_Enh");
    
    //TF1 correc("correc","TMath::Sqrt(0.14*0.14 + x)
    
    
    
    PionHijing =  (TH1F*) PionHijing->Rebin(Nbins,"PionHijing",bins);
    EtaHijing =  (TH1F*) EtaHijing->Rebin(Nbins,"EtaHijing",bins);
    
    PionHijing->SetTitle("#pi from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    EtaHijing->SetTitle("#eta from non-enhanced;p_{T} (GeV/c);Counts/bin size ");
    
    PionHijing->Scale(1./Nevents);
    EtaHijing->Scale(1./Nevents);
    
    PionHijing->Scale(1.,"width");
    EtaHijing->Scale(1.,"width");
    
    PionHijing->Scale(1./(2*TMath::Pi()));
    EtaHijing->Scale(1./(2*TMath::Pi()));
    
    PionHijing->Scale(1./1.6);
    EtaHijing->Scale(1./1.6);
    
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i,PionHijing->GetBinContent(i)/PionHijing->GetBinCenter(i));
        PionHijing->SetBinError(i,0);
        
        EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        EtaHijing->SetBinError(i,EtaHijing->GetBinError(i)/EtaHijing->GetBinCenter(i));
        
    }
    
    TF1 *PionFitFunction = new TF1("#pi^{#pm}", "[0]*(([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*pow((1 + (TMath::Sqrt([4]*[4]+x*x)-[3])/([1]*[2])),-[1])",0,100);
    //C, n, T, m(pion), m(meson)
    PionFitFunction->SetParameters(0.5*2.70114,6.79896,1.50142e-01,139.57018/1000.,139.57018/1000.);
    
    TF1 *EtaFitFunction = new TF1("#eta", "[0]*(([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+[3]*([1]-2)))*pow((1 + (TMath::Sqrt([4]*[4]+x*x)-[3])/([1]*[2])),-[1])",0,100);
    //C, n, T, m(pion), m(meson)
    EtaFitFunction->SetParameters(0.48*0.5*2.70114,6.79896,1.50142e-01,139.57018/1000.,547.862/1000.);
    
    // TF1 *EtaFitFunction =new TF1("levy1","((((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.13498*0.13498+25)-0.13498)/(7.331*0.1718),-7.331)) / (((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+25)-0.13498)/(7.331*0.1718),-7.331)))*(x/sqrt(x*x + 0.54751*0.54751 - 0.13498*0.13498))*1.245*((7.331-1.)*(7.331-2.))/(7.331*0.1718*(7.331*0.1718+0.13498*(7.331-2.)))*pow(1.+(sqrt(0.54751*0.54751+x*x)-0.13498)/(7.331*0.1718),-7.331)",0,100);
    
    PionHijing->SaveAs("Pionspectra.root");
    EtaHijing->SaveAs("EtaSpectra.root");
    PionHijing->Divide(PionFitFunction);
    EtaHijing->Divide(EtaFitFunction);
    
    //Divide by data!
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        //PionHijing->SetBinContent(i,(0.5*ybinsDATA[i-1])/PionHijing->GetBinContent(i));
        PionHijing->SetBinError(i,0);
        
        //EtaHijing->SetBinContent(i,EtaHijing->GetBinContent(i)/EtaHijing->GetBinCenter(i));
        //EtaHijing->SetBinError(i,0);
        
    }
    
    
    
    TCanvas *Weight = new TCanvas("Fit", "Fit", 400,300);
    Weight->SetLogy();
    //PionHijing->Divide(HagFunctionPi0);
    PionHijing->SetTitle("#pi ;p_{T} GeV/c;MC/Data");
    //PionHijing->GetXaxis()->SetRangeUser(0,20);
    
    
    //EtaHijing->Divide(HagFunctionEta);
    EtaHijing->SetTitle("#eta;p_{T} GeV/c;MC/Data");
    EtaHijing->GetXaxis()->SetRangeUser(0,20);
    //EtaHijing->SetLineColor(kRed);
    EtaHijing->SetMarkerColor(kRed);
    
    PionHijing->GetSumw2()->Set(0);
    EtaHijing->GetSumw2()->Set(0);
    
    
    
    for (Int_t i = 1; i <= PionHijing->GetNbinsX(); i++)
    {
        PionHijing->SetBinContent(i, 1./PionHijing->GetBinContent(i));
        EtaHijing->SetBinContent(i,1./EtaHijing->GetBinContent(i));
        
    }
    
    /*
     
     //PionHijing->Scale(1.44);
     TFile *ResultsMinjung = new TFile("$ALICE_PHYSICS/PWGHF/hfe/macros/nonHFEcorrect.root");
     TH1F* WeightsRun1 = (TH1F*) ResultsMinjung->Get("hRatio_pPb_5.023TeV_HIJING_pion");
     WeightsRun1->Draw();
     
     TH1F* WeightsRun12 = (TH1F*) ResultsMinjung->Get("hRatio_pPb_5.023TeV_HIJING_eta");
     WeightsRun12->Draw("same");
     */
    //PionHijing->Scale(1.44);
    //EtaHijing->Scale(1.44);
    
    PionHijing->SetLineColor(kGreen);
    EtaHijing->SetLineColor(kMagenta);
    PionHijing->Draw("same");
    EtaHijing->Draw("same");
    
    Weight->BuildLegend();
    
    TFile *ExportInverseW = new TFile("BackgroundWtoData_Inverse_Enhanced.root", "RECREATE");
    PionHijing->Write("Pi0W");
    EtaHijing->Write("Eta");
    ExportInverseW->Clone();
    
    
    
}


void AliHFehpPbTool::LocalMerge(TString ExportName)
{
    TFile *MergedFile = new TFile(ExportName.Data(),"RECREATE");
    MergedFile->cd();
    
    TString DirectoryName = "HFE_h_";
    DirectoryName.Append(fConfigurationName.Data());
    
    TString ListName = "eh_";
    ListName.Append(fConfigurationName.Data());
    
    TDirectoryFile *dir = new TDirectoryFile(DirectoryName.Data(),DirectoryName.Data());
    dir->cd();
    
    TList *List = new TList();
    List->SetOwner();
    
    TString NamesToCopy[4] = {"fPtElec_Inc","fPtElec_ULS","fPtElec_LS","fPtTrigger_Inc"};
    
    //Copy only
    for (Int_t i = 0; i < 4; i++) {
        TH1F *temp = (TH1F*) fInputList->FindObject(NamesToCopy[i].Data());
        List->Add(temp);
    }
    
    
    TString NamesToMerge[10] = {"fCEtaPhi_Inc", "fCEtaPhi_Inc_EM","fCEtaPhi_Inc_DiHadron", "fCEtaPhi_ULS_Weight", "fCEtaPhi_LS_Weight", "fCEtaPhi_ULS_NoP_Weight", "fCEtaPhi_LS_NoP_Weight", "fCEtaPhi_ULS_Weight_EM", "fCEtaPhi_LS_Weight_EM", "fTOFTPCnsigma_pt"};
    
    Double_t OriginalBins[] = {0.5,0.75,1.0,1.25,1.5,2.0,3.0,4.0,6.0};
    Double_t MergeBins[] = {0.5,1.0,2.0,4.0,6.0};
    
    
    for (Int_t i = 0 ; i < 10; i++ )
    {
        //1st bin
        TString NameFirstBin = NamesToMerge[i] + Form("%d",0);
        TH1F *FirstBin = (TH1F*) fInputList->FindObject(NameFirstBin.Data());
        
        for (Int_t j = 1 ; j < 2; j++)
        {
            TString Name = NamesToMerge[i] + Form("%d",j);
            TH1F *temp = (TH1F*) fInputList->FindObject(Name.Data());
            FirstBin->Add(temp);
        }
        
        FirstBin->SetTitle("0.50 < p_{T}^{e} < 1.00");
        List->Add(FirstBin);
        
        
        //2nd bin
        TString NameBin2 = NamesToMerge[i] + Form("%d",2);
        TH1F *Bin2 = (TH1F*) fInputList->FindObject(NameBin2.Data());
        
        for (Int_t j = 3 ; j < 5; j++)
        {
            TString Name = NamesToMerge[i] + Form("%d",j);
            TH1F *temp = (TH1F*) fInputList->FindObject(Name.Data());
            Bin2->Add(temp);
        }
        
        Bin2->SetTitle("1.00 < p_{T}^{e} < 2.00");
        Bin2->SetName(Form("%s1",NamesToMerge[i].Data()));
        List->Add(Bin2);
        
        //3rd bin
        TString NameBin3 = NamesToMerge[i] + Form("%d",5);
        TH1F *Bin3 = (TH1F*) fInputList->FindObject(NameBin3.Data());
        
        for (Int_t j = 6 ; j < 8; j++)
        {
            TString Name = NamesToMerge[i] + Form("%d",j);
            TH1F *temp = (TH1F*) fInputList->FindObject(Name.Data());
            Bin3->Add(temp);
        }
        Bin3->SetName(Form("%s2",NamesToMerge[i].Data()));
        
        Bin3->SetTitle("2.00 < p_{T}^{e} < 4.00");
        List->Add(Bin3);
        
        
        //3rd bin
        TString NameBin4 = NamesToMerge[i] + Form("%d",8);
        TH1F *Bin4 = (TH1F*) fInputList->FindObject(NameBin4.Data());
        Bin4->SetName(Form("%s3",NamesToMerge[i].Data()));
        
        List->Add(Bin4);
        
        
    }
    
    List->Write(ListName.Data(),TObject::kSingleKey);
    
    MergedFile->Close();
    
    
}




AliHFehpPbTool::~AliHFehpPbTool() {}

#endif

