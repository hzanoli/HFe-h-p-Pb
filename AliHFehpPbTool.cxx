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
fEffCorrectionForElectrons(kTRUE),
fPreMerge(kFALSE),
fTPCNSigmaCenter(0),
fTPCNSigmaSTD(0),
fConfigIndex(0),
fCentralityIndex(0),
fLegendTitle(0)
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
fEffCorrectionForElectrons(kTRUE),
fPreMerge(kFALSE),
fTPCNSigmaCenter(0),
fTPCNSigmaSTD(0),
fLegendTitle(0)

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
        
        
        fDiHadron[i]->Scale( fHFepT->GetBinContent(i+1)/DihadronpT->GetBinContent(i+1) * (1-fHadronContamination->GetBinContent(i+1)) );
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
            fHFEh[i]->Scale(1./fEffHFe->GetBinContent(i+1));
        
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




Bool_t AliHFehpPbTool::MergeCorrelationDistributions()
{
    fHFEhNormalized = new TH2F *[fpTBinsResults.GetSize() -1];
    
    //Assuming both always start at the same bin
    Int_t lastj = 0;
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        fHFEhNormalized[i] = (TH2F*) fHFEh[0]->Clone(Form("fHFEhNormalized%d",i));
        fHFEhNormalized[i]->Reset();
        
        for (Int_t j = lastj ; j < fpTBins.GetSize()-1 ; j++)
        {
            if (fpTBinsResults.At(i+1) >= fpTBins.At(j+1))
            {
                fHFEhNormalized[i]->Add(fHFEh[j]);
            }
            else
            {
                lastj = j;
                break;
            }
        }
        fHFEhNormalized[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f",fpTBinsResults.At(i), fpTBinsResults.At(i+1)));
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
    //gStyle->SetOptStat(0);
    
    Bool_t DrawG = kFALSE;
    //if (!DrawG) delete HadronContamination;
    //else
    //  HadronContamination->Divide(4,3);
    
    TH2F **TPCTOFNsigma = new TH2F *[fpTBins.GetSize()];
    TH1F **TPCNsigma = new TH1F *[fpTBins.GetSize()];
    
    fHadronContamination = new TH1F("fHadronContamination","Hadron Contamination", fpTBins.GetSize()-1, fpTBins.GetArray());
    
    fTPCNSigmaCenter = new TH1F("fTPCNSigmaCenter","Mean TPCNSigma Values for electrons", fpTBins.GetSize()-1, fpTBins.GetArray());
    fTPCNSigmaSTD = new TH1F("fTPCNSigmaSTD","Mean TPCNSigma Values for electrons", fpTBins.GetSize()-1, fpTBins.GetArray());
    
    for (Int_t i = 0; i < fpTBins.GetSize() -1; i++)
    {
        TPCTOFNsigma[i] = (TH2F*)fInputList->FindObject(Form("fTOFTPCnsigma_pt%d",i));
        
        if ( (i == (fpTBins.GetSize() - 2))  && (fPreMerge == kTRUE ))
        {
            TPCTOFNsigma[i]->Add( (TH2F*)fInputList->FindObject(Form("fTOFTPCnsigma_pt%d",i+1)));
        }
        TPCTOFNsigma[i]->SetTitle(Form("%1.2f < p_{T}^{e} < %1.2f",fpTBins.At(i), fpTBins.At(i+1)));
        
        if  (DrawG)
        {
            TCanvas *Plot2DTPCTOF = new TCanvas(Form("Plot2DTPCTOF%d",i),Form("Plot2DTPCTOF%d",i), 900,600 );
            Plot2DTPCTOF->SetLogz();
            TPCTOFNsigma[i]->Draw("colz");
            Plot2DTPCTOF->SaveAs(Form("Plot2DTPCTOF_Cent_%d_pT%d_Config%d.pdf",fCentralityIndex, i,fConfigIndex));
        }
        if (!TPCTOFNsigma[i])
        {
            printf("Error reading HadronContamination histogram %d\n", i);
            return kFALSE;
        }
        
        TPCNsigma[i] = (TH1F*) TPCTOFNsigma[i]->ProjectionY(Form("TPCnsigma_pt%d",i), TPCTOFNsigma[i]->GetXaxis()->FindBin(-fTOFNSigma), TPCTOFNsigma[i]->GetXaxis()->FindBin(fTOFNSigma)) ;
        TPCNsigma[i]->GetYaxis()->SetTitle("Counts");
        
        TPCNsigma[i]->Rebin(4);
        if (fpTBins.At(i) > 4.0)
            TPCNsigma[i]->Rebin(2);
        
        
        TF1 *pion = new TF1("Pions", "exp([0]*x)*landau(1)",-12,8);
        
        pion->SetParameters(-4,TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(-6)), -6.1, 2.5);
        
        if (fpTBins.At(i) > 3)
            pion->SetParameters(-4,TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(-3)), -3, 2);
        if (fpTBins.At(i) > 4)
            pion->SetParameters(-4,TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(-4)), -3, 2);
        
        if (fpTBins.At(i) >= 4)
            TPCNsigma[i]->Fit(pion,"RN","", -4, -2);
        else
            TPCNsigma[i]->Fit(pion,"RN","", -7, -2);
        //TPCNsigma[i]->Fit("landau","IR","", -7, -2);
        
        TF1 *protonkaon = new TF1("Protons and Kaons", "gaus(0)", -12,8);
        
        protonkaon->SetParLimits(1, -10,-4);
        if (fpTBins.At(i) > 1.4)
        {
            protonkaon->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(-9)), -9, 1);
            TPCNsigma[i]->Fit(protonkaon,"RN","", -12, -6);
        }
        else
        {
            protonkaon->SetParameters(0, 0, 0);
        }
        
        
        TF1 *elec = new TF1("Electrons", "gaus", -12, 8);
        elec->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(0.)), 0, 1);
        elec->SetParLimits(1, -1,1);
        
        
        TPCNsigma[i]->Fit(elec,"RN","",-0.5, 3);
        
        TF1 *TotalFit;
        
        if (fpTBins.At(i) > 1.4)
        {
            TotalFit = new TF1("Total", "gaus(0) + gaus(3) + exp([6]*x)*landau(7)",-12,8); //Elec + proton/Kaon + pion
            TotalFit->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(0.)), 0, 1, protonkaon->GetParameter(0), protonkaon->GetParameter(1), protonkaon->GetParameter(2), pion->GetParameter(0), pion->GetParameter(1), pion->GetParameter(2), pion->GetParameter(3) );
            //elec parm limi
            TotalFit->SetParLimits(1, -0.2,0.2);
            TotalFit->SetParLimits(2, 0.8, 1.2);
            
            //kaon/proton limits
            if (fpTBins.At(i) > 3.9)
                TotalFit->SetParLimits(4, -10, -6);
            else
                TotalFit->SetParLimits(4, -10, -4);
            
            //pion limits
            
            //TotalFit->SetParLimits(9, -9, -1);
            
        }
        else
        {
            TotalFit = new TF1("Total", "gaus(0) + exp([3]*x)*landau(4)",-10,3.5); //Elec + proton/Kaon + pion
            TotalFit->SetParameters(TPCNsigma[i]->GetBinContent(TPCNsigma[i]->GetXaxis()->FindBin(0.)), 0, 1, pion->GetParameter(0), pion->GetParameter(1), pion->GetParameter(2), pion->GetParameter(3) );
            TotalFit->SetParLimits(1, -0.2,0.2);
            TotalFit->SetParLimits(2, 0.8, 1.2);
            
        }
        
        
        
        TPCNsigma[i]->Fit(TotalFit,"0R");
        elec->SetParameters(TotalFit->GetParameter(0), TotalFit->GetParameter(1), TotalFit->GetParameter(2));
        
        if (fpTBins.At(i) > 1.4)
        {
            protonkaon->SetParameters(TotalFit->GetParameter(3), TotalFit->GetParameter(4), TotalFit->GetParameter(5));
            pion->SetParameters(TotalFit->GetParameter(6), TotalFit->GetParameter(7), TotalFit->GetParameter(8), TotalFit->GetParameter(9));
        }
        else
        {
            pion->SetParameters(TotalFit->GetParameter(3), TotalFit->GetParameter(4), TotalFit->GetParameter(5), TotalFit->GetParameter(6));
        }
        
        
        
        if(DrawG)
        {
            TCanvas *HadronContamination = new TCanvas(Form("HadronContamination%d",i),Form("HadronContamination%d",i), 900,600 );
            HadronContamination->SetLogy();
            TPCNsigma[i]->Draw();
            TPCNsigma[i]->SetLineColor(1);
            TPCNsigma[i]->SetMarkerColor(1);
            TPCNsigma[i]->SetMarkerSize(1);
            TPCNsigma[i]->SetMarkerStyle(20);

            if (fpTBins.At(i) > 1.4)
            {
                protonkaon->SetLineColor(kBlue);
                protonkaon->Draw("same");
                pion->SetLineColor(kMagenta);
                pion->Draw("same");
            }
            else
            {
                pion->SetLineColor(kMagenta);
                pion->Draw("same");
                
            }
            
            elec->Draw("same");
            
            elec->SetLineColor(kGreen);
            TLegend leg(0.65,0.67,0.88,0.88);
            leg.AddEntry(TPCNsigma[i],"Data", "lp");

            if (fpTBins.At(i) > 1.4)
                leg.AddEntry(protonkaon,"Protons/Kaons", "l");
            leg.SetHeader(fLegendTitle.Data());
            leg.AddEntry(pion,"Pions", "l");
            leg.AddEntry(elec,"Electrons", "l");
            //leg.AddEntry((TObject*)0,Form("#chi^2/NDF = %1.2f", TotalFit->GetChisquare()), "l");
            leg.Draw("same");
            //HadronContamination->BuildLegend();
            
            HadronContamination->Print(Form("Contamination_Cent_%d_pT%d_Config%d.pdf",fCentralityIndex, i,fConfigIndex));
            
        }
        
        Double_t ElectronIntegral = elec->Integral(-0.5,3.0);
        Double_t PionIntegral = pion->Integral(-0.5,3.0);
        //Double_t ProtonkaonIntegral= protonkaon->Integral(-0.5,3.0);
        
        Double_t Contamination = (PionIntegral)/(PionIntegral+ElectronIntegral);
        
        //REMOVE THIS FOR PRODUCTION
        fHadronContamination->SetBinContent(i+1,1-Contamination);
        fHadronContamination->SetBinError(i+1,(1-Contamination)/100000000000000000000.);
        
        //fHadronContamination->SetBinContent(i+1,1);
        //fHadronContamination->SetBinError(i+1,(1)/100000000000000000000.);

        
        fTPCNSigmaCenter->SetBinContent(i+1, TotalFit->GetParameter(1));
        fTPCNSigmaCenter->SetBinError(i+1, TotalFit->GetParError(1));
        
        fTPCNSigmaSTD->SetBinContent(i+1, TotalFit->GetParameter(2));
        fTPCNSigmaSTD->SetBinError(i+1, TotalFit->GetParError(2));
        
                                              
        //Double_t BinCounting = TPCNsigma[i]->Integrate(-0.5,3.0);
        printf("Bin: %d, Elec: %1.3f Pion = %1.3f; Total = %1.3f;Cont: = %1.3f \n", i+1, ElectronIntegral, PionIntegral, ElectronIntegral+PionIntegral, Contamination);
        
    }
    
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
        fHFEhSame[i]->GetZaxis()->SetTitle("N_{eh}^{same}");
        fHFEhSame[i]->GetXaxis()->SetTitle("#Delta#varphi(rad)");

        fHFEhSame[i]->Draw("surf1");
        print2Dplot.cd(2);
        FormatTH2Titles(fHFEhMixed[i]);
        fHFEhMixed[i]->Draw("surf1");
        fHFEhMixed[i]->GetZaxis()->SetTitle("N_{eh}^{mixed}/N(0,0)");
        fHFEhMixed[i]->GetXaxis()->SetTitle("#Delta#varphi(rad)");
        print2Dplot.cd(3);
        FormatTH2Titles(fHFEh[i]);
        fHFEh[i]->Draw("surf1");
        fHFEh[i]->GetZaxis()->SetTitle("N_{eh}/#varepsilon_{tracking}");
        fHFEh[i]->GetXaxis()->SetTitle("#Delta#varphi(rad)");
        
        print2Dplot.Print(Form("HfeEffCorrected_%d_%d_%d.pdf",fConfigIndex,fCentralityIndex,i));
    }
    
    for (int i = 0; i < fpTBinsResults.GetSize() - 1 ; i++)
    {
        TCanvas merged("merged","merged",400,300);
        merged.cd();
        FormatTH2Titles(fHFEhNormalized[i]);
        fHFEhNormalized[i]->Draw("surf1");
        fHFEhNormalized[i]->GetZaxis()->SetTitle("N_{eh}/N_{e}");
        merged.Print(Form("Hfe_merged_%d_%d_%d.pdf",fConfigIndex,fCentralityIndex,i));
    }

       //fHFEhNormalized[pT]; //
}

void AliHFehpPbTool::FormatTH2Titles(TH2F* histo)
{
    histo->GetYaxis()->SetRangeUser(-1.6,1.6);
    histo->GetXaxis()->SetTitleSize(0.06);
    histo->GetXaxis()->SetTitleOffset(1.3);
    histo->GetYaxis()->SetTitleSize(0.06);
    histo->GetYaxis()->SetTitleOffset(1.3);
    histo->GetZaxis()->SetTitleSize(0.05);
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
            FlowFunction[i]->SetParameters(10.,0.05, 0.1);
            //FlowFunction[i]->FixParameter(1,0);
        }

        JetReference[i] = Reference->GetHFeh1DSub(i);

        fFlowHistograms[i] = (TH1F*) fHFEhNormalized1D[i]->Clone(Form("FlowHistograms%d",i));
        fFlowHistograms[i]->Add(JetReference[i],-1);
        fFlowHistograms[i]->Fit(FlowFunction[i],"0");

        
    }
    
    
}


Bool_t AliHFehpPbTool::NormalizeHFeCorrrelation()
{
    TH1F *HFepTRebin = (TH1F*)fHFepT->Rebin(fpTBinsResults.GetSize()-1, "HFepTRebin", fpTBinsResults.GetArray());
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
        fHFEhNormalized[i]->Scale(1./HFepTRebin->GetBinContent(i+1));
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
    fHFEhNormalized1D = new TH1F *[fpTBinsResults.GetSize() -1];
    fHFEhNormSub1D = new TH1F *[fpTBinsResults.GetSize() -1];
    
    fBaseline = new TH1F("fBaseline", "Baseline", fpTBinsResults.GetSize() -1, fpTBinsResults.GetArray());
    
    for (Int_t i = 0 ; i < fpTBinsResults.GetSize() -1 ; i++ )
    {
        fHFEhNormalized1D[i] = (TH1F*) fHFEhNormalized[i]->ProjectionX(Form("fHFEhNormalized1D%d",i),1, fHFEhNormalized[i]->GetNbinsY());
        fHFEhNormalized1D[i]->Scale(1.,"width");
        //fHFEhNormalized1D[i]->Scale(1.0/3.2); //Full eta
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
        
        Double_t Pedestal = (MinValue[0] + MinValue[1] + MinValue[2])/3.;
        
        TF1 *PedestalFunction = new TF1("PedestalFunction",Form("%1.4f",Pedestal), -3,7);
        
        //Correlation
        
        //TF1 *CorrelationNoFlow = new TF1("CorrelationNoFlow", "[0] + gaus(1) + gaus(4)", -TMath::Pi()/2.0, (3.0/2.0)*TMath::Pi());
        //CorrelationNoFlow->SetParameters(Pedestal, fHFEhNormalized1D[i]->GetBinContent( fHFEhNormalized1D[i]->GetXaxis()->FindBin(0.)) , 0., 0.5, fHFEhNormalized1D[i]->GetBinContent( fHFEhNormalized1D[i]->GetXaxis()->FindBin(TMath::Pi())), TMath::Pi(), 0.5 );
        
        //fHFEhNormalized1D[i]->Fit(CorrelationNoFlow, "0");
        
        fHFEhNormSub1D[i]->Add(PedestalFunction,-1);
        
        fBaseline->SetBinContent(i+1,Pedestal);
        
        
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
    TF1 **FlowFunction = new TF1 *[fpTBinsResults.GetSize()];
    
    for (int pT = 0; pT < fpTBinsResults.GetSize() -1; pT++)
    {
        fHFehProjectionForV2NonSub[pT] = (TH1F*)fHFEhNormalized1D[pT]->Clone(Form("fHFehProjectionForV2NonSub%d",pT));
        fHFehProjectionForV2NonSub[pT]->Reset();
        //fBackNonIDEh[pT]->Divide(fBackMixedEh[pT]);
        //fHFEhNormalized[pT] = fBackNonIDEh[pT];
        //fHFEhNormalized[pT]->Draw("surf1");
        
        /*
         //Near side: only delta eta > 0.8 considered
         for (Int_t ix = 1 ; ix <= (int)fHFEhNormalized[pT]->GetNbinsX()/2; ix++)
         {
         Double_t Mean;
         Double_t Normalization;
         
         // -1.6 < delta eta  < -0.8
         for (Int_t iy = fHFEhNormalized[pT]->GetYaxis()->FindBin(-1.55) ; iy <= fHFEhNormalized[pT]->GetYaxis()->FindBin(-0.85) ; iy++)
         {
         Mean += fHFEhNormalized[pT]->GetBinContent(ix,iy)/( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
         Normalization += 1./( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
         
         //Mean += fHFEhNormalized[pT]->GetBinContent(ix,iy);
         //Normalization += 1;
         }
         
         // 0.8 < delta eta  < 1.6
         for (Int_t iy = fHFEhNormalized[pT]->GetYaxis()->FindBin(0.85) ; iy <= fHFEhNormalized[pT]->GetYaxis()->FindBin(1.55) ; iy++)
         {
         Mean += fHFEhNormalized[pT]->GetBinContent(ix,iy)/( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
         Normalization += 1./( pow(fHFEhNormalized[pT]->GetBinError(ix,iy),2) );
         }
         
         fHFehProjectionForV2NonSub[pT]->SetBinContent(ix, Mean/Normalization);
         fHFehProjectionForV2NonSub[pT]->SetBinError(ix, TMath::Sqrt(1./Normalization));
         
         }
         */
        
        //away side: all eta cosidered
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
        
        FlowFunction[pT] = new TF1(Form("FitV2NoSub%d",pT),"[0]*(1 + 2 * [1] * TMath::Cos(x) + 2 * [2] * TMath::Cos(2*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        //FlowFunction[pT]->FixParameter(1,0);
        fHFehProjectionForV2NonSub[pT]->Fit(FlowFunction[pT], "0");
        
        //fHFehProjectionForV2NonSub[pT]->Draw();
        
    }
    
    
}


AliHFehpPbTool::~AliHFehpPbTool() {}

#endif

