void Systematics(TString NameMCInput = "MC.root", TString NameDataInput = "DATA.root")
{
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    gROOT->LoadMacro("AliHFehpPbTool.cxx++g");
    
    //Basic Configuration:
    Double_t pTBinsCorrelation[] = {0.5,0.75,1.0,1.25,1.5,2,2.5,3,4,5,6};
    Double_t pTBinsResults[] = {0.5,1,2,4,6};
    
    Int_t NumerOfpTBins = (Int_t) sizeof(pTBinsResults)/sizeof(pTBinsResults[0]) - 1 ;
    
    //Configurations, Centrality
    AliHFehpPbTool* Data[20][3];
    TH1F* TagEff[20];
    TH1F* HFeEff[20];
    TH1F* pTHFe[20][3];
    TH1F* ElectronPurity[20][3];
    TH1F* HFeh1D[20][3][8]; //last is pT bin
    
    //Centrality Independent corrections
    AliHFehpPbTool* MC[20];
   // Int_t Array[] = {2,3,4,5};
    Int_t Array[] = {0,1,2,3,4,5};
   //Int_t Array[] = {10000};
    TArrayI ConfigurationsArray ((Int_t) sizeof(Array)/sizeof(Array[0]),Array );
    
    TCanvas *TaggingEff = new TCanvas("TagEff", "TagEff",800,600);
    TCanvas *HFeEffC = new TCanvas("HFeEff", "HFeEff",800,600);
    TCanvas *HFepT = new TCanvas("HFepT", "HFepT",800,600);
    TCanvas *Purity = new TCanvas("Purity", "Purity", 800,600);
    TCanvas *Corre1D = new TCanvas("Corre1D", "Corre1D", 1200,800);
    Corre1D->Divide(3,2);
    
    for (Int_t Config = 0 ; Config < ConfigurationsArray.GetSize() ; Config++ )
    {
        MC[Config] = new AliHFehpPbTool(Form("MC%d",Config), Form("MC%d",Config));
        MC[Config]->SetpTBins((Int_t) sizeof(pTBinsCorrelation)/sizeof(pTBinsCorrelation[0]),pTBinsCorrelation);
        MC[Config]->SetpTBinsResults((Int_t) sizeof(pTBinsResults)/sizeof(pTBinsResults[0]),pTBinsResults);
        
        
        TString ConfigName = EasyName(0,ConfigurationsArray.At(Config),kTRUE);
        //printf("%s\n",ConfigName.Data());
        MC[Config]->ConnectToInputFile(NameMCInput,ConfigName);
        MC[Config]->CalculateTaggingEfficiencyW();
        
        
        //Plot the Tagging Efficiencies for all configurations
        TagEff[Config] = MC[Config]->GetEff();
        TagEff[Config]->GetYaxis()->SetRangeUser(0,0.8);
        BasicFormat(TagEff[Config], Color(ConfigurationsArray.At(Config)));
        TaggingEff->cd();
        TagEff[Config]->Draw("same");
        
        //Plot the Reco Efficiencies for all configurations
        
        MC[Config]->CalculateHFeEfficiency();
        HFeEff[Config]= MC[Config]->GetHFeEff();
        BasicFormat(HFeEff[Config], Color(ConfigurationsArray.At(Config)));
        HFeEff[Config]->GetYaxis()->SetRangeUser(0,0.5);
        HFeEffC->cd();
        HFeEff[Config]->Draw("same");
        
        
        for (Int_t Centrality = 0 ; Centrality < 1 ; Centrality++)
        {
            Data[Config][Centrality] = new AliHFehpPbTool(Form("Data%d_%d",Config,Centrality), Form("Data%d_%d",Config,Centrality));
            Data[Config][Centrality]->SetUseEffElectrons();
            Data[Config][Centrality]->ConnectToInputFile(NameDataInput,EasyName(Centrality,ConfigurationsArray.At(Config),kFALSE));
            
            Data[Config][Centrality]->SetpTBins((Int_t) sizeof(pTBinsCorrelation)/sizeof(pTBinsCorrelation[0]),pTBinsCorrelation);
            Data[Config][Centrality]->SetpTBinsResults((Int_t) sizeof(pTBinsResults)/sizeof(pTBinsResults[0]),pTBinsResults);
            
            Data[Config][Centrality]->SetTagEff(TagEff[Config]);
            Data[Config][Centrality]->SetHFeEff(HFeEff[Config]);
            
            //Process Data...
            Data[Config][Centrality]->Process();
            /*
            Data[Config][Centrality]->ReadpTHistograms();
            Data[Config][Centrality]->CalculateHadronContamination();
            Data[Config][Centrality]->ProcesspTHistograms();
            Data[Config][Centrality]->ReadAndProcessCorrelationDistributions();
            Data[Config][Centrality]->MergeCorrelationDistributions();
            Data[Config][Centrality]->NormalizeHFeCorrrelation();
            Data[Config][Centrality]->ProjectTo1D();
            */
            
            for (Int_t pT = 0 ; pT < NumerOfpTBins; pT++)
            {
                HFeh1D[Config][Centrality][pT] = Data[Config][Centrality]->GetHFeh1D(pT);
                BasicFormat(HFeh1D[Config][Centrality][pT], Color(ConfigurationsArray.At(Config)));
                Corre1D->cd(pT+1);
                // HFeh1D[Config][Centrality][pT]->Divide(HFeh1D[Config][Centrality][pT],HFeh1D[0][Centrality][pT],1,1,"B");
                HFeh1D[Config][Centrality][pT]->Draw("same");
                
                
            }
            
            pTHFe[Config][Centrality] =  Data[Config][Centrality]->GetHFepT();
            BasicFormat(pTHFe[Config][Centrality], Color(ConfigurationsArray.At(Config)));
            
            HFepT->cd();
            pTHFe[Config][Centrality]->Draw("same");
            
            Purity->cd();
            ElectronPurity[Config][Centrality] = Data[Config][Centrality]->GetPurity();
            BasicFormat(ElectronPurity[Config][Centrality], Color(ConfigurationsArray.At(Config)));
            ElectronPurity[Config][Centrality]->Draw("same");
            
            
            

        }
        
    
        
    }
    

    
}

TString EasyName(Int_t Centrality = 0, Int_t Config = 0, Bool_t isMC = kFALSE)
{
    /*Config
     
     //Tracking!
     ITS variatins
     0 = Default (First 2)
     1 = First 4
     2 = Both 2
     3 = Both 4
     4 = Any 4
     5 = Any 2
    
     
     DCA + ITS variations
     
     10 = Both 2 No DCA cut
     11 = First 2 No DCA
     
     Inv Mass
     20 = 100 MeV
     21 = 120 MeV
     
     Min pT Asso
     
     30 = 0.20 GeV
     31 = 0.25 GeV
     32 = 0.30 GeV
     
     Hadron DCA xy variation
     use default effficiencis (hadrons are corrected online)
     40 = 0.5cm
     41 = 0.75cm
     42 = 1.00cm
     
     
     */

    Double_t CentralityMin = 0.;
    Double_t CentralityMax = 100.;
    if (Centrality == 0)
        CentralityMax = 20;
    else if (Centrality == 1)
    {
        CentralityMin = 20.;
        CentralityMax = 60.;
    }
    else if (Centrality == 2)
    {
        CentralityMin = 60.;
        CentralityMax = 100.;

    }
        TString Name;
    if (!isMC)
    {
        //DATA
        switch (Config) {
            //ITS
            case 0:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 1:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,4,1,100,60,80);
                break;
            case 2:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,0,100,60,80);
                break;
            case 3:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,4,0,100,60,80);
                break;
            case 4:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,4,2,100,60,80);
                break;
            case 5:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,2,100,60,80);
                break;
                
            //DCA + ITS
            case 10:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,2.4,3.2,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,0,100,60,80);
                break;
            case 11:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,2.4,3.2,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            
            //Inv Mass
            case 20:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.10,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 21:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.12,0.00,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
            //Min pT partner
            case 30:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.20,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 31:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.25,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 32:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.30,0.00,20.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
        }
    }
    else
    {
        //MC
        switch (Config) {
            //ITS
            
            case 0:
                Name = AliHFehpPbTool::TaskName(0,1,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 1:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,4,1,100,60,80);
                break;
            case 2:
                Name = AliHFehpPbTool::TaskName(0,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,0,100,60,80);
                break;
            case 3:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,4,0,100,60,80);
                break;
            case 4:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,4,2,100,60,80);
                break;
            case 5:
                Name = AliHFehpPbTool::TaskName(0,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,2,100,60,80);
                break;

            //DCA + ITS
            case 10:
                Name = AliHFehpPbTool::TaskName(2,1,0,1,2.4,3.2,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,0,100,60,80);
                break;
            case 11:
                Name = AliHFehpPbTool::TaskName(2,1,0,1,2.4,3.2,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            
            //Inv Mass
            case 20:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.10,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 21:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.12,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            
                
            case 30:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.2,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 31:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.25,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 32:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.3,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
            case 10000:
                Name = AliHFehpPbTool::TaskName(1,0,0,1,1.00,2.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.60,0.60,999.00,4,1,100,60,80);
                break;

                

                
        }

        
    }
    return Name;
}

Int_t Color(Int_t Config)
{
    return Config+1;
    
    switch (Config)
    {
        //ITS
        case 0:
            return TColor::GetColor("#000000");
        case 1:
            return TColor::GetColor("#99d8c9");
        case 2:
            return TColor::GetColor("#66c2a4");
        case 3:
            return TColor::GetColor("#41ae76");
        case 4:
            return TColor::GetColor("#238b45");
        case 5:
            return TColor::GetColor("#006d2c");
        case 6:
            return TColor::GetColor("#00441b");
        
        // Inv mass
        case 20:
            return TColor::GetColor("#810f7c");
        case 21:
            return TColor::GetColor("#4d004b");
        
        //Min pT Ass
        case 30:
            return TColor::GetColor("#4eb3d3");
        case 31:
            return TColor::GetColor("#2b8cbe");
        case 32:
            return TColor::GetColor("#08589e");
            
        
    }
}

void BasicFormat (TH1* histo, Int_t Color = 0, Int_t Marker = 20)
{
    histo->SetMarkerStyle(Marker);
    histo->SetLineColor(Color);
    histo->SetMarkerColor(Color);
    
}
