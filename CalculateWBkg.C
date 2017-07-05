void CalculateWBkg(TString NameDataInput)
{
    gStyle->SetOptStat(0);
    gStyle->SetLegendFont(42);
    TH1::AddDirectory(kFALSE);
    TGaxis::SetMaxDigits(4);
    gROOT->LoadMacro(Form("%s/AliHFehpPbTool.cxx++g",gSystem->Getenv("PathToGitAnalysis")));
    
    //Basic Configuration:
    Double_t pTBinsCorrelation[] = {0.5,0.75,1.0,1.25,1.5,2,2.5,3,4,5,6};
    Double_t pTBinsResults[] = {0.5,1.0,2.0,4.0,6.0};
    
    Int_t NumerOfpTBins = (Int_t) sizeof(pTBinsResults)/sizeof(pTBinsResults[0]) - 1 ;
    
    //Configurations, Centrality
    AliHFehpPbTool* Tool;
    
    Tool = new AliHFehpPbTool("Tool", "Tool");
        
    TString ConfigName = EasyName(0,ConfigurationsArray.At(Config),kTRUE);
    MC[Config]->ConnectToInputFile(NameMCInput,ConfigName);
        
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
     
     //Min Op Angle
     50 = <0.1
     51 = <0.15
     
     TPC NClus
     60 = 80 Clus
     61 = 90 Clus
     
     TPC PID
     70 = 0 to 3 sigma
     
     
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
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 1:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,4,1,100,60,80);
                break;
            case 2:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,0,100,60,80);
                break;
            case 3:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,4,0,100,60,80);
                break;
            case 4:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,4,2,100,60,80);
                break;
            case 5:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,2,100,60,80);
                break;
                
                //DCA + ITS
            case 10:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,2.4,3.2,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,0,100,60,80);
                break;
            case 11:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,2.4,3.2,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
                //Inv Mass
            case 20:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.10,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 21:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.12,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
                //Min pT partner
            case 30:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.20,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 31:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.25,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 32:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.30,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
                //Hadron DCA
            case 40:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.50,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 41:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.75,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 42:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,1.00,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
                //Op angle
            case 50:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,0.1,2,1,100,60,80);
                break;
            case 51:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,0.15,2,1,100,60,80);
                break;
                
                //TPC NClusters
            case 60:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,0.15,2,1,90,60,80 );
                break;
            case 61:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,80,60,80);
                break;
                //TPC NSgima
            case 70:
                Name = AliHFehpPbTool::TaskName(0,1,0,0,0.25,1.00,0.25,1.00,0.00,3.00,0.14,0.00,CentralityMin,CentralityMax,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
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
                
                //Min pT Asso
            case 30:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.2,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 31:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.25,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 32:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.3,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
                //Hadron DCA
            case 40:
                Name = AliHFehpPbTool::TaskName(0,1,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 41:
                Name = AliHFehpPbTool::TaskName(0,1,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
            case 42:
                Name = AliHFehpPbTool::TaskName(0,1,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
                //Op angle
            case 50:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,0.1,2,1,100,60,80);
                break;
            case 51:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,0.15,2,1,100,60,80);
                break;
                
                //TPC NClusters
            case 60:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,90,60,80);
                break;
            case 61:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,80,60,80);
                break;
                //TPC NSgima
            case 70:
                Name = AliHFehpPbTool::TaskName(2,0,0,1,0.25,1.00,0.25,1.00,0.00,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.80,0.80,999.00,2,1,100,60,80);
                break;
                
            case 10000:
                Name = AliHFehpPbTool::TaskName(1,0,0,1,1.00,2.00,0.25,1.00,-0.50,3.00,0.14,0.00,0.00,100.00,0.30,2.00,-0.60,0.60,999.00,4,1,100,60,80);
                break;
                
        }
        
        
    }
    return Name;
}

TString NamesToLeg(Int_t Config = 0)
{
    TString Name;
    switch (Config) {
        case 0:
            Name = "Default";
            break;
        case 1:
            Name = "SPD First, 4 ITS hits";
            break;
        case 2:
            Name = "SPD Both, 2 ITS hits";
            break;
        case 3:
            Name = "SPD Both, 4 ITS hits";
            break;
        case 4:
            Name = "SPD Any, 4 ITS hits";
            break;
        case 5:
            Name = "SPD Any, 2 ITS hits";
            break;
            
        case 20:
            Name = "Bkg Inv Mass < 100 GeV/c^2";
            break;
        case 21:
            Name = "Bkg Inv Mass < 120 GeV/c^2";
            break;
            
            
        case 30:
            Name = "pT partner > 0.20 GeV/c";
            break;
        case 31:
            Name = "pT partner > 0.25 GeV/c";
            break;
        case 32:
            Name = "pT partner > 0.30 GeV/c";
            break;
            
            //Hadron DCA
        case 40:
            Name = "DCAxy = 0.50cm";
            break;
        case 41:
            Name = "DCAxy = 0.75cm";
            break;
        case 42:
            Name = "DCAxy = 1.00cm";
            break;
            
            //Op angle
        case 50:
            Name = "Op. Angle < 0.1 rad";
            break;
        case 51:
            Name = "Op. Angle < 0.15 rad";
            break;
            
            //TPC NClusters
        case 60:
            Name = "90 N. Clust TPC";
            break;
        case 61:
            Name = "80 N Clust TPC";
            break;
            //TPC NSgima
        case 70:
            Name = "TPC N#Sigma from 0 to 3";
            break;
            
        default:
            Name = "No Config Name: Add please";
            break;
    }
    return Name;
}
Int_t Color(Int_t Config)
{
    /// return Config+1;
    
    switch (Config)
    {
            //ITS
        case 0:
            return 1;
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
            
        case 40:
            return 2;
        case 41:
            return 3;
        case 42:
            return 4;
            
            
    }
}

void BasicFormat (TH1* histo, Int_t Color = 0, Int_t Marker = 20)
{
    histo->SetMarkerStyle(Marker);
    histo->SetMarkerSize(0.5);
    histo->SetLineColor(Color);
    histo->SetMarkerColor(Color);
    histo->GetYaxis()->SetTitleOffset(1.2);
    
}

void BasicFormat (TGraphErrors* histo, Int_t Color = 0, Int_t Marker = 20)
{
    histo->SetMarkerStyle(Marker);
    histo->SetLineColor(Color);
    histo->SetMarkerColor(Color);
    
}

void BasicFormat (TF1* func, Int_t Color = 0, Int_t Marker = 20)
{
    func->SetLineColor(Color);
}
