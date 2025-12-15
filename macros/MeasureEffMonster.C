void MeasureEffMonster()
{
    TFile* FileMann = TFile::Open("~/phd/analysis/monster25/root_files/manhart.root", "READ");
    TFile* FileEn = TFile::Open("~/phd/analysis/monster25/root_files/En252Cf.root", "READ");
    TFile* FileEffMon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "RECREATE");
    TGraph* GraphMann = (TGraph*)FileMann->Get("g1");

    TH1D* HistEn = (TH1D*)FileEn->Get("hist_E");

    TGraphErrors* GraphEffMon = new TGraphErrors();
    
    GraphEffMon->SetTitle("GraphEffMon");

    GraphEffMon->AddPointError(0., 0., 0., 0.);

    double step = 0.;
    double Ndecays = 3881479.;

    double branching = 0.;

    double Nemitted = 0.;
    double Ndetected = 0.;

    double EffMon = 0.;

    for(int i = HistEn->GetXaxis()->FindBin(0.3); i <= HistEn->GetNbinsX(); i++)
    {
        Ndetected = HistEn->GetBinContent(i);
        
        branching = GraphMann->Eval(HistEn->GetBinCenter(i));
        
        Nemitted = branching * Ndecays;
    
        EffMon = (Ndetected / (0.01 * Nemitted));

        GraphEffMon->AddPointError(HistEn->GetBinCenter(i), EffMon, 0, 0.05*EffMon);

        if(HistEn->GetBinCenter(i) > 7.) break;
    }

    GraphEffMon->Draw("APL");

    GraphEffMon->Write("GraphEffMon");

    FileEffMon->Write();
    FileEffMon->Close();
}
