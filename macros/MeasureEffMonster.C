void MeasureEffMonster()
{
    TFile* FileMann = TFile::Open("~/phd/analysis/monster25/root_files/mannharthist.root", "READ");
    TFile* FileEn = TFile::Open("~/phd/analysis/monster25/root_files/En252Cf.root", "READ");
    TFile* FileEffMon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "RECREATE");
    TH1D* HistMannhart = (TH1D*)FileMann->Get("HistMann");

    TH1D* HistEn = (TH1D*)FileEn->Get("hist_E");

    TGraphErrors* GraphEffMon = new TGraphErrors();
    
    GraphEffMon->SetTitle("GraphEffMon");

    GraphEffMon->AddPointError(0., 0., 0., 0.);

    double step = 0.;
    double Ndecays = 13881479.;

    double branching = 0.;

    double Nemitted = 0.;
    double Ndetected = 0.;

    double EffMon = 0.;

    double Nbins = HistEn->GetNbinsX();
    double Xmax = HistEn->GetXaxis()->GetXmax();

    for(int i = HistEn->GetXaxis()->FindBin(0.3); i <= Nbins; i++)
    {
        Ndetected = HistEn->GetBinContent(i);
        
        branching = HistMannhart->GetBinContent(i);
        
        Nemitted = branching * Ndecays;
    
        EffMon = (Ndetected * (Nbins/Xmax) / (5. * Nemitted));

        GraphEffMon->AddPointError(HistEn->GetBinCenter(i), EffMon, 0, 0.05*EffMon);

        if(HistEn->GetBinCenter(i) > 7.) break;
    }

    GraphEffMon->Draw("APL");

    GraphEffMon->Write("GraphEffMon");

    FileEffMon->Write();
    FileEffMon->Close();
}
