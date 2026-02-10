void MeasureEffMonster()
{
    TFile* FileMann = TFile::Open("~/phd/analysis/monster25/root_files/mannhartgraph.root", "READ");
    TFile* FileEn = TFile::Open("~/phd/analysis/monster25/root_files/En252Cf.root", "READ");
    TFile* FileEffMon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "RECREATE");

    TGraph* GraphMannhart = (TGraph*)FileMann->Get("mannhart");

    TH1D* HistEn = (TH1D*)FileEn->Get("hist_E");

    const int Nbins = 177;
    double Xmax = HistEn->GetXaxis()->GetXmax();

    TGraphErrors* GraphEffMon = new TGraphErrors();
        
    GraphEffMon->SetTitle("GraphEffMon");

    GraphEffMon->AddPointError(0., 0., 0., 0.);

    double step = 0.;
    double Ndecays = 13881479.;

    double branching = 0.;

    double Nemitted = 0.;
    double Ndetected = 0.;

    double EffMon = 0.;

    double ErrEffMon = 0.;

    double x[Nbins]; 
    double y[Nbins]; 
    double ymin[Nbins]; 
    double ymax[Nbins]; 

    for(int i = 0; i < Nbins; i++)
    {
        Ndetected = HistEn->GetBinContent(i);
        
        branching = GraphMannhart->Eval(HistEn->GetBinCenter(i));

        Nemitted = branching * Ndecays;
    
        EffMon = (Ndetected / (0.03 * Nemitted));
        ErrEffMon = pow(HistEn->GetBinError(i)/(Nemitted*0.03),2) + pow(HistEn->GetBinContent(i)*0.3/(Nemitted*Ndecays*0.03),2);

        x[i] = HistEn->GetBinCenter(i);

        y[i] = EffMon * 100.;;
        ymin[i] = (EffMon - (sqrt(ErrEffMon)/2.)) * 100.;
        ymax[i] = (EffMon + (sqrt(ErrEffMon)/2.)) * 100.;

        GraphEffMon->AddPointError(HistEn->GetBinCenter(i), EffMon, -4.83e-3+0.0416*HistEn->GetBinCenter(i)+3.82e-3*pow(HistEn->GetBinCenter(i),2), sqrt(ErrEffMon));
    }

    TGraph* GraphEffMons = new TGraph(Nbins, x, y);
    TGraph* GraphEffMonMin = new TGraph(Nbins, x, ymin);
    TGraph* GraphEffMonMax = new TGraph(Nbins, x, ymax);
    TGraph* GraphEffMonShade = new TGraph(2 * Nbins);
    
        
    for(int i = 0; i < Nbins; i++)
    {
        GraphEffMonShade->SetPoint(i, x[i], ymax[i]);       
        GraphEffMonShade->SetPoint(Nbins+i, x[Nbins-i-1], ymin[Nbins-i-1]);       
    }

    GraphEffMonShade->SetFillColorAlpha(kRed, 0.5);
    GraphEffMonShade->Draw("AF");

    GraphEffMons->SetLineColor(2);
    GraphEffMons->SetLineWidth(2);
    GraphEffMons->SetMarkerColor(2);
    GraphEffMons->SetMarkerSize(1);
    GraphEffMons->SetFillStyle(1001);
    GraphEffMons->SetFillColorAlpha(2, 0.5);

    GraphEffMons->Draw("CP");

    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);

    legend->AddEntry(GraphEffMons,"MONSTER efficiency","plf");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);

    legend->Draw();

    GraphEffMon->Write("GraphEffMon");

    FileEffMon->Write();
    FileEffMon->Close();
}
