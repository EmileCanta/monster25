void MeasureEffMonster()
{
    TFile* FileMann = TFile::Open("~/phd/analysis/monster25/root_files/mannhartgraph.root", "READ");
    TFile* FileEn = TFile::Open("~/phd/analysis/monster25/root_files/En252Cf.root", "READ");
    TFile* FileEffMon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "RECREATE");

    TGraph* GraphMannhart = (TGraph*)FileMann->Get("mannhart");

    TH1D* HistEn = (TH1D*)FileEn->Get("hist_E");

    int Nbins = HistEn->GetNbinsX();

    TGraphErrors* GraphEffMon = new TGraphErrors();
        
    GraphEffMon->SetTitle("GraphEffMon");

    double step = 0.;
    double Ndecays = 13881479.;

    double branching = 0.;

    double Nemitted = 0.;
    double Ndetected = 0.;

    double EffMon = 0.;

    double ErrEffMon = 0.;

    vector<double> x; 
    vector<double> y; 
    vector<double> ymin; 
    vector<double> ymax; 
    
    TGraph* GraphEffMonPlot = new TGraph();
    TGraph* GraphEffMonMin = new TGraph();
    TGraph* GraphEffMonMax = new TGraph();
    TGraph* GraphEffMonShade = new TGraph();

    for(int i = 1; i <= Nbins; i++)
    {
        Ndetected = HistEn->GetBinContent(i);
        
        branching = GraphMannhart->Eval(HistEn->GetBinCenter(i));

        Nemitted = branching * Ndecays;
    
        EffMon = (Ndetected / (0.03 * Nemitted));
        ErrEffMon = pow(HistEn->GetBinError(i)/(Nemitted*0.03),2) + pow(HistEn->GetBinContent(i)*0.3/(Nemitted*Ndecays*0.03),2);

        GraphEffMon->SetPoint(i-1, HistEn->GetBinCenter(i), EffMon);
        GraphEffMon->SetPointError(i-1, -4.83e-3+0.0416*HistEn->GetBinCenter(i)+3.82e-3*pow(HistEn->GetBinCenter(i),2), sqrt(ErrEffMon));

        GraphEffMonPlot->SetPoint(i-1, HistEn->GetBinCenter(i), EffMon * 100.); 
        GraphEffMonMin->SetPoint(i-1, HistEn->GetBinCenter(i), (EffMon - (sqrt(ErrEffMon)/2.) - (0.025*EffMon)) * 100.);
        GraphEffMonMax->SetPoint(i-1, HistEn->GetBinCenter(i), (EffMon + (sqrt(ErrEffMon)/2.) + (0.025*EffMon))* 100.);
        
        x.push_back(HistEn->GetBinCenter(i));
        y.push_back(EffMon * 100.);
        ymin.push_back((EffMon - (sqrt(ErrEffMon)/2.) - (0.025*EffMon)) * 100.);
        ymax.push_back((EffMon + (sqrt(ErrEffMon)/2.) + (0.025*EffMon)) * 100.);
    }

    for(int i = 1; i <= Nbins; i++)
    {
        GraphEffMonShade->SetPoint(i-1, x[i-1], ymax[i-1]);       
        GraphEffMonShade->SetPoint(Nbins+i-1, x[Nbins-i], ymin[Nbins-i]);       
    }

    GraphEffMonShade->SetFillColorAlpha(kRed, 0.5);
    GraphEffMonShade->Draw("AF");

    GraphEffMonPlot->SetLineColor(2);
    GraphEffMonPlot->SetLineWidth(2);
    GraphEffMonPlot->SetMarkerColor(2);
    GraphEffMonPlot->SetMarkerSize(1);
    GraphEffMonPlot->SetFillStyle(1001);
    GraphEffMonPlot->SetFillColorAlpha(2, 0.5);

    GraphEffMonPlot->Draw("L");
    //GraphEffMon->Draw("SAMEP");

    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);

    legend->AddEntry(GraphEffMonPlot,"MONSTER efficiency","plf");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);

    legend->Draw();

    GraphEffMon->Write("GraphEffMon");

    FileEffMon->Write();
    FileEffMon->Close();
}
