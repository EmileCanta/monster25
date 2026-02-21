void effplastic()
{
    const int n = 14;
    double x[n]        = {8.6637, 7.952523, 7.7681, 7.3605, 6.80839, 6.11376, 5.53468, 4.89085, 4.38696, 3.38696, 2.38696, 1.38696, 0.88696, 0.38696};
    double y[n]        = {0.6657753, 0.6697565, 0.6672611, 0.6618724, 0.6562806, 0.6362569, 0.6287861, 0.6043276, 0.5592845, 0.5141187, 0.4131114, 0.2651557, 0.1461876, 0.0343465};
    double ex[n]       ={0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // erreur sur x
    double ey[n]       = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; // erreur sur y
    
    reverse(begin(x), end(x));
    reverse(begin(y), end(y));
    reverse(begin(ex), end(ex));
    reverse(begin(ey), end(ey));

    //double eytest[n];
    //fill(begin(eytest), end(eytest), 0);

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);

    TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);

    TGraph *gr0 = new TGraph();
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    TGraph *GraphShade = new TGraph();

    for(int i = 0; i<=13; i++)
    {
        gr0->AddPoint(x[i], (y[i] * 100.));
        gr1->AddPoint(x[i], (y[i] + (ey[i]/2)) * 100.);
        gr2->AddPoint(x[i], (y[i] - (ey[i]/2)) * 100.);
    }

    for(int i = 0; i < n; i++)
    {
        GraphShade->SetPoint(i, x[i], (y[i] + (ey[i]/2)) * 100.);
        GraphShade->SetPoint(n+i, x[n-i-1], (y[n-i-1] - (ey[n-i-1]/2)) * 100.);
    }

    GraphShade->SetFillColorAlpha(kBlue, 0.5);
    //GraphShade->Draw("AF");

    gr0->SetLineColor(kBlue);
    gr0->SetLineWidth(2);
    gr0->SetMarkerColor(kBlue);
    gr0->SetMarkerSize(1);
    gr0->SetFillStyle(1001);
    gr0->SetFillColorAlpha(kBlue, 0.5);
    //gr0->Draw("CP");
    
    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);

    legend->AddEntry(gr0,"Plastic efficiency","plf");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    //legend->Draw();
    gr->Draw();
}
