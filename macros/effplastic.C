void effplastic()
{
    const int n = 14;
    double x[n]        = {8.6637, 7.952523, 7.7681, 7.3605, 6.80839, 6.11376, 5.53468, 4.89085, 4.38696, 3.38696, 2.38696, 1.38696, 0.88696, 0.38696};
    double y[n]        = {0.6657753, 0.6697565, 0.6672611, 0.6618724, 0.6562806, 0.6362569, 0.6287861, 0.6043276, 0.5592845, 0.5141187, 0.4131114, 0.2651557, 0.1461876, 0.0343465};
    double ex[n]       ={0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // erreur sur x
    double ey[n]       = {0.06657753, 0.06697565, 0.06672611, 0.06618724, 0.06562806, 0.06362569, 0.06287861, 0.06043276, 0.05592845, 0.05141187, 0.04131114, 0.02651557, 0.01461876, 0.00343465}; // erreur sur y
    
    reverse(begin(x), end(x));
    reverse(begin(y), end(y));
    reverse(begin(ex), end(ex));
    reverse(begin(ey), end(ey));

    //double eytest[n];
    //fill(begin(eytest), end(eytest), 0);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);

    TGraph *gr0 = new TGraph();
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();

    for(int i = 0; i<=13; i++)
    {
        gr0->AddPoint(x[i], (y[i] * 100.));
        gr1->AddPoint(x[i], (y[i] + (ey[i]/2)) * 100.);
        gr2->AddPoint(x[i], (y[i] - (ey[i]/2)) * 100.);
    }

    //gr->Draw("ACP");
    
    gr0->Draw("APL");
    gr1->Draw("LSAME");
    gr2->Draw("LSAME");
}
