#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>

void distances()
{

    const int n = 39;
    double x[n]        = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    double y[n]        = {151,152,155,157,160,162,163,152,152,154,156,157,158,161,152,153,153,154,155,157,158,152,152,154,154,155,153,154,154,155,156,152,153,154,156,157,154,156,158};
    double ex[n]       = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // erreur sur x
    double ey[n]       = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // erreur sur y

    TCanvas *c1 = new TCanvas("c1", "Graph avec erreurs", 800, 600);
    TGraphErrors *gr = new TGraphErrors(n, x, y, ex, ey);
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->Draw("AP");

    auto fitFunc = new TF1("fitFunc","[0]",0,39);
    fitFunc->SetLineColor(kRed);
    gr->Fit(fitFunc, "R"); // "R" limite le fit à la plage spécifiée (0 à 9)

    c1->Update();
}
