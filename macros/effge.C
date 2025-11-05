#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>

void effge()
{

    const int n = 11;
    double x[n] = {121.8, 244.7, 344.3, 411.1, 444.0, 778.9, 867.4, 964.1, 1085.9, 1112.1, 1408};
    double y1[n] = {2.13E-02,1.50E-02,1.27E-02,1.10E-02,1.02E-02,7.14E-03,6.33E-03,6.24E-03,6.00E-03,5.75E-03,4.81E-03};
    double y2[n] = {1.61E-02,1.22E-02,1.07E-02,9.82E-03,9.03E-03,6.36E-03,5.79E-03,5.59E-03,5.32E-03,5.24E-03,4.47E-03};
    double y3[n] = {2.00E-02,1.82E-02,1.58E-02,1.34E-02,1.29E-02,9.79E-03,8.61E-03,8.83E-03,8.61E-03,8.06E-03,7.01E-03};
    double y4[n] = {5.501E-02,4.388E-02,3.785E-02,3.288E-02,3.094E-02,2.251E-02,1.998E-02,1.998E-02,1.932E-02,1.841E-02,1.574E-02};
    double ex[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double ey1[n] = {1.28E-03,8.99E-04,7.62E-04,6.60E-04,6.13E-04,4.29E-04,3.80E-04,3.74E-04,3.60E-04,3.45E-04,2.89E-04};
    double ey2[n] = {9.68E-04,7.32E-04,6.42E-04,5.89E-04,5.42E-04,3.82E-04,3.47E-04,3.35E-04,3.19E-04,3.15E-04,2.68E-04};
    double ey3[n] = {1.20E-03,1.09E-03,9.49E-04,8.07E-04,7.75E-04,5.87E-04,5.17E-04,5.30E-04,5.17E-04,4.84E-04,4.21E-04};
    double ey4[n] = {3.30E-03,2.63E-03,2.27E-03,1.97E-03,1.86E-03,1.35E-03,1.20E-03,1.20E-03,1.16E-03,1.10E-03,9.45E-04};

    TCanvas *c1 = new TCanvas("c1", "Graph avec erreurs", 800, 600);

    TMultiGraph *mg = new TMultiGraph();

    TGraphErrors *gr1 = new TGraphErrors(n, x, y1, ex, ey1);
    TGraphErrors *gr2 = new TGraphErrors(n, x, y2, ex, ey2);
    TGraphErrors *gr3 = new TGraphErrors(n, x, y3, ex, ey3);
    TGraphErrors *gr4 = new TGraphErrors(n, x, y4, ex, ey4);

    gr1->SetMarkerStyle(21);
    gr1->SetMarkerColor(kBlue);
    gr1->SetLineColor(kBlue);
    gr1->SetMarkerSize(2);

    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(kRed);
    gr2->SetLineColor(kRed);
    gr2->SetMarkerSize(2);

    gr3->SetMarkerStyle(21);
    gr3->SetMarkerColor(kBlack);
    gr3->SetLineColor(kBlack);
    gr3->SetMarkerSize(2);

    gr4->SetMarkerStyle(21);
    gr4->SetMarkerColor(kGreen);
    gr4->SetLineColor(kGreen);
    gr4->SetMarkerSize(2);

    mg->Add(gr3);
    mg->Add(gr1);
    mg->Add(gr2);
    //mg->Add(gr4);

    mg->Draw("AP");

    auto fitFunc1 = new TF1("fitFunc1", "[0]*TMath::Power(x, [1])", 0., 1408.);
    auto fitFunc2 = new TF1("fitFunc2", "[0]*TMath::Power(x, [1])", 0., 1408.);
    auto fitFunc3 = new TF1("fitFunc3", "[0]*TMath::Power(x, [1])", 244.7, 1408.);
    fitFunc1->SetLineColor(kBlue);
    fitFunc2->SetLineColor(kRed);
    fitFunc3->SetLineColor(kBlack);

    fitFunc1->SetLineWidth(7);
    fitFunc2->SetLineWidth(7);
    fitFunc3->SetLineWidth(7);

    gr1->Fit(fitFunc1,"R");
    gr2->Fit(fitFunc2,"R");
    gr3->Fit(fitFunc3,"R");

    c1->Update();
}
