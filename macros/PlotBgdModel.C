#include "../include/E_spectrum.h"

Double_t f1(Double_t *x, Double_t *par)
{
    return (80.278/(x[0] * x[0]))+(-0.666354*x[0])+12.8215;
}

void PlotBgdModel()
{	
    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/BgdModel.root","RECREATE");

    TF1* function = new TF1("f1", f1, 0, 20);

    TH1D *hist_model_bgd = new TH1D("hist_model_bgd", "hist_model_bgd", 500, 0, 20);

    for(int i=1; i<=500; i++)
	{
		hist_model_bgd->SetBinContent(i,function->Eval((float)i*0.04));
	}
    
    hist_model_bgd->Draw("hist");

    fileout->Write();
}
