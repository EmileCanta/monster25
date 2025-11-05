#include "../include/E_spectrum.h"

Double_t f1(Double_t *x, Double_t *par)
{
    return (19./(x[0] * x[0]))+(-3.7*x[0]+18.);
}

void PlotBgdModel()
{	
    TFile *fileout = new TFile("/Users/cantacuzene/phd/analysis/monster25/root_files/BgdModel.root","RECREATE");

    TF1* function = new TF1("f1", f1, 0, 15);

    TH1D *hist_model_bgd = new TH1D("hist_model_bgd", "hist_model_bgd", 1500, 0, 15);

    for(int i=1; i<=1500; i++)
	{
		hist_model_bgd->SetBinContent(i,function->Eval((float)i/100.));
	}
    
    hist_model_bgd->Draw("hist");

    fileout->Write();
}