using namespace std;

void PlotHist()
{
	TCanvas* c1 = new TCanvas();

    TFile *fileeff = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *filehist = TFile::Open("~/phd/analysis/monster25/root_files/E_Spectrum_252Cf.root", "READ");

    TGraph* grapheff = (TGraph*)fileeff->Get("Graph");
    TH1D* hist_all = (TH1D*)filehist->Get("hist_all");
    TH1D* hist_bgd = (TH1D*)filehist->Get("hist_bgd");

    hist_all->Add(hist_bgd,-1);

    for(int i=0; i<= hist_all->GetNbinsX(); i++)
    {
        hist_all->SetBinContent(i, hist_all->GetBinContent(i)/grapheff->Eval(hist_all->GetBinCenter(i)));
    }

    hist_all->Draw();

}
