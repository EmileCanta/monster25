using namespace std;

void PlotHist()
{
	TCanvas* c1 = new TCanvas();

    TFile *fileeff = TFile::Open("/Users/cantacuzene/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *filehist = TFile::Open("/Users/cantacuzene/phd/analysis/monster25/root_files/en_82ga_best.root", "READ");

    TGraph* grapheff = (TGraph*)fileeff->Get("Graph");
    TH1D* hist = (TH1D*)filehist->Get("hist_all");

    for(int i=23; i<= 233;i++)
    {
        hist->SetBinContent(i, hist->GetBinContent(i)/grapheff->Eval(hist->GetBinCenter(i)));
    }

    hist->Draw();

}
