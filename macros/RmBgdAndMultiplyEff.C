using namespace std;

void RmBgdAndMultiplyEff()
{
	TCanvas* c1 = new TCanvas();
    
    //TPad* pad1 = new TPad("pad1","pad1",0,0,1,0.3);
 	//TPad* pad2 = new TPad("pad2","pad2",0,0.3,1,1);
    
    //pad1->Draw();
    //pad2->Draw();
    
    TFile *fileeff = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *filehist = TFile::Open("~/phd/analysis/monster25/root_files/E_Spectrum_82Ga_1125.root", "READ");
    TFile *filemanhart = TFile::Open("~/phd/analysis/monster25/root_files/manhart.root", "READ");
    TFile *filebgd_manual = TFile::Open("~/phd/analysis/monster25/root_files/BgdModel.root", "READ");
    TFile *fileibu = TFile::Open("/home/emile/phd/analysis/mcnp_tetra/root_files/tests_regina/gaibu.root", "READ");

    TGraph* grapheff = (TGraph*)fileeff->Get("Graph");
    TGraph* graphmanhart = (TGraph*)filemanhart->Get("g1");
    TH1D* hist_all = (TH1D*)filehist->Get("hist_all");
    TH1D* hist_bgd = (TH1D*)filehist->Get("hist_bgd");
    TH1D* hist_bgd_manual = (TH1D*)filebgd_manual->Get("hist_model_bgd");
    TH1D* hist_ibu = (TH1D*)fileibu->Get("proba");
    
    hist_all->Add(hist_bgd_manual,-1);
   
    hist_all->Rebin(10);

    hist_ibu->Scale(1000000.);

    //graphmanhart->Scale(952300.);
    
    TH1D* hist_resi = (TH1D*)hist_all->Clone("hist_resi");
    
    hist_resi->Reset();
    
    for(int i=0; i<= hist_all->GetNbinsX(); i++)
    {
        hist_all->SetBinContent(i, hist_all->GetBinContent(i)/grapheff->Eval(hist_all->GetBinCenter(i)));
    }
    
    for(int i=0; i<= hist_all->GetNbinsX(); i++)
    {
        hist_resi->SetBinContent(i, (graphmanhart->Eval(hist_all->GetBinCenter(i))-hist_all->GetBinContent(i))/graphmanhart->Eval(hist_all->GetBinCenter(i)));    
    }

    //pad2->cd(); 
    
    gStyle->SetErrorX(0.);
    
    hist_ibu->Draw("hist");
    //graphmanhart->Draw("same");

    hist_all->Draw("sameshist");

    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);

    legend->AddEntry(hist_all,"Data","f");
    //legend->AddEntry(graphmanhart,"Mannhart evaluation","l");
    legend->AddEntry(hist_ibu,"Iterative Bayesian unfolding method","f");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);

    TLine* l = new TLine(0.3,0,4.8,0);

    legend->Draw();

    //pad1->cd();

    //gStyle->SetErrorX(0.);

    //hist_resi->Draw("E0 P");
    
    //l->Draw("same");
}
