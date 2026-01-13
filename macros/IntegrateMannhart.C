void IntegrateMannhart()
{
    TFile* FileMann = TFile::Open("~/phd/analysis/monster25/root_files/manhart.root", "READ");
    TGraph* GraphMann = (TGraph*)FileMann->Get("mannhart");
    
    TH1D* HistMann = new TH1D("HistMann", "HistMann", 500, 0, 20);

    for(int i = 1; i <= 500; i++)
    {
        HistMann->SetBinContent(i, GraphMann->Eval(HistMann->GetBinCenter(i) * 1.e6)  * (20./500.) * 1.e6);
    }
        
    //HistMann->Scale(1./HistMann->Integral());
    HistMann->Draw("hist");
}
