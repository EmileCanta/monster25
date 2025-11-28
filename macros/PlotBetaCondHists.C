void PlotBetaCondHists()
{
	TCanvas* c = new TCanvas();
    
    TPad* pad3 = new TPad("pad3","pad3",0,0.1,1,0.4);
 	TPad* pad2 = new TPad("pad2","pad2",0,0.4,1,0.7);
 	TPad* pad1 = new TPad("pad1","pad1",0,0.7,1,1.);
    
    //pad1->SetFillColor(kBlack);
    //pad2->SetFillColor(kBlack);
    //pad3->SetFillColor(kBlack);

    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    
    TFile *filesingle = TFile::Open("~/phd/analysis/monster25/root_files/Single.root", "READ");
    TFile *filebetacond = TFile::Open("~/phd/analysis/monster25/root_files/BetaCondGeHists.root", "READ");

    TH1D* histsingle = (TH1D*)filesingle->Get("hist3");
    TH1D* histbetacond = (TH1D*)filebetacond->Get("hist3");

    TH1D* histsingle1 = (TH1D*)histsingle->Clone("hs1");
    TH1D* histbetacond1 = (TH1D*)histbetacond->Clone("hb1");
    
    TH1D* histsingle2 = (TH1D*)histsingle->Clone("hs2");
    TH1D* histbetacond2 = (TH1D*)histbetacond->Clone("hb2");
    
    TH1D* histsingle3 = (TH1D*)histsingle->Clone("hs3");
    TH1D* histbetacond3 = (TH1D*)histbetacond->Clone("hb3");

    pad1->cd(); 
    
    histsingle1->Draw();
    histbetacond1->Draw("same");

    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);

    legend->AddEntry(histsingle1,"Single","f");
    legend->AddEntry(histbetacond1,"Beta-conditionned","f");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);

    legend->Draw();
    
    pad2->cd(); 
    
    histsingle2->Draw();
    histbetacond2->Draw("same");
    
    pad3->cd(); 
    
    histsingle3->Draw();
    histbetacond3->Draw("same");
}
