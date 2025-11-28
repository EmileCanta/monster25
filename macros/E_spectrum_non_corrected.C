#include "../include/E_spectrum.h"

Double_t GetDistanceCell(Int_t index)
{
    Double_t distance = (distances_cell[index-3]*0.01)+0.025;

    return distance;
}

void E_spectrum_non_corrected()
{	
    TFile *filein = TFile::Open("~/phd/data/monster25/sorted/82Ga/AllRuns.root", "READ");

    TFile *fileeff = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");

    TFile *cutfile = TFile::Open("~/phd/analysis/monster25/root_files/cutfile.root","READ");
    
    TCutG *cut = new TCutG();

    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/E_Spectrum_82Ga_1125.root","RECREATE");

    TTree* tree = (TTree*)filein->Get("tcoinc");

    TGraph* graph = (TGraph*)fileeff->Get("Graph");

    TString name;
    
    std::vector<double> *tdiff = 0;
    std::vector<double> *nrj2 = 0;
    std::vector<double> *nrj3 = 0;
    std::vector<double> *label = 0;
    Int_t mult;

    tree->SetBranchAddress("MonsterPlastic_tDiff", &tdiff);
    tree->SetBranchAddress("MonsterPlastic_Q3Cond",&nrj3);
    tree->SetBranchAddress("MonsterPlastic_Q2Cond",&nrj2);
    tree->SetBranchAddress("MonsterPlastic_Id",&label);
    tree->SetBranchAddress("MonsterPlastic_Mult",&mult);

    TH1D *hist_bgd = new TH1D("hist_bgd", "hist_bgd", 1500, 0, 15);
    TH1D *hist_all = new TH1D("hist_all", "hist_all", 1500, 0, 15);

    TH1D *hist_tof = new TH1D("hist_tof", "hist_tof", 2000, -1000, 1000);

    int fEntries = tree->GetEntries();

    for(int j = 0; j < fEntries; j++)
    {
        tree->GetEntry(j);

        for(ULong_t k = 0; k < tdiff->size(); ++k)
        { 
            for(int l = 3; l <= 41; l++)
            {   
                name = Form("cut%d",l);

                if(label->at(k) == l)
                {
                    cut = (TCutG*)cutfile->Get(name);

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && tdiff->at(k) >= -300 && tdiff->at(k) <=-25 && mult <= 3)
                    {
                        hist_bgd->Fill(5.21980556e3*(GetDistanceCell(label->at(k)))*(GetDistanceCell(label->at(k)))/(tdiff->at(k)*tdiff->at(k)));
                    }

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && tdiff->at(k) >=25 && tdiff->at(k) <=300 && mult <= 3)
                    {
                        hist_all->Fill(5.21980556e3*(GetDistanceCell(label->at(k)))*(GetDistanceCell(label->at(k)))/(tdiff->at(k)*tdiff->at(k)));
                    }

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && mult <= 3)
                    {
                        hist_tof->Fill(tdiff->at(k));
                    }

                    delete cut;
                }

            }


        }

        std::cout << std::setprecision(3) << std::setw(5) << (100.*j/fEntries) << " %\r";
    }

    /*for(int i=1; i<=1500; i++)
	{
		hist2->SetBinContent(i,hist2->GetBinContent(i)/(graph->Eval(hist2->GetBinCenter(i))));
	}*/
    
    hist_all->Draw();
    
    fileout->Write();
}
