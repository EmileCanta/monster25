#include "../include/E_spectrum.h"

Double_t GetDistanceCell(Int_t index)
{
    Double_t distance = (distances_cell[index-3]*0.01)+0.025;

    return distance;
}

void MakeEspectrumNonCorrected()
{	
    TFile *filein = TFile::Open("~/phd/data/monster25/sorted/82Ga/AllRuns.root", "READ");

    TFile *cutfile = TFile::Open("~/phd/analysis/monster25/root_files/cutfile.root","READ");
    
    TCutG *cut = new TCutG();

    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/TOF_82Ga_0312.root","RECREATE");

    TTree* tree = (TTree*)filein->Get("tcoinc");

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

    TH1D *hist_bgd = new TH1D("hist_bgd", "hist_bgd", 500, 0, 20);
    TH1D *hist_all = new TH1D("hist_all", "hist_all", 500, 0, 20);

    TH1D *hist_tof_all = new TH1D("hist_tof_all", "hist_tof_all", 625, -1000, 1000);
    TH1D *hist_tof_bgd = new TH1D("hist_tof_bgd", "hist_tof_bgd", 625, -1000, 1000);
    TH1D *hist_tof_forvarbins = new TH1D("hist_tof_forvarbins", "hist_tof_forvarbins", 305, 25, 1000);
    
    const int Ntof = hist_tof_forvarbins->GetNbinsX();

    vector<double> keEdges;

    double nSigma = 2.0;

    double massn = 1.67492750056e-27;
    double joultoMeV = 6.241509343260e12;
    double d = 1.575;

    for(int i = 1; i <= Ntof; ++i)
    {
        double t = hist_tof_forvarbins->GetBinCenter(i);
        
        double K = joultoMeV * massn * 0.5 * pow(d,2) * pow((t*1.e-9),-2);

        double sigmaK = joultoMeV * massn * pow(d,2) * pow((t*1.e-9),-3) * 1.6e-9;

        double right = K + nSigma*sigmaK;

        keEdges.push_back(right);
    }

    reverse(keEdges.begin(), keEdges.end());

    int nBinsKE = keEdges.size() -1;
    TH1D* hist_varbins_all = new TH1D("hist_varbins_all", "hist_varbins_all", nBinsKE, keEdges.data());
    TH1D* hist_varbins_bgd = new TH1D("hist_varbins_bgd", "hist_varbins_bgd", nBinsKE, keEdges.data());
    
    hist_bgd->Sumw2();
    hist_all->Sumw2();
    hist_varbins_all->Sumw2();
    hist_varbins_bgd->Sumw2();

    int fEntries = tree->GetEntries();

    int N_sample = 1000;

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

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && tdiff->at(k) >= -500 && tdiff->at(k) <=-50 && mult <= 3)
                    {
                        //for(int m = 0; m <= N_sample; ++m)
                        {
                            double t_smear = gRandom->Gaus(tdiff->at(k)*1.e-9, 1.6*1.e-9);

                            //hist_bgd->Fill(joultoMeV*massn*0.5*(GetDistanceCell(label->at(k)))*(GetDistanceCell(label->at(k)))/(t_smear*t_smear), 1./(double)N_sample);
                            //hist_varbins_bgd->Fill(joultoMeV*massn*0.5*(GetDistanceCell(label->at(k)))*(GetDistanceCell(label->at(k)))/(t_smear*t_smear), 1./(double)N_sample);
                        }
                        
                        hist_tof_bgd->Fill(tdiff->at(k));
                    }

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && tdiff->at(k) >=50 && tdiff->at(k) <=500 && mult <= 3)
                    {
                        //for(int m = 0; m <= N_sample; ++m)
                        {
                            double t_smear = gRandom->Gaus(tdiff->at(k)*1.e-9, 1.6*1.e-9);

                            //hist_all->Fill(joultoMeV*massn*0.5*(GetDistanceCell(label->at(k)))*(GetDistanceCell(label->at(k)))/(t_smear*t_smear), 1./(double)N_sample);
                            //hist_varbins_all->Fill(joultoMeV*massn*0.5*(GetDistanceCell(label->at(k)))*(GetDistanceCell(label->at(k)))/(t_smear*t_smear), 1./(double)N_sample);
                        }

                        hist_tof_all->Fill(tdiff->at(k));
                    }

                    delete cut;
                }

            }


        }

        std::cout << std::setprecision(3) << std::setw(5) << (100.*j/fEntries) << " %\r";
    }
    
    //hist_all->Draw("hist");
    
    fileout->Write();
}
