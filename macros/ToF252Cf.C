void ToF252Cf()
{	
    TFile *filein = TFile::Open("~/phd/data/monster25/sorted/cfruns/RUN111.root", "READ");
    
    TFile *cutfile = TFile::Open("~/phd/analysis/monster25/root_files/cutfile.root","READ");
    
    TCutG *cut = new TCutG();
    
    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/ToF252Cf.root","RECREATE");
    
    TTree* tree = (TTree*)filein->Get("tcoinc");
    
    TString name;
    
    std::vector<double> *tdiff = 0;
    std::vector<double> *nrj2 = 0;
    std::vector<double> *nrj3 = 0;
    std::vector<double> *label = 0;
    Int_t mult;
    
    double windowback = 25.; //Corresponds to Eneut = 20.7 MeV
    double windowfront = 216.; //Corresponds to Eneut = 280 kev (threshold)
    int nbins_tof = 48; //Bin width equal to 4 ns

    tree->SetBranchAddress("MonsterPlastic_tDiff", &tdiff);
    tree->SetBranchAddress("MonsterPlastic_Q3Cond",&nrj3);
    tree->SetBranchAddress("MonsterPlastic_Q2Cond",&nrj2);
    tree->SetBranchAddress("MonsterPlastic_Id",&label);
    tree->SetBranchAddress("MonsterPlastic_Mult",&mult);

    TH1D *hist_tof_all = new TH1D("hist_tof_all", "hist_tof_all", nbins_tof, windowback, windowfront);
    TH1D *hist_tof_bgd = new TH1D("hist_tof_bgd", "hist_tof_bgd", nbins_tof, windowback, windowfront);
    
    hist_tof_all->Sumw2();
    hist_tof_bgd->Sumw2();

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

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && tdiff->at(k) >= -windowfront && tdiff->at(k) <=-windowback && mult <= 3)
                    {
                        hist_tof_bgd->Fill(abs(tdiff->at(k)));
                    }

                    if(cut->IsInside(nrj2->at(k),nrj3->at(k)/nrj2->at(k)) && tdiff->at(k) >= windowback && tdiff->at(k) <= windowfront && mult <= 3)
                    {
                        hist_tof_all->Fill(tdiff->at(k));
                    } 

                    delete cut;
                }
            }
        }
        
        std::cout << std::setprecision(3) << std::setw(5) << (100.*j/fEntries) << " %\r";
    }

    TH1D* hist_tof_sub = (TH1D*)hist_tof_all->Clone("hist_tof_sub");

    hist_tof_sub->Add(hist_tof_bgd,-1);
    
    fileout->Write();
}
