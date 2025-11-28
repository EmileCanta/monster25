void MakeSingleHists()
{	
    TH1D *hist[4];

    TFile *filein = TFile::Open("/home/emile/phd/data/monster25/sorted/82Ga/AllRuns.root", "READ");

    TFile *fileout = new TFile("../root_files/Single.root","RECREATE");

    TTree* tree = (TTree*)filein->Get("tsingle");;
    
    std::vector<double> *nrj = 0;
    std::vector<double> *label = 0;

    tree->SetBranchAddress("Ge_ESingle",&nrj);
    tree->SetBranchAddress("Ge_IdSingle",&label);

    for(int i =0; i<4; i++) 
    {
        TString namehist = Form("hist%d",i);

        hist[i] = new TH1D(namehist, "Histogram", 10000, 0, 10000);
    }

    int fEntries = tree->GetEntries();
        
    for(int j = 0; j < fEntries; j++)
    {
        tree->GetEntry(j);

        for(ULong_t k = 0; k < label->size(); ++k)
        { 
            if(label->at(k) == 43)
            {
                hist[0]->Fill(nrj->at(k));
                hist[3]->Fill(nrj->at(k));
            }

            if(label->at(k) ==  44)
            {
                hist[1]->Fill(nrj->at(k));
                hist[3]->Fill(nrj->at(k));
            }

            if(label->at(k) ==  45)
            {
                hist[2]->Fill(nrj->at(k));
                hist[3]->Fill(nrj->at(k));
            }
        }

        std::cout << std::setprecision(3) << std::setw(5) << (100.*j/fEntries) << " %\r";
    }
    
    fileout->Write();
}
