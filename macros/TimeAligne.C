void TimeAligne()
{	
    
    TH1D *hist[39];

    TFile *filein = TFile::Open("/Users/cantacuzene/data/monster25/runs/sorted/RUN29_plastic_trigger.root", "READ");

    TFile *fileout = new TFile("out.root","UPDATE");

    TTree* tree = (TTree*)filein->Get("tcoinc");;
    
    std::vector<double> *var1 = 0;
    std::vector<int> *var2 = 0;
    
    for(int i = 0; i < 39; i++)
    {   
        
        TString namehist = Form("hist%d",i);

        hist[i] = new TH1D(namehist, "Histogram", 2000, -100, 100);
       
        tree->SetBranchAddress("MonsterPlastic_tDiff", &var1);
        tree->SetBranchAddress("MonsterPlastic_Id",&var2);
        
        int fEntries = tree->GetEntries();
        
        for(int j = 0; j < fEntries; j++)
        {
            tree->GetEntry(j);

            for(ULong_t k = 0; k < var1->size(); ++k)
            { 
                if(var2->at(k) == i+3) hist[i]->Fill(var1->at(k));
            }
        }

        cout << i << endl;
    }
    
    fileout->Write();
}