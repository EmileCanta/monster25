void Analyze252Cf()
{	
    TFile *fileeffmon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *filein = TFile::Open("~/phd/data/monster25/sorted/cfruns/RUN111.root", "READ");
    
    TGraph* grapheffmon = (TGraph*)fileeffmon->Get("Graph");

    TFile *cutfile = TFile::Open("~/phd/analysis/monster25/root_files/cutfile.root","READ");
    
    TCutG *cut = new TCutG();
    
    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/252Cf_ana.root","RECREATE");
    
    TTree* tree = (TTree*)filein->Get("tcoinc");
    
    TString name;
    
    std::vector<double> *tdiff = 0;
    std::vector<double> *nrj2 = 0;
    std::vector<double> *nrj3 = 0;
    std::vector<double> *label = 0;
    Int_t mult;
    
    double windowback = 25.; //Corresponds to Eneut = 20.7 MeV
    double windowfront = 216.; //Corresponds to Eneut = 280 kev (threshold)
    int nbins_tof = 174; //Bin width equal to sigma =  1.1 ns

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
    
    double massn = 1.67492750056e-27;
    double joultoMeV = 6.241509343260e12;
    double d = 1.575;

    int NsmearPerBin = 50;

    int Ntoys = 2000;
    
    TRandom3 random(0);
    
    TH1D* hSum = new TH1D("hSum", "hSum", 500, 0, 20);
    hSum->Sumw2();

    TH1D* hSum2 = (TH1D*)hSum->Clone("hSum2");

    for(int it = 0; it < Ntoys; ++it)
    {
        TH1D hToy("hToy", "hToy", 500, 0, 20);
        hToy.Sumw2();

        for(int ib=1; ib<=hist_tof_sub->GetNbinsX(); ++ib)
        {
            double Ni = hist_tof_sub->GetBinContent(ib);
            double errNi = hist_tof_sub->GetBinError(ib);
        
            if (Ni <= 0) continue;
            
            double Ni_toy;

            if(errNi > 0.0) { Ni_toy = random.Gaus(Ni, errNi); }

            else { Ni_toy = random.Poisson(Ni); }
                
            if(Ni_toy <= 0) continue;

            double t_center = hist_tof_sub->GetXaxis()->GetBinCenter(ib);
        
            double weight_per_smear = Ni_toy / double(NsmearPerBin);

            for(int ks = 0; ks<NsmearPerBin; ++ks)
            {
                double t_sample = gRandom->Gaus(t_center*1e-9, 1.1e-9/2.);

                if(t_sample <= 0) continue;

                double K = joultoMeV * 0.5 * massn * d * d / (t_sample * t_sample);

                hToy.Fill(K, weight_per_smear);
            }
        }

        for(int kb = 1; kb <= 500; ++kb)
        {
            double v = hToy.GetBinContent(kb);
            hSum->AddBinContent(kb, v);
            hSum2->AddBinContent(kb, v*v);
        }
    }
    
    TH1D* hist_E = (TH1D*)hSum->Clone("hist_E");
    hist_E->Reset();

    for(int kb = 1; kb <= 500; ++kb)
    {
        double mean = hSum->GetBinContent(kb) / double(Ntoys);
        double mean2 = hSum2->GetBinContent(kb) / double(Ntoys);
        double var = mean2 - mean * mean;

        if (var < 0 && var > -1e-18) var = 0;
        double rms = (var > 0) ? sqrt(var) : 0.0;

        hist_E->SetBinContent(kb, mean);
        hist_E->SetBinError(kb, rms);
    }

    TRandom3* rnd = new TRandom3(0);
    
    Int_t NtoysMon = 2000;

    TH1D* hist_E_corrected_mon = (TH1D*)hist_E->Clone("hist_E_corrected_mon");
    hist_E_corrected_mon->Reset();
     
    const int nBins = hist_E->GetNbinsX();
    const int nPoints = grapheffmon->GetN();

    std::vector<double> sum(nBins+1, 0.0);
    std::vector<double> sum2(nBins+1, 0.0);
    
    std::vector<double> sumInv(nBins+1, 0.0);
    std::vector<double> sumInv2(nBins+1, 0.0);

    std::vector<double> gx(nPoints), gy(nPoints), gey(nPoints);
    
    for(int i=0;i<nPoints;++i)
    {
        grapheffmon->GetPoint(i, gx[i], gy[i]);
        
        gey[i] = grapheffmon->GetErrorY(i);
    }

    for(int itoy=0; itoy<NtoysMon; ++itoy)
    {
        std::vector<double> toyY(nPoints);
        
        for(int ip=0; ip<nPoints; ++ip)
        {
            double sigma = gey[ip];
            toyY[ip] = (sigma > 0.0) ? gy[ip] + rnd->Gaus(0.0, sigma) : gy[ip];
        }

        for(int ib=1; ib<=nBins; ++ib)
        {
            double xBin = hist_E->GetXaxis()->GetBinCenter(ib);
            double scale = 0.0;
            
            if(xBin <= gx.front())
            {
                scale = toyY.front();
            }
            
            else if(xBin >= gx.back())
            {
                scale = toyY.back();
            }
            
            else
            {
                int k = 0;
                
                for(int j=0;j<nPoints-1;++j)
                
                {
                    if(xBin >= gx[j] && xBin < gx[j+1]) { k = j; break; }
                }
                
                double x0 = gx[k], x1 = gx[k+1];
                double y0 = toyY[k], y1 = toyY[k+1];
                double t = (xBin - x0) / (x1 - x0);
                
                scale = y0 + t*(y1 - y0);
            }

            double inv = 1. / scale;
            double C = hist_E->GetBinContent(ib);
            double val = C * inv;

            sum[ib]  += val;
            sum2[ib] += val*val;
            sumInv[ib]  += inv;
            sumInv2[ib] += inv*inv;
        }
    }

    for(int ib=1; ib<=nBins; ++ib)
    {
       double mean = sum[ib] / double(NtoysMon);
       double mean2 = sum2[ib] / double(NtoysMon);
       double var = mean2 - mean*mean;
       
       if(var < 0 && var > -1e-18) var = 0;
       
       double rms = (var>0) ? sqrt(var) : 0.0;

       double sigmaC = hist_E->GetBinError(ib);

       double meanInv2 = sumInv2[ib] / double(NtoysMon);   
    
       double sigma_from_original = sigmaC * sqrt(meanInv2);
       double total_err = sqrt(rms*rms + sigma_from_original*sigma_from_original);
       
       hist_E_corrected_mon->SetBinContent(ib, mean);
       hist_E_corrected_mon->SetBinError(ib, total_err);
    }
    
    delete hSum2;
    delete hSum;
    
    fileout->Write();
}
