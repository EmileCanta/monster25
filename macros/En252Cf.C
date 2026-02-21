void En252Cf()
{	
    TFile *fileeffmon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *filein = TFile::Open("~/phd/analysis/monster25/root_files/ToF252Cf.root", "READ");

    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/En252Cf.root","RECREATE");

    TH1D* hist_tof_sub = (TH1D*)filein->Get("hist_tof_sub");

    double massn = 1.67492750056e-27;
    double joultoMeV = 6.241509343260e12;
    double d = 1.575;

    int NsmearPerBin = 500;

    int Ntoys = 5000;

    int NbinsE = 500;
    double MaxE = 20.;
    double MinE = 0.;

    double SmearingVar = 2.35*0.8e-9;

    TRandom3 random(0);

    TH1D* hSum = new TH1D("hSum", "hSum", NbinsE, MinE, MaxE);

    hSum->Sumw2();

    TH1D* hSum2 = (TH1D*)hSum->Clone("hSum2");

    for(int it = 0; it < Ntoys; ++it)
    {
        TH1D hToy("hToy", "hToy", NbinsE, MinE, MaxE);
        hToy.Sumw2();

        for(int ib=1; ib<=hist_tof_sub->GetNbinsX(); ++ib)
        {
            double Ni = hist_tof_sub->GetBinContent(ib);
            double errNi = hist_tof_sub->GetBinError(ib);

            if (Ni <= 0) continue;

            double Ni_toy;

            //if(errNi > 0.0) { Ni_toy = random.Gaus(Ni, errNi); }

            //else { Ni_toy = random.Poisson(Ni); }
            
            Ni_toy = random.Poisson(Ni);

            //if(Ni_toy <= 0) continue;

            double t_center = hist_tof_sub->GetXaxis()->GetBinCenter(ib);

            double weight_per_smear = Ni_toy / double(NsmearPerBin);

            for(int ks = 0; ks<NsmearPerBin; ++ks)
            {
                double t_sample = gRandom->Gaus(t_center*1e-9, SmearingVar);

                if(t_sample <= 0) continue;

                double K = joultoMeV * 0.5 * massn * d * d / (t_sample * t_sample);

                hToy.Fill(K, weight_per_smear);
            }
        }

        for(int kb = 1; kb <= NbinsE; ++kb)
        {
            double v = hToy.GetBinContent(kb);
            hSum->AddBinContent(kb, v);
            hSum2->AddBinContent(kb, v*v);
        }
    }

    TH1D* hist_E = (TH1D*)hSum->Clone("hist_E");
    hist_E->Reset();

    for(int kb = 1; kb <= NbinsE; ++kb)
    {
        double mean = hSum->GetBinContent(kb) / double(Ntoys);
        double mean2 = hSum2->GetBinContent(kb) / double(Ntoys);
        double var = mean2 - mean * mean;

        if (var < 0 && var > -1e-18) var = 0;
        double rms = (var > 0) ? sqrt(var) : 0.0;

        hist_E->SetBinContent(kb, mean);
        hist_E->SetBinError(kb, rms);
    }

    delete hSum2;
    delete hSum;

    fileout->Write();
    fileout->Close();
}
