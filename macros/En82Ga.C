void En82Ga()
{	
    TFile *fileeffmon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *fileeffpla = TFile::Open("~/phd/analysis/monster25/root_files/effplastic.root", "READ");
    TFile *filein = TFile::Open("~/phd/analysis/monster25/root_files/ToF82Ga.root", "READ");
    
    TGraph* grapheffmon = (TGraph*)fileeffmon->Get("GraphEffMon");
    TGraph* grapheffpla = (TGraph*)fileeffpla->Get("Graph");

    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/En82Ga.root","RECREATE");
   
    TH1D* hist_tof_all = (TH1D*)filein->Get("hist_tof_all");
    TH1D* hist_tof_bgd = (TH1D*)filein->Get("hist_tof_bgd");

    double massn = 1.67492750056e-27;
    double joultoMeV = 6.241509343260e12;
    double d = 1.575;

    double sigma_bkg = 5;
    double bkgPerBin = 207;

    int NsmearPerBin = 500;

    int Ntoys = 5000;
    
    int NbinsE = 500;
    double MaxE = 20.;
    double MinE = 0.; //Corresponds to 40 keV bin width which is the uncertainty of the setup at 1 MeV neutrons
    
    double SmearingVar = 1.2e-9;
    
    TRandom3 random(0);
    
    TH1D* hSum = new TH1D("hSum", "hSum", NbinsE, MinE, MaxE);
    hSum->Sumw2();

    TH1D* hSum2 = (TH1D*)hSum->Clone("hSum2");

    for(int it = 0; it < Ntoys; ++it)
    {
        double delta_bkg = random.Gaus(0.0, sigma_bkg);
        
        TH1D hToy("hToy", "hToy", NbinsE, MinE, MaxE);
        hToy.Sumw2();

        for(int ib=1; ib<=hist_tof_all->GetNbinsX(); ++ib)
        {
            double Ndata_nom = hist_tof_all->GetBinContent(ib);
            double Ndata_toy = random.Poisson(Ndata_nom);

            double Ni_toy = Ndata_toy - (bkgPerBin + delta_bkg);
        
            if (Ni_toy == 0) continue;

            double weight_per_smear = Ni_toy / double(NsmearPerBin);

            for(int ks = 0; ks<NsmearPerBin; ++ks)
            {
                double t_center = hist_tof_all->GetXaxis()->GetBinCenter(ib);
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

    TRandom3* rnd = new TRandom3(0);
    
    Int_t NtoysMon = 5000;

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
        
        TGraph* GraphToy = new TGraph(nPoints, &gx[0], &toyY[0]); 

        TSpline3* SplineToy = new TSpline3("SplineToy", GraphToy);

        for(int ib = hist_E->GetXaxis()->FindBin(0.28); ib<=nBins; ++ib)
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
                /* int k = 0;

                   for(int j=0;j<nPoints-1;++j)

                   {
                   if(xBin >= gx[j] && xBin < gx[j+1]) { k = j; break; }
                   }

                   double x0 = gx[k], x1 = gx[k+1];
                   double y0 = toyY[k], y1 = toyY[k+1];
                   double t = (xBin - x0) / (x1 - x0);

                   scale = y0 + t*(y1 - y0); */
            
                scale = SplineToy->Eval(xBin);

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
    
    Int_t NtoysPla = 5000;
    double Qbeta = 5.290;

    TH1D* hist_E_corrected_MonAndPla = (TH1D*)hist_E_corrected_mon->Clone("hist_E_corrected_MonAndPla");
    hist_E_corrected_MonAndPla->Reset();
     
    const int nPointsPla = grapheffpla->GetN();

    std::vector<double> sumpla(nBins+1, 0.0);
    std::vector<double> sum2pla(nBins+1, 0.0);
    
    std::vector<double> sumInvpla(nBins+1, 0.0);
    std::vector<double> sumInv2pla(nBins+1, 0.0);

    std::vector<double> gxpla(nPointsPla), gypla(nPointsPla), geypla(nPointsPla);
    
    for(int i=0;i<nPointsPla;++i)
    {
        grapheffpla->GetPoint(i, gxpla[i], gypla[i]);
        
        geypla[i] = grapheffpla->GetErrorY(i);
    }

    for(int itoy=0; itoy<NtoysPla; ++itoy)
    {
        std::vector<double> toyY(nPointsPla);
        
        for(int ip=0; ip<nPointsPla; ++ip)
        {
            double sigma = geypla[ip];
            toyY[ip] = (sigma > 0.0) ? gypla[ip] + rnd->Gaus(0.0, sigma) : gypla[ip];
        }
        
        TGraph* GraphToy = new TGraph(nPointsPla, &gxpla[0], &toyY[0]); 

        TSpline3* SplineToy = new TSpline3("SplineToy", GraphToy);
        
        //if(itoy == 0) { SplineToy->Draw(); }

        //else if(itoy % 200 == 0) { SplineToy->Draw("same"); }
        
        for(int ib=1; ib<=nBins; ++ib)
        {
            double xBin = hist_E_corrected_mon->GetXaxis()->GetBinCenter(ib);
            double xBinMinusQbeta = Qbeta - xBin;
            double scale = 0.0;
            
            if(xBinMinusQbeta <= gxpla.front())
            {
                scale = toyY.front();
            }
            
            else if(xBinMinusQbeta >= gxpla.back())
            {
                scale = toyY.back();
            }
            
            else
            {
                /* int k = 0;

                   for(int j=0;j<nPointsPla-1;++j)

                   {
                   if(xBinMinusQbeta >= gxpla[j] && xBinMinusQbeta < gxpla[j+1]) { k = j; break; }
                   }

                   double x0 = gxpla[k], x1 = gxpla[k+1];
                   double y0 = toyY[k], y1 = toyY[k+1];
                   double t = (xBinMinusQbeta - x0) / (x1 - x0);

                   scale = y0 + t*(y1 - y0);
                   */           

                scale = SplineToy->Eval(xBinMinusQbeta);
            }
            
            double inv = 1. / scale;
            double C = hist_E_corrected_mon->GetBinContent(ib);
            double val = C * inv;

            sumpla[ib]  += val;
            sum2pla[ib] += val*val;
            sumInvpla[ib]  += inv;
            sumInv2pla[ib] += inv*inv;
        }
    }

    for(int ib=1; ib<=nBins; ++ib)
    {
       double mean = sumpla[ib] / double(NtoysPla);
       double mean2 = sum2pla[ib] / double(NtoysPla);
       double var = mean2 - mean*mean;

       if(var < 0 && var > -1e-18) var = 0;
       
       double rms = (var>0) ? sqrt(var) : 0.0;

       double sigmaC = hist_E_corrected_mon->GetBinError(ib);

       double meanInv2 = sumInv2pla[ib] / double(NtoysPla);   
    
       double sigma_from_original = sigmaC * sqrt(meanInv2);
       double total_err = sqrt(rms*rms + sigma_from_original*sigma_from_original);
       
       hist_E_corrected_MonAndPla->SetBinContent(ib, mean);
       hist_E_corrected_MonAndPla->SetBinError(ib, total_err);
    }
  
    //hist_E_corrected_MonAndPla->RebinX(5);

    double ResultMean;
    double ResultMeanErr;

    double meansum = 0.;
    double sumcontent = 0.;
    double sumvar = 0.;
    
    for(int i = hist_E_corrected_MonAndPla->GetXaxis()->FindBin(0.); i <= hist_E_corrected_MonAndPla->GetXaxis()->FindBin(5.2); ++i)
    {
        double x = hist_E_corrected_MonAndPla->GetBinCenter(i);
        double y = hist_E_corrected_MonAndPla->GetBinContent(i);
        double err = hist_E_corrected_MonAndPla->GetBinError(i);

        meansum += x * y;
        sumcontent += y;
        sumvar += (err * err) * (x * x);
    }

    ResultMean = meansum / sumcontent;
    ResultMeanErr = sqrt(sumvar) / fabs(sumcontent);
   
    cout << "Mean = " << ResultMean << " +/- " << ResultMeanErr << endl;

    delete hSum2;
    delete hSum;
    
    fileout->Write();
    fileout->Close();
}
