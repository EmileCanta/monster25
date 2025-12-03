Double_t f1(Double_t *x, Double_t *par)
{
    return 1.98648e2;;
}

void Test()
{	
    TFile *fileeffmon = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *fileeffpla = TFile::Open("~/phd/analysis/monster25/root_files/effplastic.root", "READ");
    //TFile *filein = TFile::Open("~/phd/data/monster25/sorted/82Ga/AllRuns.root", "READ");  //For 82Ga
    TFile *filein = TFile::Open("~/phd/data/monster25/sorted/cfruns/RUN111.root", "READ"); //For 252Cf
    
    TGraph* grapheffmon = (TGraph*)fileeffmon->Get("Graph");
    TGraph* grapheffpla = (TGraph*)fileeffpla->Get("Graph");

    TFile *cutfile = TFile::Open("~/phd/analysis/monster25/root_files/cutfile.root","READ");
    
    TCutG *cut = new TCutG();
    
    //TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/82Ga_ana.root","RECREATE");
    TFile *fileout = new TFile("~/phd/analysis/monster25/root_files/252Cf_ana.root","RECREATE");
    
    TTree* tree = (TTree*)filein->Get("tcoinc");
    
    TF1* function = new TF1("f1", f1, 0, 20);

    TH1D *hist_model_bgd = new TH1D("hist_model_bgd", "hist_model_bgd", 625, -1000, 1000);
    hist_model_bgd->Sumw2();

    for(int i=1; i<=625; i++)
	{
		hist_model_bgd->SetBinContent(i,function->Eval((float)i*0.3125));
		hist_model_bgd->SetBinError(i,1.33777);
	}

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

    TH1D *hist_tof_all = new TH1D("hist_tof_all", "hist_tof_all", 625, -1000, 1000);
    TH1D *hist_tof_bgd = new TH1D("hist_tof_bgd", "hist_tof_bgd", 625, -1000, 1000);

    int fEntries = tree->GetEntries();

    //double windowback = 50; //For 82Ga
    double windowback = 25; //For 252Cf
    double windowfront = 200;

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
                        hist_tof_bgd->Fill(tdiff->at(k));
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
    
    //hist_tof_sub->Add(hist_model_bgd,-1); //For 82Ga
    hist_tof_sub->Add(hist_tof_bgd,-1); //For 252Cf

    TH1D* hist_E = new TH1D("hist_E", "hist_E", 500, 0, 20);

    hist_E->Sumw2();
    
    double massn = 1.67492750056e-27;
    double joultoMeV = 6.241509343260e12;
    double d = 1.575;

    int NtoysPerBin = 2000;

    for(int ib=1; ib<=hist_tof_all->GetNbinsX(); ++ib)
    {
        double Ni = hist_tof_all->GetBinContent(ib);
        double errNi = hist_tof_all->GetBinError(ib);
        
        if (Ni <= 0) continue;

        double t_center = hist_tof_all->GetXaxis()->GetBinCenter(ib);

        double weight = Ni / double(NtoysPerBin);

        for(int k=0; k<NtoysPerBin; ++k)
        {
            double t_sample = gRandom->Gaus(t_center*1e-9, 2.1e-9);
            
            double K = joultoMeV * 0.5 * massn * d * d / (t_sample * t_sample);
            
            hist_E->Fill(K, weight);
        }
    }

    TGraph* graphtoysmon = new TGraph();
    graphtoysmon->SetTitle("graphtoysmon");
    
    TRandom3* rnd = new TRandom3(0);

    Int_t NtoyMon = 2000;

    TH1D* hist_E_corrected_mon = (TH1D*)hist_E->Clone("hist_E_corrected_mon");
    hist_E_corrected_mon->Reset(); 

    const int nBins = hist_E->GetNbinsX();
    const int nPoints = grapheffmon->GetN();

    std::vector<double> sum(nBins+1, 0.0);
    std::vector<double> sum2(nBins+1, 0.0);
    
    std::vector<double> sumInv(nBins+1, 0.0);
    std::vector<double> sumInv2(nBins+1, 0.0);

    std::vector<double> gx(nPoints), gy(nPoints), gey(nPoints);
    for (int i=0;i<nPoints;++i) {
        grapheffmon->GetPoint(i, gx[i], gy[i]);
        
        gey[i] = grapheffmon->GetErrorY(i);
        //cout << gy[i] << endl;
    }

    for (int itoy=0; itoy<NtoyMon; ++itoy) {
        std::vector<double> toyY(nPoints);
        for (int ip=0; ip<nPoints; ++ip) {
            double sigma = gey[ip];
            toyY[ip] = (sigma > 0.0) ? gy[ip] + rnd->Uniform(-sigma, sigma) : gy[ip];
            graphtoysmon->AddPoint(gx[ip],toyY[ip]);
        }

        for (int ib=1; ib<=nBins; ++ib) {
            double xBin = hist_E->GetXaxis()->GetBinCenter(ib);
            //cout << xBin << endl;
            double scale = 0.0;
            if (xBin <= gx.front()) {
                scale = toyY.front();

                //cout << "IIII" << endl;
            } else if (xBin >= gx.back()) {
                scale = toyY.back();
            } else {
                int k = 0;
                for (int j=0;j<nPoints-1;++j) {
                    if (xBin >= gx[j] && xBin < gx[j+1]) { k = j; break; }
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

    for (int ib=1; ib<=nBins; ++ib) {
       double mean = sum[ib] / double(NtoyMon);
       double mean2 = sum2[ib] / double(NtoyMon);
       double var = mean2 - mean*mean;
       if (var < 0 && var > -1e-18) var = 0;
       double rms = (var>0) ? sqrt(var) : 0.0;

       double sigmaC = hist_E->GetBinError(ib);

       double meanInv2 = sumInv2[ib] / double(NtoyMon);   
    
       double sigma_from_original = sigmaC * sqrt(meanInv2);
       double total_err = sqrt(rms*rms + sigma_from_original*sigma_from_original);
          
       //cout << sigmaC << " " << sigma_from_original << " " << rms << " " << total_err << endl;
       //cout << total_err/mean << endl;
       
       hist_E_corrected_mon->SetBinContent(ib, mean);
       hist_E_corrected_mon->SetBinError(ib, total_err);
    }
    
    fileout->Write();
}
