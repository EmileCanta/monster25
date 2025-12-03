void MultiplyHistByGraph_Toy() 
{
    TFile *fileeff = TFile::Open("~/phd/analysis/monster25/root_files/effmonster.root", "READ");
    TFile *fileeffplastic = TFile::Open("~/phd/analysis/monster25/root_files/effplastic.root", "READ");
    TFile *filehist = TFile::Open("~/phd/analysis/monster25/root_files/E_Spectrum_82Ga_0312.root", "READ");
    TFile *filemanhart = TFile::Open("~/phd/analysis/monster25/root_files/manhart.root", "READ");
    TFile *filebgd_manual = TFile::Open("~/phd/analysis/monster25/root_files/BgdModel.root", "READ");
    TFile *fileibu = TFile::Open("/home/emile/phd/analysis/mcnp_tetra/root_files/tests_regina/gaibu.root", "READ");

    TGraph* grapheff = (TGraph*)fileeff->Get("Graph");
    TGraph* grapheffpla = (TGraph*)fileeffplastic->Get("Graph");
    TGraph* graphmanhart = (TGraph*)filemanhart->Get("g1");
    TH1D* hist_all = (TH1D*)filehist->Get("hist_all");
    TH1D* hist_bgd = (TH1D*)filehist->Get("hist_bgd");
    TH1D* hist_bgd_manual = (TH1D*)filebgd_manual->Get("hist_model_bgd");
    TH1D* hist_ibu = (TH1D*)fileibu->Get("proba");

    TGraph* graphtoys = new TGraph();
    graphtoys->SetTitle("graphtoys");

    //hist_all->Add(hist_bgd_manual, -1);
    //hist_all->Rebin(10);

    TRandom3* rnd = new TRandom3(0);

    Int_t Ntoy = 2000;

    TH1D* hOut = (TH1D*)hist_all->Clone(Form("%s_x%s_MC", hist_all->GetName(), grapheff->GetName()));
    hOut->Reset(); 
    //hOut->Sumw2();

    const int nBins = hist_all->GetNbinsX();
    const int nPoints = grapheff->GetN();

    std::vector<double> sum(nBins+1, 0.0);
    std::vector<double> sum2(nBins+1, 0.0);
    
    std::vector<double> sumInv(nBins+1, 0.0);
    std::vector<double> sumInv2(nBins+1, 0.0);

    std::vector<double> gx(nPoints), gy(nPoints), gey(nPoints);
    for (int i=0;i<nPoints;++i) {
        grapheff->GetPoint(i, gx[i], gy[i]);
        
        gey[i] = grapheff->GetErrorY(i);
        //cout << gy[i] << endl;
    }

    for (int itoy=0; itoy<Ntoy; ++itoy) {
        std::vector<double> toyY(nPoints);
        for (int ip=0; ip<nPoints; ++ip) {
            double sigma = gey[ip];
            toyY[ip] = (sigma > 0.0) ? gy[ip] + rnd->Uniform(-sigma, sigma) : gy[ip];
            graphtoys->AddPoint(gx[ip],toyY[ip]);
        }

        for (int ib=1; ib<=nBins; ++ib) {
            double xBin = hist_all->GetXaxis()->GetBinCenter(ib);
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
            double C = hist_all->GetBinContent(ib);
            double val = C * inv;

            sum[ib]  += val;
            sum2[ib] += val*val;
            sumInv[ib]  += inv;
            sumInv2[ib] += inv*inv;
        }
    }

    for (int ib=1; ib<=nBins; ++ib) {
       double mean = sum[ib] / double(Ntoy);
       double mean2 = sum2[ib] / double(Ntoy);
       double var = mean2 - mean*mean;
       if (var < 0 && var > -1e-18) var = 0;
       double rms = (var>0) ? sqrt(var) : 0.0;

       double sigmaC = hist_all->GetBinError(ib);

       double meanInv2 = sumInv2[ib] / double(Ntoy);   
    
       double sigma_from_original = sigmaC * sqrt(meanInv2);
       double total_err = sqrt(rms*rms + sigma_from_original*sigma_from_original);
          
       //cout << sigmaC << " " << sigma_from_original << " " << rms << " " << total_err << endl;
       //cout << total_err/mean << endl;
       
       hOut->SetBinContent(ib, mean);
       hOut->SetBinError(ib, total_err);
    }
   
    //hOut->Scale(1./hOut->Integral());
    
    hOut->Draw("E");
    //TCanvas* c2 = new TCanvas();
    //graphtoys->Draw("AL");
    //hist_all->Draw("sameE");
    //graphmanhart->Draw("same");
}
