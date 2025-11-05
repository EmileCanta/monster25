void TransformXAxis() {
    // Ouvrir un fichier ROOT (si nécessaire)
    // TFile *file = TFile::Open("mon_fichier.root");

    // Histogramme original (par exemple généré ou lu d'un fichier)
    TFile *filein = TFile::Open("~/figures_these/C/tof_252cf_nosub.root", "READ");

    TFile *fileeff = TFile::Open("/Users/cantacuzene/monster25/divers/effmonster.root", "READ");

    TH1D* h_original = (TH1D*)filein->Get("hbis__1");

    TGraph* graph = (TGraph*)fileeff->Get("Graph");

    int nBins = h_original->GetNbinsX();

    std::vector<double> new_bins(nBins + 1);

    for (int i = 22000; i <= 23000; ++i) {
        double x_edge = h_original->GetBinLowEdge(i + 1);
        if (x_edge != 0) new_bins[i] = 1.21393103e4 / (x_edge * x_edge);
    }

    std::sort(new_bins.begin(), new_bins.end());

    // Créer un nouvel histogramme avec les nouveaux bins
    TH1D* h_transformed = new TH1D("h_transformed", "Histogram with x -> 1/x^{2};1/x^{2};Entries", nBins, new_bins.data());

    // Remplir le nouvel histogramme
    for (int i = 22000; i <= 23000; ++i) {
        double content = h_original->GetBinContent(i);
        double x = h_original->GetBinCenter(i);
        if (x != 0) 
        {
            double new_x = 1.21393103e4 / (x * x);
            //if(x >= 10. && x <=500.)
            {
                h_transformed->Fill(new_x, content);
            }
        }
    }

    // Dessiner
    TCanvas* c = new TCanvas("c", "Transformed Axis", 800, 600);
    h_transformed->Draw("hist");
}