#include "TFile.h"
#include "TF1.h"
#include "fstream"

#include <math.h>

using namespace std;

void PlotData()
{
	TCanvas* c1 = new TCanvas();

    TFile *fileeff = TFile::Open("/Users/cantacuzene/phd/analysis/monster25/root_files/effmonster.root", "READ");

    TGraph* grapheff = (TGraph*)fileeff->Get("Graph");

    fstream file;
    double x;
    double y;
    double z;
    TGraph* graph = new TGraph();

    file.open("../divers/mannhart.dat");

    TGraph* g1 = new TGraph();

    while(1)
    {
        file >> x >> y >> z;

        if(file.eof()) break;

        g1->AddPoint(x,y*grapheff->Eval(x));
    }

g1->Draw();
g1->SetName("mannhart");

}
