#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TCutG.h"
#include "TRandom.h"

void unfolding()
{
  double energy;
  double ring;
  double energyinit;

  double testring;
  
  double sum;
  double arg;
  double addeff;

  double chi2 = 1.;

  vector<double> exp;
  vector<double> eff;

  TFile* data_file = TFile::Open("response.root");
  TFile* test_file = TFile::Open("252Cf.root");
    
  TTree* tree = (TTree*)data_file->Get("Rings");
  TTree* treetest = (TTree*)test_file->Get("Rings");
  TTree* treeinit = (TTree*)test_file->Get("PrimaryEnergy");
  
  tree->SetBranchAddress("NeutronEnergy", &energy);
  tree->SetBranchAddress("RingN", &ring);

  treetest->SetBranchAddress("RingN", &testring);
  treeinit->SetBranchAddress("Energy", &energyinit);
  
  double entries = tree->GetEntries();
  
  TH1D* p0 = new TH1D("p0","p0", 500,0.,5.01);
  
  TH2D* response = new TH2D("response","response", 4,1,5,500,0.,5.01);
  
  TH2D* inverse = new TH2D("inverse","inverse", 4,1,5,500,0.,5.01);
  
  TH1D* causes = new TH1D("causes", "causes", 500, 0., 5.01);
  TH1D* probas = new TH1D("probas", "probas", 500, 0., 5.01);

  TH1D* testrings = new TH1D("testrings","testrings", 4,1,5);

  TH1D* inithist = new TH1D("inithist", "inithist", 1000, 0., 10.01);

  //To get experimental ring numbers from root file

  for(int i=0; i<treetest->GetEntries(); i++)
  {
    treetest->GetEntry(i);
    testrings->Fill(testring);
  }

  for(int i=0; i<treeinit->GetEntries(); i++)
  { 
    treeinit->GetEntry(i);
    inithist->Fill(energyinit);
  }

  inithist->Scale(1./inithist->GetEntries());

  inithist->Draw("hist");

  for(int i=1; i<5; i++)
  {
    exp.push_back(testrings->GetBinContent(i));
  }

  //To define ring numbers by hand

  /*exp.push_back(301397.);
  exp.push_back(215293.);
  exp.push_back(103025.);
  exp.push_back(51063.);*/

  //Defining first p0 : Uniform distribution

  for(int i=1; i<501; i++)
  { 
    p0->SetBinContent(i, 1./500.);
    
    //p0->Draw();
  }

  //Constructing response matrix (TH2D)

  for(int i=0; i <= entries; i++)
  {
    tree->GetEntry(i);
    
    response->Fill(ring, energy);
  }

  response->Scale(1./10000.);

  //Constructing efficiency vector
  
  for(int i=1; i<501; i++)
  {
    addeff = 0;
    
    for(int j=1; j<5; j++)
    { 
      addeff = addeff + response->GetBinContent(j,i);
    }
    
    eff.push_back(addeff/(10000.));
  }

  //Starting while for Chi2 reduction

  while(chi2 > 1e-10)
  {
    double oug = 0.;
    
    for(int i=1; i<5; i++)
    {
      sum = 0.;

      for(int j=1; j<501; j++)
      { 
        sum = sum + (response->GetBinContent(i,j))*(p0->GetBinContent(j));
      }
      
      for(int k=1; k<501; k++)
      {
        double numerator = (response->GetBinContent(i,k))*(p0->GetBinContent(k));
        
        inverse->SetBinContent(i,k,(numerator/sum));
      }
    }
    
    //inverse->Draw("colz");
    
    for(int i=1; i<501; i++)
    {
      arg=0.;

      for(int j=1; j<5; j++)
      { 
        arg = arg + (exp[j-1]*inverse->GetBinContent(j,i));
      }
      
      causes->SetBinContent(i, (arg/eff[i-1]));
    }

    //causes->Draw("hist");

    for(int i=1; i<501; i++)
    {
      oug = oug + causes->GetBinContent(i);
    }

    for(int i=1; i<501; i++)
    {
      probas->SetBinContent(i,causes->GetBinContent(i)/oug); 
    }

    //probas->Draw("sameshist");

    chi2 = probas->Chi2Test(p0,"WWCHI2");

    cout << chi2 << endl;

    for(int i=1; i<501; i++)
    {
      p0->Clear();
      p0->SetBinContent(i,probas->GetBinContent(i));
    }

    inverse->Clear();
    causes->Clear();
    probas->Clear();
  }

  //causes->ResetStats();
  //causes->Draw("sameshist");

  p0->ResetStats();
  p0->Draw("sameshist");
  p0->SetLineColor(kRed);
}
