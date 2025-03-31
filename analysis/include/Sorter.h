#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>

using namespace std;

class Sorter{

	private:

	TFile *data_file;
	TFile *output_file;

	TTree *raw_data_tree;
	TTree *output_tree_single;
	TTree *output_tree_coinc;

	Int_t raw_energy1;
	Int_t raw_energy2;
	Int_t raw_energy3;
	ULong64_t raw_time;
	UShort_t raw_det_nbr;

	Long64_t fEntries;
	Long64_t fEntry;

	Long64_t lastevent;
	Long64_t eventafter;
	Long64_t eventbefore;

    Double_t fGe_E;
    UShort_t fGe_Id;
    Double_t fGe_Time;

    UShort_t fPlastic_Id;
    Double_t fPlastic_Time;

    UShort_t fMonster_Id;
    Double_t fMonster_Time;

    UShort_t fMrBig_Id;
    Double_t fMrBig_Time;

    //Gamma calibration coefficients
    Double_t a[4] = {0.018830012, 0.021904491, 0.034686097, 0.0119};
    Double_t b[4] = {0.155951777, -0.00662418, 1.089112811, -0.785};

    //Monster charge calibration coefficients
    std::string MonsterCalibFile = "monster_calib_params.txt";
    ifstream calib_monster_coeff;
    Double_t monster_a[39];
    Double_t monster_b[39];

    //Monster time alignement
    Double_t time_offset[39]={25.59,28.32,28.2,26.71,25.19,24.23,22.58,23.87,24.19,27.8,28.63,25.64,23.01,36.17,24.08,23.46,22.54,27.17,28.59,23.87,25.81,42.35,31.12,24.75,21.19,19.5,28.11,24.58,27.54,26.32,12.42,17.19,19.57,18.8,22.22,28.02,22.69,22.22,27.86};



        //Single

    std::vector<Double_t> fPlastic_tSingle;

    std::vector<Double_t> fGe_ESingle;
    
    std::vector<UShort_t> fGe_IdSingle;
        
        //Coincidences

    std::vector<Long64_t> double_events_check;

    Int_t fEventMult;

    std::vector<Double_t> fGePlastic_ECond;
    std::vector<Double_t> fGePlastic_tCond;
    std::vector<Double_t> fGePlastic_tDiff;
    std::vector<Double_t> fGePlastic_Id;
    Int_t fGePlastic_Mult;

    std::vector<Double_t> fMonsterPlastic_tDiff;
    std::vector<Double_t> fMonsterPlastic_tCond;
    std::vector<Double_t> fMonsterPlastic_Q1Cond;
    std::vector<Double_t> fMonsterPlastic_Q2Cond;
    std::vector<Double_t> fMonsterPlastic_Q3Cond;
    std::vector<Double_t> fMonsterPlastic_Id;
    Int_t fMonsterPlastic_Mult;

    std::vector<Double_t> fMonsterPlastic_Q2Cond_calibrated;

    std::vector<Double_t> fMrBigPlastic_ECond;
    std::vector<Double_t> fMrBigPlastic_tCond;
    std::vector<Double_t> fMrBigPlastic_tDiff;
    Int_t fMrBigPlastic_Mult;

	public:

	Sorter(const char*, const char*);
	virtual ~Sorter();

	void Load_rawdata(const char*);

	protected:

	virtual void LoadMonsterCalibrationCoefficients();
	virtual void SetTreesAndBranches(const char*);
	virtual void FillSingleBranches();
	virtual void FillCoincBranches(Double_t, Double_t);
	virtual Double_t Ge_alignement(Int_t, Int_t);
	virtual Double_t MrBig_alignement(Int_t);
	virtual Double_t Monster_alignement(Int_t, Int_t);
	virtual Double_t Time_alignement(Int_t);
	virtual void ResetVar();
	virtual void ClearVectors();
	virtual void SetVar(UShort_t);
};
