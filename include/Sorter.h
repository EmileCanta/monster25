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

    Double_t fMrBig_E;

    UShort_t fPlastic_Id;
    Double_t fPlastic_Time;

    UShort_t fMonster_Id;
    Double_t fMonster_Time;

    UShort_t fMrBig_Id;
    Double_t fMrBig_Time;

    //Gamma calibration coefficients
    Double_t a[4] = {0.018830012, 0.021904491, 0.034686097, 1.};
    Double_t b[4] = {0.155951777, -0.00662418, 1.089112811, 0.};

    //Monster charge calibration coefficients
    std::string MonsterCalibFile = "monster_calib_params.txt";
    ifstream calib_monster_coeff;
    Double_t monster_a[39];
    Double_t monster_b[39];

    //Monster time alignement
    Double_t time_offset_plastic[39]={-2.68552, -5.241, -5.06835, -3.73247, -1.96254, -1.28285, 0.622969, -0.665302, -0.626229, -3.03628, -5.67668, -2.3431, 0.104158, -12.9122, -1.31232, -0.788829, 0.932808, -3.86213, -5.08842, -0.407243, -2.24909, -19.0685, -7.79083, -1.36452, 2.31159, 3.81558, -5.06251, -1.41587, -4.54087, -3.01263, 10.7251, 5.92087, 3.38479, 4.29565, 1.28501, -4.43648, 0.442536, 1.18371, -4.16247};
    Double_t time_offset_mrbig[39]={-27.3, -29.5, -29.4, -28.1, -26.7, -25.8, -23.7, -25.2, -25.3, -27.6, -30.4, -26.6, -24.3, -37.1, -25.8, -25, -23.4, -28.2, -29.9, -25.1, -26.5, -43.5, -32.1, -26, -22.1, -20.7, -29.5, -25.7, -29.1, -27.6, -14, -18.7, -21, -20.3, -23, -29, -24.1, -23.3, -29};
        
        //Single

    std::vector<Double_t> fPlastic_tSingle;

    std::vector<Double_t> fGe_ESingle;
    std::vector<Double_t> fMrBig_ESingle;
    
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
    virtual void ForwardCoinc(Double_t);
    virtual void BackwardCoinc(Double_t);
	virtual Double_t Ge_alignement(Int_t, Int_t);
	virtual Double_t MrBig_alignement(Int_t);
	virtual Double_t Monster_alignement(Int_t, Int_t);
	virtual Double_t Time_alignement(Int_t);
	virtual void ResetVar();
	virtual void ClearVectors();
	virtual void SetVar(UShort_t);
};
