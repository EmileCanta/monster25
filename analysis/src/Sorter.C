#include "../include/Sorter.h"

Sorter::Sorter(const char* InputFileName, const char* OutputFileName)
{
    cout << "Starting Sorter" << endl;

    gROOT->ProcessLine("#include <vector>");

    LoadMonsterCalibrationCoefficients();

    Load_rawdata(InputFileName);

    SetTreesAndBranches(OutputFileName);

    FillSingleBranches();

    output_tree_single->Write();

    FillCoincBranches(500e3,1500e3);

    output_tree_coinc->Write();

    data_file->Close();
}

//******************************************************************************
//******************************************************************************

Sorter::~Sorter()
{
    ;
}

//******************************************************************************
//******************************************************************************

Double_t Sorter::Ge_alignement(Int_t channel, Int_t index)
{
    Double_t channel_d = (Double_t)channel;
    Double_t channel_aligned = 0.;
    Double_t ch = (channel_d + gRandom->Uniform(1.0) - 0.5);
    channel_aligned = ch * a[index-43] + b[index-43];
    return channel_aligned;

}//******************************************************************************
//******************************************************************************

Double_t Sorter::Time_alignement(Int_t index)
{
    Double_t time_aligned = time_offset[index-3] - 28.2;

    return time_aligned;
}



//******************************************************************************
//******************************************************************************

void Sorter::LoadMonsterCalibrationCoefficients()
{
    calib_monster_coeff.open(MonsterCalibFile);
    calib_monster_coeff.ignore(1000, '\n');
    calib_monster_coeff.ignore(1000, '\n');
    for (unsigned int k=0;k<39;k++)
    {
        calib_monster_coeff >> monster_a[k] >> monster_b[k];
    }
    calib_monster_coeff.close();
    return;
}

//******************************************************************************
//******************************************************************************

Double_t Sorter::Monster_alignement(Int_t channel, Int_t index)
{
    Double_t channel_d = (Double_t)channel;
    Double_t channel_aligned = 0.;
    Double_t ch = (channel_d + gRandom->Uniform(1.0) - 0.5);
    channel_aligned = ch * monster_a[index-3] + monster_b[index-3];
    return channel_aligned;

}
//******************************************************************************
//******************************************************************************

Double_t Sorter::MrBig_alignement(Int_t channel)
{
    Double_t channel_d = (Double_t)channel;
    Double_t channel_aligned = 0.;
    Double_t ch = (channel_d + gRandom->Uniform(1.0) - 0.5);
    channel_aligned = ch * a[3] + b[3];
    return channel_aligned;
}

//******************************************************************************
//******************************************************************************

void Sorter::ResetVar()
{
    fPlastic_Id = -666;
    fPlastic_Time = -666;
    fMonster_Id = -666;
    fMonster_Time = -666;
    fGe_Id = -666;
    fGe_Time = -666;
}

//******************************************************************************
//******************************************************************************

void Sorter::ClearVectors()
{
        //Singles

    //Time
    fPlastic_tSingle.clear();

    //Energy
    fGe_ESingle.clear();

    //Ge Id
    fGe_IdSingle.clear();

        //Coincidences

    //Ge-Plastic
    fGePlastic_tDiff.clear();
    fGePlastic_ECond.clear();
    fGePlastic_tCond.clear();
    fGePlastic_Id.clear();

    //Monster-Plastic
    fMonsterPlastic_tDiff.clear();
    fMonsterPlastic_tCond.clear();
    fMonsterPlastic_Q1Cond.clear();
    fMonsterPlastic_Q2Cond.clear();
    fMonsterPlastic_Q3Cond.clear();
    fMonsterPlastic_Id.clear();

    //Monster-MrBig
    fMrBigPlastic_tDiff.clear();
    fMrBigPlastic_ECond.clear();
    fMrBigPlastic_tCond.clear();

}

//******************************************************************************
//******************************************************************************

void Sorter::SetVar(UShort_t det)
{
    if(det == 1)
    {
        fPlastic_Id = det;
        fPlastic_Time = (Double_t)raw_time;
    }

    if(det == 2)
    {
        fMrBig_Id = det;
        fMrBig_Time = (Double_t)raw_time;
    }

    if(det > 2 && det < 42)
    {
        fMonster_Id = det;
        fMonster_Time = (Double_t)raw_time;
    }

    if(det >= 43 && det <= 45)
    {
        fGe_Id = det;
        fGe_Time = (Double_t)raw_time;
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::Load_rawdata (const char* InputFileName)
{
    cout << "Loading Narval Tree" << endl;

    data_file = TFile::Open(InputFileName);

    raw_data_tree = (TTree*)data_file->Get("DataTree");

    raw_data_tree->SetBranchAddress("nrj", &raw_energy1);
    raw_data_tree->SetBranchAddress("nrj2", &raw_energy2);
    raw_data_tree->SetBranchAddress("nrj3", &raw_energy3);
    raw_data_tree->SetBranchAddress("time", &raw_time);
    raw_data_tree->SetBranchAddress("label", &raw_det_nbr);

    fEntries = raw_data_tree->GetEntries();

    std::cout << "Entries to be sorted : " << fEntries << std::endl;
}

//******************************************************************************
//******************************************************************************

void Sorter::SetTreesAndBranches(const char* OutputFileName)
{
    std::cout << "Initialization of output TTree" << std::endl;

    std::cout << "Creating branches" << std::endl;

    output_file = new TFile(OutputFileName, "recreate");

    output_tree_single = new TTree("tsingle", "Tree with singles");
    output_tree_coinc = new TTree("tcoinc", "Tree with coincidences");

        //Single

    //Time
    output_tree_single->Branch("Plastic_tSingle", &fPlastic_tSingle);

    //Energy
    output_tree_single->Branch("Ge_ESingle", &fGe_ESingle);

    //Ge Id
     output_tree_single->Branch("Ge_IdSingle", &fGe_IdSingle);

        //Coincidences

    //Ge-Plastic
    output_tree_coinc->Branch("GePlastic_ECond", &fGePlastic_ECond);
    output_tree_coinc->Branch("GePlastic_tCond", &fGePlastic_tCond);
    output_tree_coinc->Branch("GePlastic_tDiff", &fGePlastic_tDiff);
    output_tree_coinc->Branch("GePlastic_Id", &fGePlastic_Id);
    output_tree_coinc->Branch("GePlastic_Mult", &fGePlastic_Mult);

    //Monster-Plastic
    output_tree_coinc->Branch("MonsterPlastic_Q1Cond", &fMonsterPlastic_Q1Cond);
    output_tree_coinc->Branch("MonsterPlastic_Q2Cond", &fMonsterPlastic_Q2Cond);
    output_tree_coinc->Branch("MonsterPlastic_Q3Cond", &fMonsterPlastic_Q3Cond);
    output_tree_coinc->Branch("MonsterPlastic_tCond", &fMonsterPlastic_tCond);
    output_tree_coinc->Branch("MonsterPlastic_tDiff", &fMonsterPlastic_tDiff);
    output_tree_coinc->Branch("MonsterPlastic_Id", &fMonsterPlastic_Id);
    output_tree_coinc->Branch("MonsterPlastic_Mult", &fMonsterPlastic_Mult);

    //Monster-MrBig
    output_tree_coinc->Branch("MrBigPlastic_ECond", &fMrBigPlastic_ECond);
    output_tree_coinc->Branch("MrBigPlastic_tCond", &fMrBigPlastic_tCond);
    output_tree_coinc->Branch("MrBigPlastic_tDiff", &fMrBigPlastic_tDiff);
    output_tree_coinc->Branch("MrBigPlastic_Mult", &fMrBigPlastic_Mult);

    output_tree_coinc->Branch("EventMult", &fEventMult);
}

//******************************************************************************
//******************************************************************************

void Sorter::FillSingleBranches()
{
    std::cout << "Sorting single data" << std::endl;

    ResetVar();

    for(fEntry = 0; fEntry < fEntries; fEntry++)
    {
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();
        
        if(raw_det_nbr == 1)
        {
            fPlastic_tSingle.push_back((Double_t)raw_time / 1.e9);
        }

        if(raw_det_nbr >= 43 && raw_det_nbr <= 45)
        {
            fGe_ESingle.push_back(Ge_alignement(raw_energy1, raw_det_nbr));
            fGe_IdSingle.push_back(raw_det_nbr);
        }

        output_tree_single->Fill();

        if(fEntry%100==0)
        {
            std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
        }
    }
}

//******************************************************************************
//******************************************************************************

void Sorter::FillCoincBranches(Double_t backward_window, Double_t forward_window)
{
    std::cout << "Sorting Plastic coinc data" << std::endl;

    lastevent = 0;

    ResetVar();

    for(fEntry = 0; fEntry < fEntries; fEntry++)
    {
        raw_data_tree->GetEntry(fEntry);

        ClearVectors();

        SetVar(raw_det_nbr);

        eventafter = 0;
        eventbefore = 0;

        fMrBigPlastic_Mult = 1;
        fMonsterPlastic_Mult = 1;
        fGePlastic_Mult = 1;
        fEventMult = 1;

        if(fEntry < lastevent) continue;

        while(TMath::Abs((Double_t)raw_time -fPlastic_Time) < backward_window) //backward
        {
            raw_data_tree->GetEntry(fEntry - eventbefore++);

            if(fEntry - eventbefore < 0) break;

            if(TMath::Abs((Double_t)raw_time - fPlastic_Time) > backward_window) continue;
            if(TMath::Abs((Double_t)raw_time - fPlastic_Time) < backward_window)
            {
                if(std::find(double_events_check.begin(), double_events_check.end(), (fEntry - eventbefore)) != double_events_check.end()) continue;

                if(raw_det_nbr == 2)
                {
                    fMrBigPlastic_tDiff.push_back(((Double_t)raw_time - fPlastic_Time) / 1.e3); //nano
                    fMrBigPlastic_ECond.push_back(MrBig_alignement(raw_energy1));
                    fMrBigPlastic_tCond.push_back((Double_t)raw_time / 1.e9); //mili;
                    fMrBigPlastic_Mult++;
                    fEventMult++;
                }

                if(raw_det_nbr >= 3 && raw_det_nbr <= 41)
                {
                    fMonsterPlastic_tDiff.push_back((((Double_t)raw_time - fPlastic_Time / 1e3)) - Time_alignement(raw_det_nbr)); //nano
                    fMonsterPlastic_Q1Cond.push_back(raw_energy1);
                    fMonsterPlastic_Q2Cond.push_back(raw_energy2);
                    fMonsterPlastic_Q3Cond.push_back(raw_energy3);
                    fMonsterPlastic_tCond.push_back((Double_t)raw_time / 1.e9); //mili;
                    fMonsterPlastic_Id.push_back(raw_det_nbr);
                    fMonsterPlastic_Mult++;
                    fEventMult++;
                }

                if(raw_det_nbr >= 43 && raw_det_nbr <= 45)
                {
                    fGePlastic_tDiff.push_back(((Double_t)raw_time - fPlastic_Time) / 1.e3); //nano
                    fGePlastic_ECond.push_back(Ge_alignement(raw_energy1, raw_det_nbr));
                    fGePlastic_tCond.push_back((Double_t)raw_time / 1.e9); //mili;
                    fGePlastic_Id.push_back(raw_det_nbr);
                    fGePlastic_Mult++;
                    fEventMult++;
                }
            }    
        }
        
        double_events_check.clear();

        raw_data_tree->GetEntry(fEntry);

        while(TMath::Abs((Double_t)raw_time -fPlastic_Time) < forward_window) //forward
        {
            raw_data_tree->GetEntry(fEntry + eventafter++);

            if(fEntry + eventafter > fEntries) break;

            if(TMath::Abs((Double_t)raw_time - fPlastic_Time) > forward_window) continue;
            if(TMath::Abs((Double_t)raw_time - fPlastic_Time) < forward_window)
            {
                if(raw_det_nbr == 2)
                {
                    fMrBigPlastic_tDiff.push_back(((Double_t)raw_time - fPlastic_Time) / 1.e3); //nano
                    fMrBigPlastic_ECond.push_back(MrBig_alignement(raw_energy1));
                    fMrBigPlastic_tCond.push_back((Double_t)raw_time / 1.e9); //mili;
                    fMrBigPlastic_Mult++;
                    fEventMult++;
                }

                if(raw_det_nbr >= 3 && raw_det_nbr <= 41)
                {
                    fMonsterPlastic_tDiff.push_back((((Double_t)raw_time - fPlastic_Time) / 1.e3) - Time_alignement(raw_det_nbr)); //nano
                    fMonsterPlastic_Q1Cond.push_back(raw_energy1);
                    fMonsterPlastic_Q2Cond.push_back(raw_energy2);
                    fMonsterPlastic_Q3Cond.push_back(raw_energy3);
                    fMonsterPlastic_tCond.push_back((Double_t)raw_time / 1.e9); //mili;
                    fMonsterPlastic_Id.push_back(raw_det_nbr);
                    fMonsterPlastic_Mult++;
                    fEventMult++;
                }

                if(raw_det_nbr >= 43 && raw_det_nbr <= 45)
                {
                    fGePlastic_tDiff.push_back(((Double_t)raw_time - fPlastic_Time) / 1.e3); //nano
                    fGePlastic_ECond.push_back(Ge_alignement(raw_energy1, raw_det_nbr));
                    fGePlastic_tCond.push_back((Double_t)raw_time / 1.e9); //mili;
                    fGePlastic_Id.push_back(raw_det_nbr);
                    fGePlastic_Mult++;
                    fEventMult++;
                }

                double_events_check.push_back(fEntry + eventafter);
            }   
        }

        lastevent = fEntry + eventafter;

        output_tree_coinc->Fill();

        std::cout << std::setprecision(3) << std::setw(5) << (100.*fEntry/fEntries) << " %\r";
    }
}