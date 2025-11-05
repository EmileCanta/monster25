{
    TChain* ch = new TChain("DataTree");
    for(unsigned int i=40;i<41;i++)
        ch->Add(Form("/data_tampon/MONSTER25/MONSTER25/root_files/run%i_raw_kept.root", i));
    ch->Draw("multiplicity:time@.at(0)","","colz");
}

