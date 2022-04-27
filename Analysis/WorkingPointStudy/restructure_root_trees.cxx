#include <vector>
#include <string>

void print_status(const unsigned int index, const unsigned int entries) {

  if(index % 1000000 == 0) cout << "Processing entry " << index << " / " << entries << endl;
}

void restructure_root_trees_qcd(const string & year) {

  cout << "Working on QCD" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_qcd = "QCD_HT200toInf_"+year+".root";

  const string infile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd;
  cout << "Open " << infile_path_qcd << endl;
  TFile *infile_qcd = TFile::Open(infile_path_qcd.c_str(), "READ");
  TTree *infile_tree_qcd = (TTree*)infile_qcd->Get("AnalysisTree");

  float event_weight;
  vector<float> *jets_pt;
  vector<float> *jets_msd;
  vector<float> *jets_subdeepcsv;
  vector<float> *jets_subdeepjet;
  vector<float> *jets_tau32;
  vector<float> *jets_tau21;
  vector<float> *jets_deepak8_TvsQCD;
  vector<float> *jets_deepak8_WvsQCD;
  vector<float> *jets_partnet_TvsQCD;
  vector<float> *jets_partnet_WvsQCD;

  infile_tree_qcd->SetBranchAddress("event_weight", &event_weight);
  infile_tree_qcd->SetBranchAddress("jets_pt", &jets_pt);
  infile_tree_qcd->SetBranchAddress("jets_msd", &jets_msd);
  infile_tree_qcd->SetBranchAddress("jets_subjets_deepcsv_max", &jets_subdeepcsv);
  infile_tree_qcd->SetBranchAddress("jets_subjets_deepjet_max", &jets_subdeepjet);
  infile_tree_qcd->SetBranchAddress("jets_tau32", &jets_tau32);
  infile_tree_qcd->SetBranchAddress("jets_tau21", &jets_tau21);
  infile_tree_qcd->SetBranchAddress("jets_deepak8_TvsQCD", &jets_deepak8_TvsQCD);
  infile_tree_qcd->SetBranchAddress("jets_deepak8_WvsQCD", &jets_deepak8_WvsQCD);
  infile_tree_qcd->SetBranchAddress("jets_partnet_TvsQCD", &jets_partnet_TvsQCD);
  infile_tree_qcd->SetBranchAddress("jets_partnet_WvsQCD", &jets_partnet_WvsQCD);

  const string outfile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd+".restructured";
  TFile *outfile_qcd = new TFile(outfile_path_qcd.c_str(), "RECREATE");
  TTree *outfile_tree_qcd = new TTree("my_tree", "new flat tree incorporating all QCD jets");

  float jet_weight;
  float jet_pt;
  float jet_msd;
  float jet_subdeepcsv;
  float jet_subdeepjet;
  float jet_tau32;
  float jet_tau21;
  float jet_deepak8_TvsQCD;
  float jet_deepak8_WvsQCD;
  float jet_partnet_TvsQCD;
  float jet_partnet_WvsQCD;

  outfile_tree_qcd->Branch("weight", &jet_weight, "weight/F");
  outfile_tree_qcd->Branch("pt", &jet_pt, "pt/F");
  outfile_tree_qcd->Branch("msd", &jet_msd, "msd/F");
  outfile_tree_qcd->Branch("subdeepcsv", &jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_qcd->Branch("subdeepjet", &jet_subdeepjet, "subdeepjet/F");
  outfile_tree_qcd->Branch("tau32", &jet_tau32, "tau32/F");
  outfile_tree_qcd->Branch("tau21", &jet_tau21, "tau21/F");
  outfile_tree_qcd->Branch("deepak8_TvsQCD", &jet_deepak8_TvsQCD, "deepak8_TvsQCD/F");
  outfile_tree_qcd->Branch("deepak8_WvsQCD", &jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_qcd->Branch("partnet_TvsQCD", &jet_partnet_TvsQCD, "partnet_TvsQCD/F");
  outfile_tree_qcd->Branch("partnet_WvsQCD", &jet_partnet_WvsQCD, "partnet_WvsQCD/F");

  for(int i = 0; i < infile_tree_qcd->GetEntries(); i++) {
    print_status(i, infile_tree_qcd->GetEntries());
    infile_tree_qcd->GetEntry(i);
    for(int j = 0; j < jets_pt->size(); j++) { // loop over jets in vector
      jet_weight = event_weight;
      jet_pt = jets_pt->at(j);
      jet_msd = jets_msd->at(j);
      jet_subdeepcsv = jets_subdeepcsv->at(j);
      jet_subdeepjet = jets_subdeepjet->at(j);
      jet_tau32 = jets_tau32->at(j);
      jet_tau21 = jets_tau21->at(j);
      jet_deepak8_TvsQCD = jets_deepak8_TvsQCD->at(j);
      jet_deepak8_WvsQCD = jets_deepak8_WvsQCD->at(j);
      jet_partnet_TvsQCD = jets_partnet_TvsQCD->at(j);
      jet_partnet_WvsQCD = jets_partnet_WvsQCD->at(j);
      outfile_tree_qcd->Fill();
    }
  }

  outfile_tree_qcd->Print();
  outfile_qcd->Write();
}


void restructure_root_trees_wjets(const string & year) {

  cout << "Working on wjets" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_wjets = "WJetsToQQ_HT200toInf"+year+".root";

  const string infile_path_wjets = sframe_output_path+file_prefix+file_postfix_wjets;
  cout << "Open " << infile_path_wjets << endl;
  TFile *infile_wjets = TFile::Open(infile_path_wjets.c_str(), "READ");
  TTree *infile_tree_wjets = (TTree*)infile_wjets->Get("AnalysisTree");

  float event_weight;
  float wnearestjet_pt;
  float wnearestjet_msd;
  float wnearestjet_tau21;
  float wnearestjet_deepak8_WvsQCD;
  float wnearestjet_partnet_WvsQCD;
  float wnearestjet_dr;

  infile_tree_wjets->SetBranchAddress("event_weight", &event_weight);
  infile_tree_wjets->SetBranchAddress("wnearestjet_pt", &wnearestjet_pt);
  infile_tree_wjets->SetBranchAddress("wnearestjet_msd", &wnearestjet_msd);
  infile_tree_wjets->SetBranchAddress("wnearestjet_tau21", &wnearestjet_tau21);
  infile_tree_wjets->SetBranchAddress("wnearestjet_deepak8_WvsQCD", &wnearestjet_deepak8_WvsQCD);
  infile_tree_wjets->SetBranchAddress("wnearestjet_partnet_WvsQCD", &wnearestjet_partnet_WvsQCD);
  infile_tree_wjets->SetBranchAddress("wnearestjet_dr", &wnearestjet_dr);

  const string outfile_path_wjets = sframe_output_path+file_prefix+file_postfix_wjets+".restructured";
  TFile *outfile_wjets = new TFile(outfile_path_wjets.c_str(), "RECREATE");
  TTree *outfile_tree_wjets = new TTree("my_tree", "new flat tree incorporating jets matched to W bosons");

  float w_jet_weight;
  float w_jet_pt;
  float w_jet_msd;
  float w_jet_tau21;
  float w_jet_deepak8_WvsQCD;
  float w_jet_partnet_WvsQCD;
  float w_jet_dr;

  outfile_tree_wjets->Branch("weight", &w_jet_weight, "weight/F");
  outfile_tree_wjets->Branch("pt", &w_jet_pt, "pt/F");
  outfile_tree_wjets->Branch("msd", &w_jet_msd, "msd/F");
  outfile_tree_wjets->Branch("tau21", &w_jet_tau21, "tau21/F");
  outfile_tree_wjets->Branch("deepak8_WvsQCD", &w_jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_wjets->Branch("partnet_WvsQCD", &w_jet_partnet_WvsQCD, "partnet_WvsQCD/F");
  outfile_tree_wjets->Branch("dr", &w_jet_dr, "dr/F");

  const float max_dr = 0.6;

  for(int i = 0; i < infile_tree_wjets->GetEntries(); i++) {
    print_status(i, infile_tree_wjets->GetEntries());
    infile_tree_wjets->GetEntry(i);

    w_jet_weight = event_weight;
    w_jet_pt = wnearestjet_pt;
    w_jet_msd = wnearestjet_msd;
    w_jet_tau21 = wnearestjet_tau21;
    w_jet_deepak8_WvsQCD = wnearestjet_deepak8_WvsQCD;
    w_jet_partnet_WvsQCD = wnearestjet_partnet_WvsQCD;
    w_jet_dr = wnearestjet_dr;
    if(w_jet_dr <= max_dr) outfile_tree_wjets->Fill();
  }

  outfile_tree_wjets->Print();
  outfile_wjets->Write();
}


void restructure_root_trees_ttbar(const string & year) {

  cout << "Working on ttbar" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_ttbar = "TTbarToHadronic_"+year+".root";

  const string infile_path_ttbar = sframe_output_path+file_prefix+file_postfix_ttbar;
  cout << "Open " << infile_path_ttbar << endl;
  TFile *infile_ttbar = TFile::Open(infile_path_ttbar.c_str(), "READ");
  TTree *infile_tree_ttbar = (TTree*)infile_ttbar->Get("AnalysisTree");

  float event_weight;

  float tnearestjet_pt;
  float tnearestjet_msd;
  float tnearestjet_subdeepcsv;
  float tnearestjet_subdeepjet;
  float tnearestjet_tau32;
  float tnearestjet_tau21;
  float tnearestjet_deepak8_TvsQCD;
  float tnearestjet_deepak8_WvsQCD;
  float tnearestjet_partnet_TvsQCD;
  float tnearestjet_partnet_WvsQCD;
  float tnearestjet_dr;

  float antitnearestjet_pt;
  float antitnearestjet_msd;
  float antitnearestjet_subdeepcsv;
  float antitnearestjet_subdeepjet;
  float antitnearestjet_tau32;
  float antitnearestjet_tau21;
  float antitnearestjet_deepak8_TvsQCD;
  float antitnearestjet_deepak8_WvsQCD;
  float antitnearestjet_partnet_TvsQCD;
  float antitnearestjet_partnet_WvsQCD;
  float antitnearestjet_dr;

  bool same_t_jets;

  float wplusnearestjet_pt;
  float wplusnearestjet_msd;
  float wplusnearestjet_subdeepcsv;
  float wplusnearestjet_subdeepjet;
  float wplusnearestjet_tau32;
  float wplusnearestjet_tau21;
  float wplusnearestjet_deepak8_TvsQCD;
  float wplusnearestjet_deepak8_WvsQCD;
  float wplusnearestjet_partnet_TvsQCD;
  float wplusnearestjet_partnet_WvsQCD;
  float wplusnearestjet_dr;
  float wplusnearestjet_dr_b;

  float wminusnearestjet_pt;
  float wminusnearestjet_msd;
  float wminusnearestjet_subdeepcsv;
  float wminusnearestjet_subdeepjet;
  float wminusnearestjet_tau32;
  float wminusnearestjet_tau21;
  float wminusnearestjet_deepak8_TvsQCD;
  float wminusnearestjet_deepak8_WvsQCD;
  float wminusnearestjet_partnet_TvsQCD;
  float wminusnearestjet_partnet_WvsQCD;
  float wminusnearestjet_dr;
  float wminusnearestjet_dr_antib;

  bool same_w_jets;

  infile_tree_ttbar->SetBranchAddress("event_weight", &event_weight);

  infile_tree_ttbar->SetBranchAddress("tnearestjet_pt", &tnearestjet_pt);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_msd", &tnearestjet_msd);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_subjets_deepcsv_max", &tnearestjet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_subjets_deepjet_max", &tnearestjet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_tau32", &tnearestjet_tau32);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_tau21", &tnearestjet_tau21);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_deepak8_TvsQCD", &tnearestjet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_deepak8_WvsQCD", &tnearestjet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_partnet_TvsQCD", &tnearestjet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_partnet_WvsQCD", &tnearestjet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_dr", &tnearestjet_dr);

  infile_tree_ttbar->SetBranchAddress("antitnearestjet_pt", &antitnearestjet_pt);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_msd", &antitnearestjet_msd);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_subjets_deepcsv_max", &antitnearestjet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_subjets_deepjet_max", &antitnearestjet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_tau32", &antitnearestjet_tau32);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_tau21", &antitnearestjet_tau21);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_deepak8_TvsQCD", &antitnearestjet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_deepak8_WvsQCD", &antitnearestjet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_partnet_TvsQCD", &antitnearestjet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_partnet_WvsQCD", &antitnearestjet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_dr", &antitnearestjet_dr);

  infile_tree_ttbar->SetBranchAddress("the_two_t_jets_are_the_same", &same_t_jets);

  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_pt", &wplusnearestjet_pt);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_msd", &wplusnearestjet_msd);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_subjets_deepcsv_max", &wplusnearestjet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_subjets_deepjet_max", &wplusnearestjet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_tau32", &wplusnearestjet_tau32);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_tau21", &wplusnearestjet_tau21);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_deepak8_TvsQCD", &wplusnearestjet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_deepak8_WvsQCD", &wplusnearestjet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_partnet_TvsQCD", &wplusnearestjet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_partnet_WvsQCD", &wplusnearestjet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_dr", &wplusnearestjet_dr);
  infile_tree_ttbar->SetBranchAddress("wplusnearestjet_dr_b", &wplusnearestjet_dr_b);

  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_pt", &wminusnearestjet_pt);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_msd", &wminusnearestjet_msd);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_subjets_deepcsv_max", &wminusnearestjet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_subjets_deepjet_max", &wminusnearestjet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_tau32", &wminusnearestjet_tau32);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_tau21", &wminusnearestjet_tau21);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_deepak8_TvsQCD", &wminusnearestjet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_deepak8_WvsQCD", &wminusnearestjet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_partnet_TvsQCD", &wminusnearestjet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_partnet_WvsQCD", &wminusnearestjet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_dr", &wminusnearestjet_dr);
  infile_tree_ttbar->SetBranchAddress("wminusnearestjet_dr_antib", &wminusnearestjet_dr_antib);

  infile_tree_ttbar->SetBranchAddress("the_two_w_jets_are_the_same", &same_w_jets);

  const string outfile_path_ttbar_t = sframe_output_path+file_prefix+file_postfix_ttbar+".restructured_t";
  TFile *outfile_ttbar_t = new TFile(outfile_path_ttbar_t.c_str(), "RECREATE");
  TTree *outfile_tree_ttbar_t = new TTree("my_tree", "new flat tree incorporating jets matched to top quarks");

  const string outfile_path_ttbar_w = sframe_output_path+file_prefix+file_postfix_ttbar+".restructured_w";
  TFile *outfile_ttbar_w = new TFile(outfile_path_ttbar_w.c_str(), "RECREATE");
  TTree *outfile_tree_ttbar_w = new TTree("my_tree", "new flat tree incorporating jets matched to W bosons");

  float t_jet_weight;
  float t_jet_pt;
  float t_jet_msd;
  float t_jet_subdeepcsv;
  float t_jet_subdeepjet;
  float t_jet_tau32;
  float t_jet_tau21;
  float t_jet_deepak8_TvsQCD;
  float t_jet_deepak8_WvsQCD;
  float t_jet_partnet_TvsQCD;
  float t_jet_partnet_WvsQCD;
  float t_jet_dr;

  float w_jet_weight;
  float w_jet_pt;
  float w_jet_msd;
  float w_jet_subdeepcsv;
  float w_jet_subdeepjet;
  float w_jet_tau32;
  float w_jet_tau21;
  float w_jet_deepak8_TvsQCD;
  float w_jet_deepak8_WvsQCD;
  float w_jet_partnet_TvsQCD;
  float w_jet_partnet_WvsQCD;
  float w_jet_dr;
  float w_jet_dr_b;

  outfile_tree_ttbar_t->Branch("weight", &t_jet_weight, "weight/F");
  outfile_tree_ttbar_t->Branch("pt", &t_jet_pt, "pt/F");
  outfile_tree_ttbar_t->Branch("msd", &t_jet_msd, "msd/F");
  outfile_tree_ttbar_t->Branch("subdeepcsv", &t_jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_ttbar_t->Branch("subdeepjet", &t_jet_subdeepcsv, "subdeepjet/F");
  outfile_tree_ttbar_t->Branch("tau32", &t_jet_tau32, "tau32/F");
  outfile_tree_ttbar_t->Branch("tau21", &t_jet_tau21, "tau21/F");
  outfile_tree_ttbar_t->Branch("deepak8_TvsQCD", &t_jet_deepak8_TvsQCD, "deepak8_TvsQCD/F");
  outfile_tree_ttbar_t->Branch("deepak8_WvsQCD", &t_jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_ttbar_t->Branch("partnet_TvsQCD", &t_jet_partnet_TvsQCD, "partnet_TvsQCD/F");
  outfile_tree_ttbar_t->Branch("partnet_WvsQCD", &t_jet_partnet_WvsQCD, "partnet_WvsQCD/F");
  outfile_tree_ttbar_t->Branch("dr", &t_jet_dr, "dr/F");

  outfile_tree_ttbar_w->Branch("weight", &w_jet_weight, "weight/F");
  outfile_tree_ttbar_w->Branch("pt", &w_jet_pt, "pt/F");
  outfile_tree_ttbar_w->Branch("msd", &w_jet_msd, "msd/F");
  outfile_tree_ttbar_w->Branch("subdeepcsv", &w_jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_ttbar_w->Branch("subdeepjet", &w_jet_subdeepcsv, "subdeepjet/F");
  outfile_tree_ttbar_w->Branch("tau32", &w_jet_tau32, "tau32/F");
  outfile_tree_ttbar_w->Branch("tau21", &w_jet_tau21, "tau21/F");
  outfile_tree_ttbar_w->Branch("deepak8_TvsQCD", &w_jet_deepak8_TvsQCD, "deepak8_TvsQCD/F");
  outfile_tree_ttbar_w->Branch("deepak8_WvsQCD", &w_jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_ttbar_w->Branch("partnet_TvsQCD", &w_jet_partnet_TvsQCD, "partnet_TvsQCD/F");
  outfile_tree_ttbar_w->Branch("partnet_WvsQCD", &w_jet_partnet_WvsQCD, "partnet_WvsQCD/F");
  outfile_tree_ttbar_w->Branch("dr", &w_jet_dr, "dr/F");
  outfile_tree_ttbar_w->Branch("dr_b", &w_jet_dr_b, "dr_b/F");

  const float max_dr = 0.6;

  for(int i = 0; i < infile_tree_ttbar->GetEntries(); i++) {
    print_status(i, infile_tree_ttbar->GetEntries());
    infile_tree_ttbar->GetEntry(i);

    if(!same_t_jets) {
      t_jet_weight = event_weight;
      t_jet_pt = tnearestjet_pt;
      t_jet_msd = tnearestjet_msd;
      t_jet_subdeepcsv = tnearestjet_subdeepcsv;
      t_jet_subdeepjet = tnearestjet_subdeepjet;
      t_jet_tau32 = tnearestjet_tau32;
      t_jet_tau21 = tnearestjet_tau21;
      t_jet_deepak8_TvsQCD = tnearestjet_deepak8_TvsQCD;
      t_jet_deepak8_WvsQCD = tnearestjet_deepak8_WvsQCD;
      t_jet_partnet_TvsQCD = tnearestjet_partnet_TvsQCD;
      t_jet_partnet_WvsQCD = tnearestjet_partnet_WvsQCD;
      t_jet_dr = tnearestjet_dr;
      if(t_jet_dr <= max_dr) outfile_tree_ttbar_t->Fill();

      t_jet_weight = event_weight;
      t_jet_pt = antitnearestjet_pt;
      t_jet_msd = antitnearestjet_msd;
      t_jet_subdeepcsv = antitnearestjet_subdeepcsv;
      t_jet_subdeepjet = antitnearestjet_subdeepjet;
      t_jet_tau32 = antitnearestjet_tau32;
      t_jet_tau21 = antitnearestjet_tau21;
      t_jet_deepak8_TvsQCD = antitnearestjet_deepak8_TvsQCD;
      t_jet_deepak8_WvsQCD = antitnearestjet_deepak8_WvsQCD;
      t_jet_partnet_TvsQCD = antitnearestjet_partnet_TvsQCD;
      t_jet_partnet_WvsQCD = antitnearestjet_partnet_WvsQCD;
      t_jet_dr = antitnearestjet_dr;
      if(t_jet_dr <= max_dr) outfile_tree_ttbar_t->Fill();
    }

    if(!same_w_jets) {
      w_jet_weight = event_weight;
      w_jet_pt = wplusnearestjet_pt;
      w_jet_msd = wplusnearestjet_msd;
      w_jet_subdeepcsv = wplusnearestjet_subdeepcsv;
      w_jet_subdeepjet = wplusnearestjet_subdeepjet;
      w_jet_tau32 = wplusnearestjet_tau32;
      w_jet_tau21 = wplusnearestjet_tau21;
      w_jet_deepak8_TvsQCD = wplusnearestjet_deepak8_TvsQCD;
      w_jet_deepak8_WvsQCD = wplusnearestjet_deepak8_WvsQCD;
      w_jet_partnet_TvsQCD = wplusnearestjet_partnet_TvsQCD;
      w_jet_partnet_WvsQCD = wplusnearestjet_partnet_WvsQCD;
      w_jet_dr = wplusnearestjet_dr;
      w_jet_dr_b = wplusnearestjet_dr_b;
      if(w_jet_dr <= max_dr) outfile_tree_ttbar_w->Fill();

      w_jet_weight = event_weight;
      w_jet_pt = wminusnearestjet_pt;
      w_jet_msd = wminusnearestjet_msd;
      w_jet_subdeepcsv = wminusnearestjet_subdeepcsv;
      w_jet_subdeepjet = wminusnearestjet_subdeepjet;
      w_jet_tau32 = wminusnearestjet_tau32;
      w_jet_tau21 = wminusnearestjet_tau21;
      w_jet_deepak8_TvsQCD = wminusnearestjet_deepak8_TvsQCD;
      w_jet_deepak8_WvsQCD = wminusnearestjet_deepak8_WvsQCD;
      w_jet_partnet_TvsQCD = wminusnearestjet_partnet_TvsQCD;
      w_jet_partnet_WvsQCD = wminusnearestjet_partnet_WvsQCD;
      w_jet_dr = wminusnearestjet_dr;
      w_jet_dr_b = wminusnearestjet_dr_antib;
      if(w_jet_dr <= max_dr) outfile_tree_ttbar_w->Fill();
    }
  }

  outfile_tree_ttbar_t->Print();
  outfile_ttbar_t->Write();
  outfile_tree_ttbar_w->Print();
  outfile_ttbar_w->Write();
}


void restructure_root_trees(const string & year) {

  restructure_root_trees_qcd(year);
  // restructure_root_trees_wjets(year);
  // restructure_root_trees_ttbar(year);
}
