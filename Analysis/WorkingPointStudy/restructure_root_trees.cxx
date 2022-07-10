#include <vector>
#include <string>

void print_status(const unsigned int index, const unsigned int entries) {

  if(index % 5000000 == 0) cout << "Processing entry " << index << " / " << entries << endl;
}

void restructure_root_trees_qcd_ak8(const string & year, const string & ht_cutoff) {

  cout << "Working on QCD" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/nominal/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_qcd = "QCD_HT"+ht_cutoff+"toInf_"+year+".root";

  const string infile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd;
  cout << "Open " << infile_path_qcd << endl;
  TFile *infile_qcd = TFile::Open(infile_path_qcd.c_str(), "READ");
  TTree *infile_tree_qcd = (TTree*)infile_qcd->Get("AnalysisTree");

  float event_weight;
  vector<float> *ak8jets_pt = 0; // sometimes I *slightly* dislike ROOT ... why the heck do you need to write "= 0" here??? https://root-forum.cern.ch/t/segmentation-fault-in-a-program-reading-a-ttree/6770/5
  vector<float> *ak8jets_msd = 0;
  vector<float> *ak8jets_subdeepcsv = 0;
  vector<float> *ak8jets_subdeepjet = 0;
  vector<float> *ak8jets_tau32 = 0;
  vector<float> *ak8jets_tau21 = 0;
  vector<float> *ak8jets_deepak8_TvsQCD = 0;
  vector<float> *ak8jets_deepak8_WvsQCD = 0;
  vector<float> *ak8jets_partnet_TvsQCD = 0;
  vector<float> *ak8jets_partnet_WvsQCD = 0;

  infile_tree_qcd->SetBranchAddress("event_weight", &event_weight);
  infile_tree_qcd->SetBranchAddress("ak8jets_pt", &ak8jets_pt);
  infile_tree_qcd->SetBranchAddress("ak8jets_msd", &ak8jets_msd);
  infile_tree_qcd->SetBranchAddress("ak8jets_subjets_deepcsv_max", &ak8jets_subdeepcsv);
  infile_tree_qcd->SetBranchAddress("ak8jets_subjets_deepjet_max", &ak8jets_subdeepjet);
  infile_tree_qcd->SetBranchAddress("ak8jets_tau32", &ak8jets_tau32);
  infile_tree_qcd->SetBranchAddress("ak8jets_tau21", &ak8jets_tau21);
  infile_tree_qcd->SetBranchAddress("ak8jets_deepak8_TvsQCD", &ak8jets_deepak8_TvsQCD);
  infile_tree_qcd->SetBranchAddress("ak8jets_deepak8_WvsQCD", &ak8jets_deepak8_WvsQCD);
  infile_tree_qcd->SetBranchAddress("ak8jets_partnet_TvsQCD", &ak8jets_partnet_TvsQCD);
  infile_tree_qcd->SetBranchAddress("ak8jets_partnet_WvsQCD", &ak8jets_partnet_WvsQCD);

  const string outfile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd+".restructured.AK8";
  TFile *outfile_qcd = new TFile(outfile_path_qcd.c_str(), "RECREATE");
  TTree *outfile_tree_qcd = new TTree("my_tree", "new flat tree incorporating all AK8 jets");

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
    for(int j = 0; j < ak8jets_pt->size(); j++) { // loop over jets in vector
      jet_weight = event_weight;
      jet_pt = ak8jets_pt->at(j);
      jet_msd = ak8jets_msd->at(j);
      jet_subdeepcsv = ak8jets_subdeepcsv->at(j);
      jet_subdeepjet = ak8jets_subdeepjet->at(j);
      jet_tau32 = ak8jets_tau32->at(j);
      jet_tau21 = ak8jets_tau21->at(j);
      jet_deepak8_TvsQCD = ak8jets_deepak8_TvsQCD->at(j);
      jet_deepak8_WvsQCD = ak8jets_deepak8_WvsQCD->at(j);
      jet_partnet_TvsQCD = ak8jets_partnet_TvsQCD->at(j);
      jet_partnet_WvsQCD = ak8jets_partnet_WvsQCD->at(j);
      outfile_tree_qcd->Fill();
    }
  }

  outfile_tree_qcd->Print();
  outfile_qcd->Write();
}


void restructure_root_trees_qcd_hotvr(const string & year, const string & ht_cutoff) {

  cout << "Working on QCD HOTVR" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/nominal/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_qcd = "QCD_HT"+ht_cutoff+"toInf_"+year+".root";

  const string infile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd;
  cout << "Open " << infile_path_qcd << endl;
  TFile *infile_qcd = TFile::Open(infile_path_qcd.c_str(), "READ");
  TTree *infile_tree_qcd = (TTree*)infile_qcd->Get("AnalysisTree");

  float event_weight;
  vector<float> *jets_reff = 0;
  vector<float> *jets_pt = 0;
  vector<float> *jets_mass = 0;
  vector<int> *jets_nsub = 0;
  vector<float> *jets_mpair = 0;
  vector<float> *jets_fpt1 = 0;
  vector<float> *jets_tau32 = 0;

  infile_tree_qcd->SetBranchAddress("event_weight", &event_weight);
  infile_tree_qcd->SetBranchAddress("hotvrjets_reff", &jets_reff);
  infile_tree_qcd->SetBranchAddress("hotvrjets_pt", &jets_pt);
  infile_tree_qcd->SetBranchAddress("hotvrjets_mass", &jets_mass);
  infile_tree_qcd->SetBranchAddress("hotvrjets_nsubjets", &jets_nsub);
  infile_tree_qcd->SetBranchAddress("hotvrjets_mpair", &jets_mpair);
  infile_tree_qcd->SetBranchAddress("hotvrjets_fpt1", &jets_fpt1);
  infile_tree_qcd->SetBranchAddress("hotvrjets_tau32", &jets_tau32);

  const string outfile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd+".restructured.HOTVR";
  TFile *outfile_qcd = new TFile(outfile_path_qcd.c_str(), "RECREATE");
  TTree *outfile_tree_qcd = new TTree("my_tree", "new flat tree incorporating all HOTVR jets");

  float jet_weight;
  float jet_reff;
  float jet_pt;
  float jet_mass;
  float jet_nsub;
  float jet_mpair;
  float jet_fpt1;
  float jet_tau32;

  outfile_tree_qcd->Branch("weight", &jet_weight, "weight/F");
  outfile_tree_qcd->Branch("reff", &jet_reff, "reff/F");
  outfile_tree_qcd->Branch("pt", &jet_pt, "pt/F");
  outfile_tree_qcd->Branch("mass", &jet_mass, "mass/F");
  outfile_tree_qcd->Branch("nsub", &jet_nsub, "nsub/F");
  outfile_tree_qcd->Branch("mpair", &jet_mpair, "mpair/F");
  outfile_tree_qcd->Branch("fpt1", &jet_fpt1, "fpt1/F");
  outfile_tree_qcd->Branch("tau32", &jet_tau32, "tau32/F");

  for(int i = 0; i < infile_tree_qcd->GetEntries(); i++) {
    print_status(i, infile_tree_qcd->GetEntries());
    infile_tree_qcd->GetEntry(i);
    for(int j = 0; j < jets_pt->size(); j++) { // loop over jets in vector
      jet_weight = event_weight;
      jet_reff = jets_reff->at(j);
      jet_pt = jets_pt->at(j);
      jet_mass = jets_mass->at(j);
      jet_nsub = float(jets_nsub->at(j));
      jet_mpair = jets_mpair->at(j);
      jet_fpt1 = jets_fpt1->at(j);
      jet_tau32 = jets_tau32->at(j);
      outfile_tree_qcd->Fill();
    }
  }

  outfile_tree_qcd->Print();
  outfile_qcd->Write();
}


void restructure_root_trees_wjets_ak8(const string & year, const string & ht_cutoff) {

  cout << "Working on wjets" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/nominal/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_wjets = "WJetsToQQ_HT"+ht_cutoff+"toInf_"+year+".root";

  const string infile_path_wjets = sframe_output_path+file_prefix+file_postfix_wjets;
  cout << "Open " << infile_path_wjets << endl;
  TFile *infile_wjets = TFile::Open(infile_path_wjets.c_str(), "READ");
  TTree *infile_tree_wjets = (TTree*)infile_wjets->Get("AnalysisTree");

  float event_weight;
  float wnearestak8jet_pt;
  float wnearestak8jet_msd;
  float wnearestak8jet_tau21;
  float wnearestak8jet_deepak8_WvsQCD;
  float wnearestak8jet_partnet_WvsQCD;
  float wnearestak8jet_dr;

  infile_tree_wjets->SetBranchAddress("event_weight", &event_weight);
  infile_tree_wjets->SetBranchAddress("wnearestak8jet_pt", &wnearestak8jet_pt);
  infile_tree_wjets->SetBranchAddress("wnearestak8jet_msd", &wnearestak8jet_msd);
  infile_tree_wjets->SetBranchAddress("wnearestak8jet_tau21", &wnearestak8jet_tau21);
  infile_tree_wjets->SetBranchAddress("wnearestak8jet_deepak8_WvsQCD", &wnearestak8jet_deepak8_WvsQCD);
  infile_tree_wjets->SetBranchAddress("wnearestak8jet_partnet_WvsQCD", &wnearestak8jet_partnet_WvsQCD);
  infile_tree_wjets->SetBranchAddress("wnearestak8jet_dr", &wnearestak8jet_dr);

  const string outfile_path_wjets = sframe_output_path+file_prefix+file_postfix_wjets+".restructured.AK8";
  TFile *outfile_wjets = new TFile(outfile_path_wjets.c_str(), "RECREATE");
  TTree *outfile_tree_wjets = new TTree("my_tree", "new flat tree incorporating AK8 jets matched to W bosons");

  float w_ak8jet_weight;
  float w_ak8jet_pt;
  float w_ak8jet_msd;
  float w_ak8jet_tau21;
  float w_ak8jet_deepak8_WvsQCD;
  float w_ak8jet_partnet_WvsQCD;
  float w_ak8jet_dr;

  outfile_tree_wjets->Branch("weight", &w_ak8jet_weight, "weight/F");
  outfile_tree_wjets->Branch("pt", &w_ak8jet_pt, "pt/F");
  outfile_tree_wjets->Branch("msd", &w_ak8jet_msd, "msd/F");
  outfile_tree_wjets->Branch("tau21", &w_ak8jet_tau21, "tau21/F");
  outfile_tree_wjets->Branch("deepak8_WvsQCD", &w_ak8jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_wjets->Branch("partnet_WvsQCD", &w_ak8jet_partnet_WvsQCD, "partnet_WvsQCD/F");
  outfile_tree_wjets->Branch("dr", &w_ak8jet_dr, "dr/F");

  // const float max_dr_ak8 = 0.6;
  const float max_dr_ak8 = 99.;

  for(int i = 0; i < infile_tree_wjets->GetEntries(); i++) {
    print_status(i, infile_tree_wjets->GetEntries());
    infile_tree_wjets->GetEntry(i);

    w_ak8jet_weight = event_weight;
    w_ak8jet_pt = wnearestak8jet_pt;
    w_ak8jet_msd = wnearestak8jet_msd;
    w_ak8jet_tau21 = wnearestak8jet_tau21;
    w_ak8jet_deepak8_WvsQCD = wnearestak8jet_deepak8_WvsQCD;
    w_ak8jet_partnet_WvsQCD = wnearestak8jet_partnet_WvsQCD;
    w_ak8jet_dr = wnearestak8jet_dr;
    if(w_ak8jet_dr <= max_dr_ak8 && w_ak8jet_pt >= 0.f) outfile_tree_wjets->Fill();
  }

  outfile_tree_wjets->Print();
  outfile_wjets->Write();
}


void restructure_root_trees_ttbar(const string & year) {

  cout << "Working on ttbar" << endl;

  const string sframe_output_path = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/nominal/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_ttbar = "TTbarToHadronic_"+year+".root";

  const string infile_path_ttbar = sframe_output_path+file_prefix+file_postfix_ttbar;
  cout << "Open " << infile_path_ttbar << endl;
  TFile *infile_ttbar = TFile::Open(infile_path_ttbar.c_str(), "READ");
  TTree *infile_tree_ttbar = (TTree*)infile_ttbar->Get("AnalysisTree");

  float event_weight;

  float tnearestak8jet_pt;
  float tnearestak8jet_msd;
  float tnearestak8jet_subdeepcsv;
  float tnearestak8jet_subdeepjet;
  float tnearestak8jet_tau32;
  float tnearestak8jet_tau21;
  float tnearestak8jet_deepak8_TvsQCD;
  float tnearestak8jet_deepak8_WvsQCD;
  float tnearestak8jet_partnet_TvsQCD;
  float tnearestak8jet_partnet_WvsQCD;
  float tnearestak8jet_dr;

  float antitnearestak8jet_pt;
  float antitnearestak8jet_msd;
  float antitnearestak8jet_subdeepcsv;
  float antitnearestak8jet_subdeepjet;
  float antitnearestak8jet_tau32;
  float antitnearestak8jet_tau21;
  float antitnearestak8jet_deepak8_TvsQCD;
  float antitnearestak8jet_deepak8_WvsQCD;
  float antitnearestak8jet_partnet_TvsQCD;
  float antitnearestak8jet_partnet_WvsQCD;
  float antitnearestak8jet_dr;

  bool same_t_ak8jets;

  float wplusnearestak8jet_pt;
  float wplusnearestak8jet_msd;
  float wplusnearestak8jet_subdeepcsv;
  float wplusnearestak8jet_subdeepjet;
  float wplusnearestak8jet_tau32;
  float wplusnearestak8jet_tau21;
  float wplusnearestak8jet_deepak8_TvsQCD;
  float wplusnearestak8jet_deepak8_WvsQCD;
  float wplusnearestak8jet_partnet_TvsQCD;
  float wplusnearestak8jet_partnet_WvsQCD;
  float wplusnearestak8jet_dr;
  float wplusnearestak8jet_dr_b;

  float wminusnearestak8jet_pt;
  float wminusnearestak8jet_msd;
  float wminusnearestak8jet_subdeepcsv;
  float wminusnearestak8jet_subdeepjet;
  float wminusnearestak8jet_tau32;
  float wminusnearestak8jet_tau21;
  float wminusnearestak8jet_deepak8_TvsQCD;
  float wminusnearestak8jet_deepak8_WvsQCD;
  float wminusnearestak8jet_partnet_TvsQCD;
  float wminusnearestak8jet_partnet_WvsQCD;
  float wminusnearestak8jet_dr;
  float wminusnearestak8jet_dr_antib;

  bool same_w_ak8jets;

  float tnearesthotvrjet_reff;
  float tnearesthotvrjet_pt;
  float tnearesthotvrjet_mass;
  int tnearesthotvrjet_nsub;
  float tnearesthotvrjet_mpair;
  float tnearesthotvrjet_fpt1;
  float tnearesthotvrjet_tau32;
  float tnearesthotvrjet_dr;

  float antitnearesthotvrjet_reff;
  float antitnearesthotvrjet_pt;
  float antitnearesthotvrjet_mass;
  int antitnearesthotvrjet_nsub;
  float antitnearesthotvrjet_mpair;
  float antitnearesthotvrjet_fpt1;
  float antitnearesthotvrjet_tau32;
  float antitnearesthotvrjet_dr;

  bool same_t_hotvrjets;

  infile_tree_ttbar->SetBranchAddress("event_weight", &event_weight);

  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_pt", &tnearestak8jet_pt);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_msd", &tnearestak8jet_msd);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_subjets_deepcsv_max", &tnearestak8jet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_subjets_deepjet_max", &tnearestak8jet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_tau32", &tnearestak8jet_tau32);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_tau21", &tnearestak8jet_tau21);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_deepak8_TvsQCD", &tnearestak8jet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_deepak8_WvsQCD", &tnearestak8jet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_partnet_TvsQCD", &tnearestak8jet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_partnet_WvsQCD", &tnearestak8jet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("tnearestak8jet_dr", &tnearestak8jet_dr);

  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_pt", &antitnearestak8jet_pt);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_msd", &antitnearestak8jet_msd);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_subjets_deepcsv_max", &antitnearestak8jet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_subjets_deepjet_max", &antitnearestak8jet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_tau32", &antitnearestak8jet_tau32);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_tau21", &antitnearestak8jet_tau21);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_deepak8_TvsQCD", &antitnearestak8jet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_deepak8_WvsQCD", &antitnearestak8jet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_partnet_TvsQCD", &antitnearestak8jet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_partnet_WvsQCD", &antitnearestak8jet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("antitnearestak8jet_dr", &antitnearestak8jet_dr);

  infile_tree_ttbar->SetBranchAddress("the_two_t_ak8jets_are_the_same", &same_t_ak8jets);

  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_pt", &wplusnearestak8jet_pt);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_msd", &wplusnearestak8jet_msd);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_subjets_deepcsv_max", &wplusnearestak8jet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_subjets_deepjet_max", &wplusnearestak8jet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_tau32", &wplusnearestak8jet_tau32);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_tau21", &wplusnearestak8jet_tau21);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_deepak8_TvsQCD", &wplusnearestak8jet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_deepak8_WvsQCD", &wplusnearestak8jet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_partnet_TvsQCD", &wplusnearestak8jet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_partnet_WvsQCD", &wplusnearestak8jet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_dr", &wplusnearestak8jet_dr);
  infile_tree_ttbar->SetBranchAddress("wplusnearestak8jet_dr_b", &wplusnearestak8jet_dr_b);

  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_pt", &wminusnearestak8jet_pt);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_msd", &wminusnearestak8jet_msd);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_subjets_deepcsv_max", &wminusnearestak8jet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_subjets_deepjet_max", &wminusnearestak8jet_subdeepjet);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_tau32", &wminusnearestak8jet_tau32);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_tau21", &wminusnearestak8jet_tau21);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_deepak8_TvsQCD", &wminusnearestak8jet_deepak8_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_deepak8_WvsQCD", &wminusnearestak8jet_deepak8_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_partnet_TvsQCD", &wminusnearestak8jet_partnet_TvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_partnet_WvsQCD", &wminusnearestak8jet_partnet_WvsQCD);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_dr", &wminusnearestak8jet_dr);
  infile_tree_ttbar->SetBranchAddress("wminusnearestak8jet_dr_antib", &wminusnearestak8jet_dr_antib);

  infile_tree_ttbar->SetBranchAddress("the_two_w_ak8jets_are_the_same", &same_w_ak8jets);

  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_reff", &tnearesthotvrjet_reff);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_pt", &tnearesthotvrjet_pt);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_mass", &tnearesthotvrjet_mass);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_nsubjets", &tnearesthotvrjet_nsub);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_mpair", &tnearesthotvrjet_mpair);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_fpt1", &tnearesthotvrjet_fpt1);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_tau32", &tnearesthotvrjet_tau32);
  infile_tree_ttbar->SetBranchAddress("tnearesthotvrjet_dr", &tnearesthotvrjet_dr);

  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_reff", &antitnearesthotvrjet_reff);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_pt", &antitnearesthotvrjet_pt);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_mass", &antitnearesthotvrjet_mass);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_nsubjets", &antitnearesthotvrjet_nsub);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_mpair", &antitnearesthotvrjet_mpair);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_fpt1", &antitnearesthotvrjet_fpt1);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_tau32", &antitnearesthotvrjet_tau32);
  infile_tree_ttbar->SetBranchAddress("antitnearesthotvrjet_dr", &antitnearesthotvrjet_dr);

  infile_tree_ttbar->SetBranchAddress("the_two_t_hotvrjets_are_the_same", &same_t_hotvrjets);

  const string outfile_path_ttbar_t_ak8 = sframe_output_path+file_prefix+file_postfix_ttbar+".restructured_t.AK8";
  TFile *outfile_ttbar_t_ak8 = new TFile(outfile_path_ttbar_t_ak8.c_str(), "RECREATE");
  TTree *outfile_tree_ttbar_t_ak8 = new TTree("my_tree", "new flat tree incorporating AK8 jets matched to top quarks");

  const string outfile_path_ttbar_w_ak8 = sframe_output_path+file_prefix+file_postfix_ttbar+".restructured_w.AK8";
  TFile *outfile_ttbar_w_ak8 = new TFile(outfile_path_ttbar_w_ak8.c_str(), "RECREATE");
  TTree *outfile_tree_ttbar_w_ak8 = new TTree("my_tree", "new flat tree incorporating AK8 jets matched to W bosons");

  const string outfile_path_ttbar_t_hotvr = sframe_output_path+file_prefix+file_postfix_ttbar+".restructured_t.HOTVR";
  TFile *outfile_ttbar_t_hotvr = new TFile(outfile_path_ttbar_t_hotvr.c_str(), "RECREATE");
  TTree *outfile_tree_ttbar_t_hotvr = new TTree("my_tree", "new flat tree incorporating HOTVR jets matched to top quarks");

  float t_ak8jet_weight;
  float t_ak8jet_pt;
  float t_ak8jet_msd;
  float t_ak8jet_subdeepcsv;
  float t_ak8jet_subdeepjet;
  float t_ak8jet_tau32;
  float t_ak8jet_tau21;
  float t_ak8jet_deepak8_TvsQCD;
  float t_ak8jet_deepak8_WvsQCD;
  float t_ak8jet_partnet_TvsQCD;
  float t_ak8jet_partnet_WvsQCD;
  float t_ak8jet_dr;

  float w_ak8jet_weight;
  float w_ak8jet_pt;
  float w_ak8jet_msd;
  float w_ak8jet_subdeepcsv;
  float w_ak8jet_subdeepjet;
  float w_ak8jet_tau32;
  float w_ak8jet_tau21;
  float w_ak8jet_deepak8_TvsQCD;
  float w_ak8jet_deepak8_WvsQCD;
  float w_ak8jet_partnet_TvsQCD;
  float w_ak8jet_partnet_WvsQCD;
  float w_ak8jet_dr;
  float w_ak8jet_dr_b;

  float t_hotvrjet_weight;
  float t_hotvrjet_reff;
  float t_hotvrjet_pt;
  float t_hotvrjet_mass;
  float t_hotvrjet_nsub;
  float t_hotvrjet_mpair;
  float t_hotvrjet_fpt1;
  float t_hotvrjet_tau32;
  float t_hotvrjet_dr;

  outfile_tree_ttbar_t_ak8->Branch("weight", &t_ak8jet_weight, "weight/F");
  outfile_tree_ttbar_t_ak8->Branch("pt", &t_ak8jet_pt, "pt/F");
  outfile_tree_ttbar_t_ak8->Branch("msd", &t_ak8jet_msd, "msd/F");
  outfile_tree_ttbar_t_ak8->Branch("subdeepcsv", &t_ak8jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_ttbar_t_ak8->Branch("subdeepjet", &t_ak8jet_subdeepcsv, "subdeepjet/F");
  outfile_tree_ttbar_t_ak8->Branch("tau32", &t_ak8jet_tau32, "tau32/F");
  outfile_tree_ttbar_t_ak8->Branch("tau21", &t_ak8jet_tau21, "tau21/F");
  outfile_tree_ttbar_t_ak8->Branch("deepak8_TvsQCD", &t_ak8jet_deepak8_TvsQCD, "deepak8_TvsQCD/F");
  outfile_tree_ttbar_t_ak8->Branch("deepak8_WvsQCD", &t_ak8jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_ttbar_t_ak8->Branch("partnet_TvsQCD", &t_ak8jet_partnet_TvsQCD, "partnet_TvsQCD/F");
  outfile_tree_ttbar_t_ak8->Branch("partnet_WvsQCD", &t_ak8jet_partnet_WvsQCD, "partnet_WvsQCD/F");
  outfile_tree_ttbar_t_ak8->Branch("dr", &t_ak8jet_dr, "dr/F");

  outfile_tree_ttbar_w_ak8->Branch("weight", &w_ak8jet_weight, "weight/F");
  outfile_tree_ttbar_w_ak8->Branch("pt", &w_ak8jet_pt, "pt/F");
  outfile_tree_ttbar_w_ak8->Branch("msd", &w_ak8jet_msd, "msd/F");
  outfile_tree_ttbar_w_ak8->Branch("subdeepcsv", &w_ak8jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_ttbar_w_ak8->Branch("subdeepjet", &w_ak8jet_subdeepcsv, "subdeepjet/F");
  outfile_tree_ttbar_w_ak8->Branch("tau32", &w_ak8jet_tau32, "tau32/F");
  outfile_tree_ttbar_w_ak8->Branch("tau21", &w_ak8jet_tau21, "tau21/F");
  outfile_tree_ttbar_w_ak8->Branch("deepak8_TvsQCD", &w_ak8jet_deepak8_TvsQCD, "deepak8_TvsQCD/F");
  outfile_tree_ttbar_w_ak8->Branch("deepak8_WvsQCD", &w_ak8jet_deepak8_WvsQCD, "deepak8_WvsQCD/F");
  outfile_tree_ttbar_w_ak8->Branch("partnet_TvsQCD", &w_ak8jet_partnet_TvsQCD, "partnet_TvsQCD/F");
  outfile_tree_ttbar_w_ak8->Branch("partnet_WvsQCD", &w_ak8jet_partnet_WvsQCD, "partnet_WvsQCD/F");
  outfile_tree_ttbar_w_ak8->Branch("dr", &w_ak8jet_dr, "dr/F");
  outfile_tree_ttbar_w_ak8->Branch("dr_b", &w_ak8jet_dr_b, "dr_b/F");

  outfile_tree_ttbar_t_hotvr->Branch("weight", &t_hotvrjet_weight, "weight/F");
  outfile_tree_ttbar_t_hotvr->Branch("reff", &t_hotvrjet_pt, "reff/F");
  outfile_tree_ttbar_t_hotvr->Branch("pt", &t_hotvrjet_pt, "pt/F");
  outfile_tree_ttbar_t_hotvr->Branch("mass", &t_hotvrjet_mass, "mass/F");
  outfile_tree_ttbar_t_hotvr->Branch("nsub", &t_hotvrjet_nsub, "nsub/F");
  outfile_tree_ttbar_t_hotvr->Branch("mpair", &t_hotvrjet_mpair, "mpair/F");
  outfile_tree_ttbar_t_hotvr->Branch("fpt1", &t_hotvrjet_fpt1, "fpt1/F");
  outfile_tree_ttbar_t_hotvr->Branch("tau32", &t_hotvrjet_tau32, "tau32/F");
  outfile_tree_ttbar_t_hotvr->Branch("dr", &t_hotvrjet_dr, "dr/F");

  // const float max_dr_ak8 = 0.6;
  const float max_dr_ak8 = 99.;
  const float max_dr_hotvr = 99.;

  for(int i = 0; i < infile_tree_ttbar->GetEntries(); i++) {
    print_status(i, infile_tree_ttbar->GetEntries());
    infile_tree_ttbar->GetEntry(i);

    if(!same_t_ak8jets) {
      t_ak8jet_weight = event_weight;
      t_ak8jet_pt = tnearestak8jet_pt;
      t_ak8jet_msd = tnearestak8jet_msd;
      t_ak8jet_subdeepcsv = tnearestak8jet_subdeepcsv;
      t_ak8jet_subdeepjet = tnearestak8jet_subdeepjet;
      t_ak8jet_tau32 = tnearestak8jet_tau32;
      t_ak8jet_tau21 = tnearestak8jet_tau21;
      t_ak8jet_deepak8_TvsQCD = tnearestak8jet_deepak8_TvsQCD;
      t_ak8jet_deepak8_WvsQCD = tnearestak8jet_deepak8_WvsQCD;
      t_ak8jet_partnet_TvsQCD = tnearestak8jet_partnet_TvsQCD;
      t_ak8jet_partnet_WvsQCD = tnearestak8jet_partnet_WvsQCD;
      t_ak8jet_dr = tnearestak8jet_dr;
      if(t_ak8jet_dr <= max_dr_ak8 && t_ak8jet_pt >= 0.f) outfile_tree_ttbar_t_ak8->Fill();

      t_ak8jet_weight = event_weight;
      t_ak8jet_pt = antitnearestak8jet_pt;
      t_ak8jet_msd = antitnearestak8jet_msd;
      t_ak8jet_subdeepcsv = antitnearestak8jet_subdeepcsv;
      t_ak8jet_subdeepjet = antitnearestak8jet_subdeepjet;
      t_ak8jet_tau32 = antitnearestak8jet_tau32;
      t_ak8jet_tau21 = antitnearestak8jet_tau21;
      t_ak8jet_deepak8_TvsQCD = antitnearestak8jet_deepak8_TvsQCD;
      t_ak8jet_deepak8_WvsQCD = antitnearestak8jet_deepak8_WvsQCD;
      t_ak8jet_partnet_TvsQCD = antitnearestak8jet_partnet_TvsQCD;
      t_ak8jet_partnet_WvsQCD = antitnearestak8jet_partnet_WvsQCD;
      t_ak8jet_dr = antitnearestak8jet_dr;
      if(t_ak8jet_dr <= max_dr_ak8 && t_ak8jet_pt >= 0.f) outfile_tree_ttbar_t_ak8->Fill();
    }

    if(!same_w_ak8jets) {
      w_ak8jet_weight = event_weight;
      w_ak8jet_pt = wplusnearestak8jet_pt;
      w_ak8jet_msd = wplusnearestak8jet_msd;
      w_ak8jet_subdeepcsv = wplusnearestak8jet_subdeepcsv;
      w_ak8jet_subdeepjet = wplusnearestak8jet_subdeepjet;
      w_ak8jet_tau32 = wplusnearestak8jet_tau32;
      w_ak8jet_tau21 = wplusnearestak8jet_tau21;
      w_ak8jet_deepak8_TvsQCD = wplusnearestak8jet_deepak8_TvsQCD;
      w_ak8jet_deepak8_WvsQCD = wplusnearestak8jet_deepak8_WvsQCD;
      w_ak8jet_partnet_TvsQCD = wplusnearestak8jet_partnet_TvsQCD;
      w_ak8jet_partnet_WvsQCD = wplusnearestak8jet_partnet_WvsQCD;
      w_ak8jet_dr = wplusnearestak8jet_dr;
      w_ak8jet_dr_b = wplusnearestak8jet_dr_b;
      if(w_ak8jet_dr <= max_dr_ak8 && w_ak8jet_pt >= 0.f) outfile_tree_ttbar_w_ak8->Fill();

      w_ak8jet_weight = event_weight;
      w_ak8jet_pt = wminusnearestak8jet_pt;
      w_ak8jet_msd = wminusnearestak8jet_msd;
      w_ak8jet_subdeepcsv = wminusnearestak8jet_subdeepcsv;
      w_ak8jet_subdeepjet = wminusnearestak8jet_subdeepjet;
      w_ak8jet_tau32 = wminusnearestak8jet_tau32;
      w_ak8jet_tau21 = wminusnearestak8jet_tau21;
      w_ak8jet_deepak8_TvsQCD = wminusnearestak8jet_deepak8_TvsQCD;
      w_ak8jet_deepak8_WvsQCD = wminusnearestak8jet_deepak8_WvsQCD;
      w_ak8jet_partnet_TvsQCD = wminusnearestak8jet_partnet_TvsQCD;
      w_ak8jet_partnet_WvsQCD = wminusnearestak8jet_partnet_WvsQCD;
      w_ak8jet_dr = wminusnearestak8jet_dr;
      w_ak8jet_dr_b = wminusnearestak8jet_dr_antib;
      if(w_ak8jet_dr <= max_dr_ak8 && w_ak8jet_pt >= 0.f) outfile_tree_ttbar_w_ak8->Fill();
    }

    if(!same_t_hotvrjets) {
      t_hotvrjet_weight = event_weight;
      t_hotvrjet_reff = tnearesthotvrjet_reff;
      t_hotvrjet_pt = tnearesthotvrjet_pt;
      t_hotvrjet_mass = tnearesthotvrjet_mass;
      t_hotvrjet_nsub = float(tnearesthotvrjet_nsub);
      t_hotvrjet_mpair = tnearesthotvrjet_mpair;
      t_hotvrjet_fpt1 = tnearesthotvrjet_fpt1;
      t_hotvrjet_tau32 = tnearesthotvrjet_tau32;
      t_hotvrjet_dr = tnearesthotvrjet_dr;
      if(t_hotvrjet_dr <= max_dr_hotvr && t_hotvrjet_pt >= 0.f) outfile_tree_ttbar_t_hotvr->Fill();

      t_hotvrjet_weight = event_weight;
      t_hotvrjet_reff = antitnearesthotvrjet_reff;
      t_hotvrjet_pt = antitnearesthotvrjet_pt;
      t_hotvrjet_mass = antitnearesthotvrjet_mass;
      t_hotvrjet_nsub = float(antitnearesthotvrjet_nsub);
      t_hotvrjet_mpair = antitnearesthotvrjet_mpair;
      t_hotvrjet_fpt1 = antitnearesthotvrjet_fpt1;
      t_hotvrjet_tau32 = antitnearesthotvrjet_tau32;
      t_hotvrjet_dr = antitnearesthotvrjet_dr;
      if(t_hotvrjet_dr <= max_dr_hotvr && t_hotvrjet_pt >= 0.f) outfile_tree_ttbar_t_hotvr->Fill();
    }
  }

  outfile_tree_ttbar_t_ak8->Print();
  outfile_ttbar_t_ak8->Write();
  outfile_tree_ttbar_w_ak8->Print();
  outfile_ttbar_w_ak8->Write();
  outfile_tree_ttbar_t_hotvr->Print();
  outfile_ttbar_t_hotvr->Write();
}


void restructure_root_trees(const string & year, const string & option) {

  if(option == "qcd_ak8_200") restructure_root_trees_qcd_ak8(year, "200");
  if(option == "qcd_ak8_300") restructure_root_trees_qcd_ak8(year, "300");
  if(option == "qcd_hotvr_200") restructure_root_trees_qcd_hotvr(year, "200");
  if(option == "wjets_ak8_200") restructure_root_trees_wjets_ak8(year, "200");
  if(option == "ttbar") restructure_root_trees_ttbar(year);
}
