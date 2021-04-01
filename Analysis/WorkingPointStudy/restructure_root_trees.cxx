#include <vector>
#include <string>

void restructure_root_trees_qcd() {

  const string sframe_output_path = "/nfs/dust/cms/user/matthies/LegacyTopTagging/WorkingPointStudy/UL17/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_qcd = "QCD_HT300toInf_UL17.root";

  const string infile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd;
  TFile *infile_qcd = TFile::Open(infile_path_qcd.c_str(), "READ");
  TTree *infile_tree_qcd = (TTree*)infile_qcd->Get("AnalysisTree");

  float event_weight;
  vector<float> *jets_pt;
  vector<float> *jets_msd;
  vector<float> *jets_subdeepcsv;
  vector<float> *jets_tau32;

  infile_tree_qcd->SetBranchAddress("event_weight", &event_weight);
  infile_tree_qcd->SetBranchAddress("jets_pt", &jets_pt);
  infile_tree_qcd->SetBranchAddress("jets_msd", &jets_msd);
  infile_tree_qcd->SetBranchAddress("jets_subjets_deepcsv_max", &jets_subdeepcsv);
  infile_tree_qcd->SetBranchAddress("jets_tau32", &jets_tau32);

  const string outfile_path_qcd = sframe_output_path+file_prefix+file_postfix_qcd+".restructured";
  TFile *outfile_qcd = new TFile(outfile_path_qcd.c_str(), "RECREATE");
  TTree *outfile_tree_qcd = new TTree("my_tree", "new flat tree incorporating all QCD jets");

  float jet_weight;
  float jet_pt;
  float jet_msd;
  float jet_subdeepcsv;
  float jet_tau32;

  // outfile_tree_ttbar->Branch("jet", &jet, "weight:pt:msd:subdeepcsv:tau32:dr/F");
  outfile_tree_qcd->Branch("weight", &jet_weight, "weight/F");
  outfile_tree_qcd->Branch("pt", &jet_pt, "pt/F");
  outfile_tree_qcd->Branch("msd", &jet_msd, "msd/F");
  outfile_tree_qcd->Branch("subdeepcsv", &jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_qcd->Branch("tau32", &jet_tau32, "tau32/F");

  for(int i = 0; i < infile_tree_qcd->GetEntries(); i++) {
    infile_tree_qcd->GetEntry(i);
    for(int j = 0; j < jets_pt->size(); j++) { // loop over jets in vector
      jet_weight = event_weight;
      jet_pt = jets_pt->at(j);
      jet_msd = jets_msd->at(j);
      jet_subdeepcsv = jets_subdeepcsv->at(j);
      jet_tau32 = jets_tau32->at(j);
      outfile_tree_qcd->Fill();
    }
  }

  outfile_tree_qcd->Print();
  outfile_qcd->Write();
}


void restructure_root_trees_ttbar() {

  const string sframe_output_path = "/nfs/dust/cms/user/matthies/LegacyTopTagging/WorkingPointStudy/UL17/";
  const string file_prefix = "uhh2.AnalysisModuleRunner.MC.";
  const string file_postfix_ttbar = "TTbarToHadronic_UL17.root";

  const string infile_path_ttbar = sframe_output_path+file_prefix+file_postfix_ttbar;
  TFile *infile_ttbar = TFile::Open(infile_path_ttbar.c_str(), "READ");
  TTree *infile_tree_ttbar = (TTree*)infile_ttbar->Get("AnalysisTree");

  float event_weight;
  bool same_jets;
  float tnearestjet_pt;
  float tnearestjet_msd;
  float tnearestjet_subdeepcsv;
  float tnearestjet_tau32;
  float tnearestjet_dr;
  float antitnearestjet_pt;
  float antitnearestjet_msd;
  float antitnearestjet_subdeepcsv;
  float antitnearestjet_tau32;
  float antitnearestjet_dr;

  infile_tree_ttbar->SetBranchAddress("event_weight", &event_weight);
  infile_tree_ttbar->SetBranchAddress("the_two_jets_are_the_same", &same_jets);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_pt", &tnearestjet_pt);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_msd", &tnearestjet_msd);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_subjets_deepcsv_max", &tnearestjet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_tau32", &tnearestjet_tau32);
  infile_tree_ttbar->SetBranchAddress("tnearestjet_dr", &tnearestjet_dr);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_pt", &antitnearestjet_pt);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_msd", &antitnearestjet_msd);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_subjets_deepcsv_max", &antitnearestjet_subdeepcsv);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_tau32", &antitnearestjet_tau32);
  infile_tree_ttbar->SetBranchAddress("antitnearestjet_dr", &antitnearestjet_dr);

  const string outfile_path_ttbar = sframe_output_path+file_prefix+file_postfix_ttbar+".restructured";
  TFile *outfile_ttbar = new TFile(outfile_path_ttbar.c_str(), "RECREATE");
  TTree *outfile_tree_ttbar = new TTree("my_tree", "new flat tree incorporating jets matched to top quarks");

  // typedef struct {
  //   float weight;
  //   float pt;
  //   float msd;
  //   float subdeepcsv;
  //   float tau32;
  //   float dr;
  // } Jet;
  // Jet jet;

  float jet_weight;
  float jet_pt;
  float jet_msd;
  float jet_subdeepcsv;
  float jet_tau32;
  float jet_dr;

  // outfile_tree_ttbar->Branch("jet", &jet, "weight:pt:msd:subdeepcsv:tau32:dr/F");
  outfile_tree_ttbar->Branch("weight", &jet_weight, "weight/F");
  outfile_tree_ttbar->Branch("pt", &jet_pt, "pt/F");
  outfile_tree_ttbar->Branch("msd", &jet_msd, "msd/F");
  outfile_tree_ttbar->Branch("subdeepcsv", &jet_subdeepcsv, "subdeepcsv/F");
  outfile_tree_ttbar->Branch("tau32", &jet_tau32, "tau32/F");
  outfile_tree_ttbar->Branch("dr", &jet_dr, "dr/F");

  float max_dr = 0.6;

  for(int i = 0; i < infile_tree_ttbar->GetEntries(); i++) {
    infile_tree_ttbar->GetEntry(i);
    if(same_jets) continue;

     jet_weight = event_weight;
     jet_pt = tnearestjet_pt;
     jet_msd = tnearestjet_msd;
     jet_subdeepcsv = tnearestjet_subdeepcsv;
     jet_tau32 = tnearestjet_tau32;
     jet_dr = tnearestjet_dr;
    if( jet_dr <= max_dr) {
      outfile_tree_ttbar->Fill();
    }

     jet_weight = event_weight;
     jet_pt = antitnearestjet_pt;
     jet_msd = antitnearestjet_msd;
     jet_subdeepcsv = antitnearestjet_subdeepcsv;
     jet_tau32 = antitnearestjet_tau32;
     jet_dr = antitnearestjet_dr;
    if( jet_dr <= max_dr) {
      outfile_tree_ttbar->Fill();
    }
  }

  outfile_tree_ttbar->Print();
  outfile_ttbar->Write();
}


void restructure_root_trees() {

  restructure_root_trees_qcd();
  // restructure_root_trees_ttbar();
}
