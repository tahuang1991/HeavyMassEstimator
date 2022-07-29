using namespace std;
#include "../interface/heavyMassEstimator.h"

#include <iostream>
//ROOT

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

struct event {
    
   TLorentzVector lep1_p4;
   TLorentzVector lep2_p4;
   TLorentzVector jet1_p4;
   TLorentzVector jet2_p4;
   TLorentzVector met_p4;
   TLorentzVector totjets_p4;
   void printparticle(TLorentzVector tmp, string message){ std::cout << message <<" ( "<< tmp.Px()<<", "<< tmp.Py() <<", "<< tmp.Pz() <<", "<< tmp.Energy() <<") "<< std::endl;}
   void print(){
       printparticle(lep1_p4, "lepton1");
       printparticle(lep2_p4, "lepton2");
       printparticle(jet1_p4, "jet1");
       printparticle(jet2_p4, "jet2");
       printparticle(met_p4, "MET");
   }

};


//To compile the file: go to main folder, aka `cd HeavyMassEstimator/ and make `
//To run the 1000 events sample: cd exec/ and ./runHME_ntuple ../data/Radion_M700_1kevents.root
//this file woudl read in the events from TTree in inputfile (../data/Radion_M700_1kevents.root)
//and compute the HME for each events
int main(int argc, char *argv[])
{                 

// agrv[0]: inputfile
// argv[1]: outputfile
  std::string inputfile = "../data/Radion_M700_1kevents.root";
  std::string outputfile = "out_runHME_ntuple.root";
  if (argc < 2){
    std::cerr <<"Number of arguments is "<< argc <<"! Would use default inputfile "<< inputfile << std::endl;
  }else if (argc == 2){
    inputfile = std::string(argv[1]);
  }else if (argc == 3){
    inputfile  = std::string(argv[1]);
    outputfile = std::string(argv[2]);
  }

  std::ifstream infile(inputfile);
  if (not infile.good()) {
    std::cerr <<"input file does not exist: "<< inputfile << std::endl;
    exit(1);
  }

  //std::string treename = "Friends";
  std::string treename = "Events";
  TChain chain(treename.c_str());
  chain.Add(inputfile.c_str());

  int totalEv = chain.GetEntries();
  std::cout <<"inputfile "<< inputfile <<" treename:" << treename <<" Entries "<< totalEv << std::endl;
  std::cout <<"outputfile "<< outputfile << std::endl;
  if (totalEv <= 0){
    std::cerr <<"No events in found in TTree in root file "<< inputfile << std::endl;
    exit(1); 
  }


  TFile* out = new TFile(outputfile.c_str(), "recreate");
  out->cd();
  TTree* tree(chain.CloneTree(0));

  //float lep1_pt, lep1_phi, lep1_eta;
  //float lep2_pt, lep2_phi, lep2_eta;
  //float jet1_pt, jet1_phi, jet1_eta, jet1_E;
  //float jet2_pt, jet2_phi, jet2_eta, jet2_E;
  //float met_pt, met_phi;
  //chain.SetBranchAddress("lep1_pt", &lep1_pt);
  //chain.SetBranchAddress("lep1_phi",&lep1_phi);
  //chain.SetBranchAddress("lep1_eta",&lep1_eta);
  //chain.SetBranchAddress("lep2_pt", &lep2_pt);
  //chain.SetBranchAddress("lep2_phi",&lep2_phi);
  //chain.SetBranchAddress("lep2_eta",&lep2_eta);
  //chain.SetBranchAddress("jet1_pt", &jet1_pt);
  //chain.SetBranchAddress("jet1_phi",&jet1_phi);
  //chain.SetBranchAddress("jet1_eta",&jet1_eta);
  //chain.SetBranchAddress("jet1_E",  &jet1_E);
  //chain.SetBranchAddress("jet2_pt", &jet2_pt);
  //chain.SetBranchAddress("jet2_phi",&jet2_phi);
  //chain.SetBranchAddress("jet2_eta",&jet2_eta);
  //chain.SetBranchAddress("jet2_E",  &jet2_E);
  //chain.SetBranchAddress("MET_pt_nom",  &met_pt);
  //chain.SetBranchAddress("MET_phi_nom",  &met_phi);

  float l1_px, l1_py, l1_pz, l1_E;
  float l2_px, l2_py, l2_pz, l2_E;
  float j1_px, j1_py, j1_pz, j1_E;
  float j2_px, j2_py, j2_pz, j2_E;
  float met_px, met_py, met_pz, met_E;
  chain.SetBranchAddress("l1_px",  &l1_px);
  chain.SetBranchAddress("l1_py",  &l1_py);
  chain.SetBranchAddress("l1_pz",  &l1_pz);
  chain.SetBranchAddress("l1_E",   &l1_E);
  chain.SetBranchAddress("l2_px",  &l2_px);
  chain.SetBranchAddress("l2_py",  &l2_py);
  chain.SetBranchAddress("l2_pz",  &l2_pz);
  chain.SetBranchAddress("l2_E",   &l2_E);
  chain.SetBranchAddress("j1_px",  &j1_px);
  chain.SetBranchAddress("j1_py",  &j1_py);
  chain.SetBranchAddress("j1_pz",  &j1_pz);
  chain.SetBranchAddress("j1_E",   &j1_E);
  chain.SetBranchAddress("j2_px",  &j2_px);
  chain.SetBranchAddress("j2_py",  &j2_py);
  chain.SetBranchAddress("j2_pz",  &j2_pz);
  chain.SetBranchAddress("j2_E",   &j2_E);
  chain.SetBranchAddress("met_px",  &met_px);
  chain.SetBranchAddress("met_py",  &met_py);
  chain.SetBranchAddress("met_pz",  &met_pz);
  chain.SetBranchAddress("met_E",   &met_E);




  float heavyMassEstimator_h2mass_prob;
  float heavyMassEstimator_h2massweight1_prob;
  float heavyMassEstimator_h2massweight4_prob;
  float heavyMassEstimator_h2mass_RMS;
  float heavyMassEstimator_h2massweight1_RMS;
  float heavyMassEstimator_h2massweight4_RMS;

  tree->Branch("hme_h2mass_prob",        &heavyMassEstimator_h2mass_prob,        "heavyMassEstimator_h2mass_prob/F");
  tree->Branch("hme_h2massweight1_prob", &heavyMassEstimator_h2massweight1_prob, "heavyMassEstimator_h2massweight1_prob/F");
  tree->Branch("hme_h2massweight4_prob", &heavyMassEstimator_h2massweight4_prob, "heavyMassEstimator_h2massweight4_prob/F");
  tree->Branch("hme_h2mass_rms",         &heavyMassEstimator_h2mass_RMS,         "heavyMassEstimator_h2mass_RMS/F");
  tree->Branch("hme_h2massweight1_rms",  &heavyMassEstimator_h2massweight1_RMS,  "heavyMassEstimator_h2massweight1_RMS/F");
  tree->Branch("hme_h2massweight2_rms",  &heavyMassEstimator_h2massweight4_RMS,  "heavyMassEstimator_h2massweight4_RMS/F");


  //HME configuration
  bool weightfromonshellnupt_func = false;
  bool weightfromonshellnupt_hist = true;
  bool weightfromonoffshellWmass_hist = true;
  string RefPDFfile = "../data/REFPDFPU40.root";//the root file contains histogram for weighting 
  std::ifstream reffile(RefPDFfile);
  if (not reffile.good()) {
      std::cerr <<"Path to root file with profiled distribution: ../data/REFPDFPU40.root NOT found!" << std::endl;
  
      const char* cmssw_base = std::getenv("CMSSW_BASE");
      if (!cmssw_base) {
          std::cerr << "Error! Environmental variable CMSSW_BASE not set!\n"
                    << "Please run cmsenv first.\n"
                    << "When running without CMSSW, you still need this variable so that the\n"
                    << "data/REFPDFPU40.root can be found.\n";
          exit(1);            
      }
      RefPDFfile = std::string(cmssw_base).append("/src/hhAnalysis/Heavymassestimator/data/REFPDFPU40.root");//the root file contains histogram for weighting 
  }
  
  int nevent = totalEv; //totalEv
  bool PUSample = true;//whether event is from PU sample or not.
  int iterations = 10000;
  bool useMET = true;//use MET or totjets_p4 to estimate kinematic sum of two nuetrino
  int bjetrescaleAlgo = 2;//jet correction
  int metcorrection = 5;//met correction
  std::cout <<"HME config: nevent "<< nevent<<" iterations "<< iterations <<(useMET ? " useMET ":" notuseMET ") 
  << " bjetrescaleAlgo " << bjetrescaleAlgo  << " metcorrection " << metcorrection << (PUSample ? " PU40 ": " PU0 ") << std::endl;

  //float h2tohhmass = 400.0;//signal benmark M=400, narrow width. in other word, HME output of all above events should be close to 400.0
  
  for (int ievent=0; ievent < nevent; ievent++){
    chain.GetEntry(ievent);
    event thiseve;
    //thiseve.lep1_p4.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, 0);
    //thiseve.lep2_p4.SetPtEtaPhiM(lep2_pt, lep2_eta, lep2_phi, 0);
    //thiseve.jet1_p4.SetPtEtaPhiE(jet1_pt, jet1_eta, jet1_phi, jet1_E);
    //thiseve.jet2_p4.SetPtEtaPhiE(jet2_pt, jet2_eta, jet2_phi, jet2_E);
    //  thiseve.met_p4.SetPtEtaPhiE(met_pt, 0.0, met_phi, met_pt);
    thiseve.lep1_p4.SetPxPyPzE(l1_px, l1_py, l1_pz, l1_E);
    thiseve.lep2_p4.SetPxPyPzE(l2_px, l2_py, l2_pz, l2_E);
    thiseve.jet1_p4.SetPxPyPzE(j1_px, j1_py, j1_pz, j1_E);
    thiseve.jet2_p4.SetPxPyPzE(j2_px, j2_py, j2_pz, j2_E);
     thiseve.met_p4.SetPxPyPzE(met_px, met_py, met_pz, met_E);
    heavyMassEstimator hme(PUSample, weightfromonshellnupt_func, weightfromonshellnupt_hist, weightfromonoffshellWmass_hist,
    iterations, RefPDFfile, useMET, bjetrescaleAlgo, metcorrection);
    hme.set_inputs(thiseve.lep1_p4, thiseve.lep2_p4, 
      thiseve.jet1_p4, thiseve.jet2_p4, 
      thiseve.totjets_p4, thiseve.met_p4, 
      ievent);
    if (ievent == 1 ){
      //use MET covariance matrix smearing 
      hme.setMETCovMatrix(20.0, 20.0, 10.0, true);
    }
    bool runheavyMassEstimatorok = hme.runheavyMassEstimator();
    if (runheavyMassEstimatorok) {
      //heavyMassEstimatortree =  (hme.getheavyMassEstimatorTree())->CloneTree();
        std::stringstream ss;
      ss <<"heavyMassEstimator_h2mass_event"<< ievent;
      const std::string histname(ss.str());
      std::stringstream ss1;
      ss1 <<"heavyMassEstimator_h2massweight1_event"<< ievent;
      const std::string histname1(ss1.str());
      std::stringstream ss4;
      ss4 <<"heavyMassEstimator_h2massweight4_event"<< ievent;
      const std::string histname4(ss4.str());
      TH1F* heavyMassEstimator_h2mass =(TH1F*)(hme.getheavyMassEstimatorh2()).Clone(histname.c_str());
      TH1F* heavyMassEstimator_h2mass_weight1 =(TH1F*)(hme.getheavyMassEstimatorh2weight1()).Clone(histname1.c_str());
      TH1F* heavyMassEstimator_h2mass_weight4 =(TH1F*)(hme.getheavyMassEstimatorh2weight4()).Clone(histname4.c_str());
      //std::cout <<" Mass_h2mass in Analyzer " << std::endl;
      heavyMassEstimator_h2mass_prob = (heavyMassEstimator_h2mass->GetXaxis())->GetBinCenter(heavyMassEstimator_h2mass->GetMaximumBin());
      heavyMassEstimator_h2massweight1_prob = (heavyMassEstimator_h2mass_weight1->GetXaxis())->GetBinCenter(heavyMassEstimator_h2mass_weight1->GetMaximumBin());
      heavyMassEstimator_h2massweight4_prob = (heavyMassEstimator_h2mass_weight4->GetXaxis())->GetBinCenter(heavyMassEstimator_h2mass_weight4->GetMaximumBin());
      heavyMassEstimator_h2mass_RMS = heavyMassEstimator_h2mass->GetRMS();
      heavyMassEstimator_h2massweight1_RMS = heavyMassEstimator_h2mass_weight1->GetRMS();
      heavyMassEstimator_h2massweight4_RMS = heavyMassEstimator_h2mass_weight4->GetRMS();
      //float heavyMassEstimator_h2mass_Entries = heavyMassEstimator_h2mass->GetEntries();
      //float heavyMassEstimator_h2mass_Mean = heavyMassEstimator_h2mass->GetMean();
      //int nbin=(heavyMassEstimator_h2mass->GetXaxis())->GetNbins();
      //float heavyMassEstimator_h2mass_overflow = heavyMassEstimator_h2mass->GetBinContent(nbin+1);
      //float heavyMassEstimator_h2mass_underflow = heavyMassEstimator_h2mass->GetBinContent(-1);
      //float heavyMassEstimator_h2mass_MaxBin = heavyMassEstimator_h2mass->GetBinContent( heavyMassEstimator_h2mass->GetMaximumBin() );
      //float heavyMassEstimator_h2mass_weight1_MaxBin = heavyMassEstimator_h2mass_weight1->GetBinContent( heavyMassEstimator_h2mass_weight1->GetMaximumBin() );
      //float heavyMassEstimator_h2mass_weight4_MaxBin = heavyMassEstimator_h2mass_weight4->GetBinContent( heavyMassEstimator_h2mass_weight4->GetMaximumBin() );
      //std::cout <<"True HH mass "<< h2tohh_mass <<"; Reconstructed HH mass " << heavyMassEstimator_h2mass_prob <<" +/- "<< heavyMassEstimator_h2mass_RMS << "; reconstructed HH mass with type1 weight "<< heavyMassEstimator_h2massweight1_prob <<" +/- " << heavyMassEstimator_h2massweight1_RMS <<"; Reconstructed HH mass with type2 weight "<< heavyMassEstimator_h2massweight4_prob <<" +/- "<< heavyMassEstimator_h2massweight4_RMS <<std::endl;
      std::cout <<"iEvent "<< ievent <<" Reconstructed HH mass " << heavyMassEstimator_h2mass_prob <<" +/- "<< heavyMassEstimator_h2mass_RMS << "; reconstructed HH mass with type1 weight "<< heavyMassEstimator_h2massweight1_prob <<" +/- " << heavyMassEstimator_h2massweight1_RMS <<"; Reconstructed HH mass with type2 weight "<< heavyMassEstimator_h2massweight4_prob <<" +/- "<< heavyMassEstimator_h2massweight4_RMS <<std::endl;
      if (heavyMassEstimator_h2mass_prob < 245 ) {
        std::cerr <<" error !!! heavyMassEstimator_h2mass_prob < 250, heavyMassEstimator_h2mass_prob: " <<heavyMassEstimator_h2mass_prob 
          <<" eventid "<< ievent <<std::endl;
        heavyMassEstimator_h2mass->Print("ALL");
      }
      if (heavyMassEstimator_h2massweight1_prob < 245 ) {
        std::cerr <<" error !!! heavyMassEstimator_h2mass_prob < 250, heavyMassEstimator_h2massweight1_prob: " <<heavyMassEstimator_h2massweight1_prob 
        <<" eventid "<< ievent <<std::endl;
        heavyMassEstimator_h2mass_weight1->Print("ALL");
      }
      //##!! NOT ALL HISTOS 
      //if (keepheavyMassEstimatorhist)
      //	file->WriteObject(heavyMassEstimator_h2mass, histname.c_str());
    }//runheavyMassEstimatorok
    else{
      //no HME is found for this event, fill default values
      heavyMassEstimator_h2mass_prob  = 0.0;
      heavyMassEstimator_h2massweight1_prob  = 0.0;
      heavyMassEstimator_h2massweight4_prob  = 0.0;
      heavyMassEstimator_h2mass_RMS  = -1.0;
      heavyMassEstimator_h2massweight1_RMS  = -1.0;
      heavyMassEstimator_h2massweight4_RMS  = -1.0;
    }
    tree->Fill();
  }
  //tree->Write();
  out->WriteObject(tree, treename.c_str());
  out->Close();
  std::cout << std::endl;
  std::cout << "*****************************************************************************************************************************************" << std::endl;
  std::cout << "if you want to use this code, please cite:                                                                                               " << std::endl;
  std::cout << "T. Huang, J. M. No, L. Pernié, M. Ramsey-Musolf, A. Safonov, M. Spannowsky, and P. Winslow                                               " << std::endl;
  std::cout << "\" Resonant di-Higgs boson production in the bbWW channel: Probing the electroweak phase transition at the LHC\"                         " << std::endl;
  std::cout << "Phys. Rev. D 96, 035007 – Published 11 August 2017                                                                                       " << std::endl;
  std::cout << "*****************************************************************************************************************************************" << std::endl;
  std:: cout << std::endl;

}
