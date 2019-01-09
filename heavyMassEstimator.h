// -*- C++ -*-
//
//// Package:    DiHiggsWWAnalyzer
//// class heavyMassEstimator
//
// Description: [one line class summary]
// use to run heavyMassEstimator algorithm. 
//  Implementation:
//       [Notes on implementation]
//       */
//       //
//       // Original Author:  tao huang
//       //         Created:  Thurs, 04 06 2015
//       // $Id$
//       //
//       //

#ifndef heavyMassEstimator_h
#define heavyMassEstimator_h

//std lib
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"  


//root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"


typedef std::pair<float, float> EtaPhi;

class heavyMassEstimator{

    public:
    //constructor
    heavyMassEstimator(TLorentzVector* lep1_lorentz, TLorentzVector* lep2_lorentz, TLorentzVector* b1jet_lorentz, TLorentzVector* b2jet_lorentz, 
	TLorentzVector* totjets_lorentz,TLorentzVector* met_lorentz, TLorentzVector* nu1_lorentz, TLorentzVector* nu2_lorentz,
	TLorentzVector* b_genp_lorentz, TLorentzVector* bbar_genp_lorentz, TLorentzVector* h2tohh_lorentz, int onshellMarker, bool simulation,	       bool PUsample_,
	int ievent, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
        int iterations, std::string RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose_=0
	);
    //constructor 
    heavyMassEstimator(TLorentzVector* lep1_lorentz, TLorentzVector* lep2_lorentz, TLorentzVector* b1jet_lorentz, TLorentzVector* b2jet_lorentz, 
	TLorentzVector* totjets_lorentz,TLorentzVector* met_lorentz,
	bool PUsample_, int ievent, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
        int iterations, std::string RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose_=0
	);
    heavyMassEstimator();
    ~heavyMassEstimator();

    private:
      TTree *hmetree;
      TH1F heavyMassEstimator_h2Mass;
      TH1F heavyMassEstimator_h2Massweight1;
      TH1F heavyMassEstimator_h2Massweight4;
      TFile *file;
 
   //runheavyMassEstimator
   public: 
      bool runheavyMassEstimator();
      TH1F getheavyMassEstimatorh2();
      TH1F getheavyMassEstimatorh2weight1();
      TH1F getheavyMassEstimatorh2weight4();
      TTree* getheavyMassEstimatorTree();
      //TH1F* getheavyMassEstimatorNeutrio_onshell1();
      //TH1F* getheavyMassEstimatorNeutrio_onshell2();
      //TH1F* getheavyMassEstimatorNeutrio_offshell1();
      //TH1F* getheavyMassEstimatorNeutrio_offshell2();
   private:
      void metCorrection();
      bool bjetsCorrection();

   private:
      void initTree(TTree* hmetree);
        
     
      float genEtaGuass(float mean, float rms);
      float genPhiFlat();
      EtaPhi generatenu1_etaphi();
      float nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* lep1lorentz, float wMass);
      bool  nulorentz_offshellW(TLorentzVector* jetlorentz, TLorentzVector* lep1lorentz, 
			       TLorentzVector* lep2lorentz, TLorentzVector* nu1lorentz, 
 			       TLorentzVector* nu2lorentz, int control, float hMass);
      bool  nulorentz_offshellW(TVector2* met, TLorentzVector* lep1lorentz, 
			       TLorentzVector* lep2lorentz, TLorentzVector* nu1lorentz, 
 			       TLorentzVector* nu2lorentz, int control, float hMass);
      bool checkSolution(TLorentzVector* jetslorentz,
                          TLorentzVector* lep1lorentz,
                          TLorentzVector* lep2lorentz,
                          TLorentzVector* nu1lorentz, int control, float hMass); 
      bool cutsCheck();
      void assignMuLorentzVec(int control);  
          
    private:
      float onshellWMassRandomWalk(float x0, float step, float random);
      float onshellWMassRandomWalk(float x0, float step, float random, TH1F* hist);
      float onshellWMassPDF(float wmass);
  
    private:
      TH1F* readoutonshellWMassPDF();
      TH1F* readoutoffshellWMassPDF();
      TH2F* readoutonoffshellWMassPDF();
      TH1F* readoutonshellnuptPDF();
      TH1F* readoutbjetrescalec1PDF();
      TH1F* readoutbjetrescalec2PDF();
      TH2F* readoutbjetrescalec1c2PDF();
 
    private:
      float weightfromhist(TH1F* pdf, float x); 
      float weightfromhist(TH2F* pdf, float x, float y, bool whole=true); 
      float weightfromonshellnupt(float nupt); 
   
    private:
      bool weightfromonshellnupt_func_;
      bool weightfromonshellnupt_hist_;
      bool weightfromonoffshellWmass_hist_;
      bool weightfrombjetrescalec1c2_hist_;
      bool useMET_;   
      bool writehmetree_;

    private:
      TLorentzVector calculateMET(); 
    public:
      void printTrueLorentz();
      void printheavyMassEstimatorresult(); 

    private:
      int iev;
      int onshellMarker;
      bool simulation;
      bool PUsample;
      int iterations_;
      int seed_;
      std::string RefPDFfile_;
      int verbose;
      int metcorrection_;
      int bjetrescale_;
      float b1rescalefactor;
      float b2rescalefactor;
      float rescalec1;
      float rescalec2;
      bool heavyMassEstimatordebug;   

    private:
      TLorentzVector* hme_lep1_lorentz;
      TLorentzVector* hme_lep2_lorentz;
      TLorentzVector* hme_bjets_lorentz;
      TLorentzVector* hme_b1jet_lorentz;
      TLorentzVector* hme_b2jet_lorentz;
      TLorentzVector* hme_totjets_lorentz;
      TVector2* hmemet_vec2;

      TLorentzVector* nu1_lorentz_true;
      TLorentzVector* nu2_lorentz_true;
      TLorentzVector* onshellW_lorentz_true;
      TLorentzVector* offshellW_lorentz_true;
      TLorentzVector* b1_lorentz;
      TLorentzVector* b2_lorentz;
      TLorentzVector* htoWW_lorentz_true;
      TLorentzVector* htoBB_lorentz_true;
      TLorentzVector* h2tohh_lorentz_true;
      
      TLorentzVector* mu_onshellW_lorentz;
      TLorentzVector* mu_offshellW_lorentz;
      TLorentzVector* jets_lorentz;
      TVector2* met_vec2;
      TLorentzVector* nu_onshellW_lorentz;
      TLorentzVector* nu_offshellW_lorentz;
      TLorentzVector* offshellW_lorentz;
      TLorentzVector* onshellW_lorentz;
      TLorentzVector* htoWW_lorentz;
      TLorentzVector* htoBB_lorentz;
      TLorentzVector* h2tohh_lorentz;


    public:
      void setlepton1kinematic(float px, float py, float pz, float E) {hme_lep1_lorentz->SetPxPyPzE(px, py,pz, E);}
      void setlepton2kinematic(float px, float py, float pz, float E) {hme_lep2_lorentz->SetPxPyPzE(px, py,pz, E);}
      void   setb1jetkinematic(float px, float py, float pz, float E) {hme_b1jet_lorentz->SetPxPyPzE(px, py,pz, E);}
      void   setb2jetkinematic(float px, float py, float pz, float E) {hme_b2jet_lorentz->SetPxPyPzE(px, py,pz, E);}
      void settotjetskinematic(float px, float py, float pz, float E) {hme_totjets_lorentz->SetPxPyPzE(px, py,pz, E);}
      void setMET(float px,float py) {hmemet_vec2->Set(px, py); }


    private:
	
      TLorentzVector ideal_met_lorentz;
      TLorentzVector h2tohh_expect_lorentz;
    private:
      //branches
      float eta_mean;
      float eta_rms;
      float eta_gen; 
      float phi_gen;
      float wmass_gen;
      float hmass_gen;
      float metpx_gen;
      float metpy_gen;
       
      
      int control;
      float weight;
      float weight1;//extra weight
      float weight2;//extra weight
      float weight3;//extra weight
      float weight4;//extra weight
 
      float mu_onshellW_Eta;
      float mu_onshellW_Phi;
      float mu_onshellW_Pt;
      float mu_onshellW_E;
      float mu_offshellW_Eta;
      float mu_offshellW_Phi;
      float mu_offshellW_Pt;
      float mu_offshellW_E;
      float nu_onshellW_Eta;
      float nu_onshellW_Phi;
      float nu_onshellW_Pt;
      float nu_onshellW_E;
      float nu_offshellW_Eta;
      float nu_offshellW_Phi;
      float nu_offshellW_Pt;
      float nu_offshellW_E;
      
      float onshellW_Eta;
      float onshellW_Phi;
      float onshellW_Pt;
      float onshellW_E;
      float onshellW_Mass;
      float offshellW_Eta;
      float offshellW_Phi;
      float offshellW_Pt;
      float offshellW_E;
      float offshellW_Mass;
    
      float b1_Pt;
      float b1_Px;
      float b1_Py;
      float b1_E;
      float b1_Eta;
      float b1_Phi;
      float b1_Mass;
      float b2_Pt;
      float b2_Px;
      float b2_Py;
      float b2_E;
      float b2_Eta;
      float b2_Phi;
      float b2_Mass;

      float b1jet_Pt;
      float b1jet_Px;
      float b1jet_Py;
      float b1jet_Energy;
      float b1jet_Eta;
      float b1jet_Phi;
      float b1jet_Mass;
      float b2jet_Pt;
      float b2jet_Px;
      float b2jet_Py;
      float b2jet_Energy;
      float b2jet_Eta;
      float b2jet_Phi;
      float b2jet_Mass;
      float b1jet_dR;
      float b2jet_dR;
      float b1rescalefactor_true;
      float b2rescalefactor_true;
      float rescalec1_true;
      float rescalec2_true;
       
      float htoBB_jets_Eta;
      float htoBB_jets_Phi;
      float htoBB_jets_Pt;
      float htoBB_jets_E;
      float htoBB_jets_Mass;
      float htoBB_Eta;
      float htoBB_Phi;
      float htoBB_Pt;
      float htoBB_E;
      float htoBB_Mass;
      float htoWW_Eta;
      float htoWW_Phi;
      float htoWW_Pt;
      float htoWW_E;
      float htoWW_Mass;

      float heavyMassEstimatormet_E;
      float heavyMassEstimatormet_Phi;
      float heavyMassEstimatormet_Px;
      float heavyMassEstimatormet_Py;
      float ideal_met_E;
      float ideal_met_Px;
      float ideal_met_Py;

      float h2tohh_Eta;
      float h2tohh_Phi;
      float h2tohh_Pt;
      float h2tohh_E;
      float h2tohh_Mass;


      float met;
      float met_phi;
      float met_px;
      float met_py;

      float eta_nuoffshellW_true;
      float phi_nuoffshellW_true;
      float pt_nuoffshellW_true;
      float px_nuoffshellW_true;
      float py_nuoffshellW_true;
      float eta_nuonshellW_true;
      float phi_nuonshellW_true;
      float pt_nuonshellW_true;
      float px_nuonshellW_true;
      float py_nuonshellW_true;
      float mass_offshellW_true;
      float mass_onshellW_true;
      float mass_htoWW_true;
      float pt_h2tohh_true;
      float mass_h2tohh_true;
      float mass_h2_expect;

};
   






#endif
