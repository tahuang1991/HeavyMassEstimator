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
    heavyMassEstimator(bool PUsample, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
        int iterations, const std::string& RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose=0);
    heavyMassEstimator();
    ~heavyMassEstimator();

    void set_inputs(const TLorentzVector& lep1_lorentz, const TLorentzVector& lep2_lorentz, 
        const TLorentzVector& b1jet_lorentz, const TLorentzVector& b2jet_lorentz, 
	const TLorentzVector& totjets_lorentz, const TLorentzVector& met_lorentz, 
        const TLorentzVector* nu1_lorentz, const TLorentzVector* nu2_lorentz, 
        const TLorentzVector* b_genp_lorentz, const TLorentzVector* bbar_genp_lorentz, 
        const TLorentzVector* h2tohh_lorentz, 
        int onshellMarker, 
        bool simulation, 
        int ievent);
    void set_inputs(const TLorentzVector& lep1_lorentz, const TLorentzVector& lep2_lorentz, 
        const TLorentzVector& b1jet_lorentz, const TLorentzVector& b2jet_lorentz, 
	const TLorentzVector& totjets_lorentz, const TLorentzVector& met_lorentz, 
        int ievent);

    private:
      TTree* hmetree_;
      TH1F* heavyMassEstimator_h2Mass_;
      TH1F* heavyMassEstimator_h2Massweight1_;
      TH1F* heavyMassEstimator_h2Massweight4_;
      TFile* file_;
      const TH1F* wmasshist_;
      const TH2F* onoffshellWmass_hist_;
      const TH1F* onshellnupt_hist_;
      const TH1F* bjetrescalec1_hist_;
      const TH1F* bjetrescalec2_hist_;
      TRandom* rnd_;
 
   //runheavyMassEstimator
   public: 
      bool runheavyMassEstimator();
      const TH1F& getheavyMassEstimatorh2();
      const TH1F& getheavyMassEstimatorh2weight1();
      const TH1F& getheavyMassEstimatorh2weight4();
      const TTree* getheavyMassEstimatorTree();
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
      float nu1pt_onshellW(EtaPhi nu1_etaphi, const TLorentzVector& lep1lorentz, float wMass);
      bool  nulorentz_offshellW(const TLorentzVector& jetlorentz, const TLorentzVector& lep1lorentz, 
			        const TLorentzVector& lep2lorentz, TLorentzVector& nu1lorentz, 
 			        TLorentzVector& nu2lorentz, int control, float hMass);
      bool  nulorentz_offshellW(const TVector2& met, const TLorentzVector& lep1lorentz, 
			       const TLorentzVector& lep2lorentz, TLorentzVector& nu1lorentz, 
 			       TLorentzVector& nu2lorentz, int control, float hMass);
      bool checkSolution(const TLorentzVector& jetslorentz,
                         const TLorentzVector& lep1lorentz,
                         const TLorentzVector& lep2lorentz,
                         const TLorentzVector& nu1lorentz, int control, float hMass); 
      bool cutsCheck();
      void assignMuLorentzVec(int control);  
          
    private:
      float onshellWMassRandomWalk(float x0, float step, float random);
      float onshellWMassRandomWalk(float x0, float step, float random, const TH1F* hist);
      float onshellWMassPDF(float wmass);
  
    private:
      const TH1F* readoutonshellWMassPDF();
      const TH1F* readoutoffshellWMassPDF();
      const TH2F* readoutonoffshellWMassPDF();
      const TH1F* readoutonshellnuptPDF();
      const TH1F* readoutbjetrescalec1PDF();
      const TH1F* readoutbjetrescalec2PDF();
      const TH2F* readoutbjetrescalec1c2PDF();
 
    private:
      float weightfromhist(const TH1F* pdf, float x); 
      float weightfromhist(const TH2F* pdf, float x, float y, bool whole=true); 
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
      int iev_;
      int onshellMarker_;
      bool simulation_;
      bool PUsample_;
      int iterations_;
      int seed_;
      std::string RefPDFfile_;
      int verbose_;
      int metcorrection_;
      int bjetrescale_;
      float b1rescalefactor_;
      float b2rescalefactor_;
      float rescalec1_;
      float rescalec2_;
      bool heavyMassEstimatordebug_;   

    private:
      TLorentzVector hme_lep1_lorentz_;
      TLorentzVector hme_lep2_lorentz_;
      TLorentzVector hme_bjets_lorentz_;
      TLorentzVector hme_b1jet_lorentz_;
      TLorentzVector hme_b2jet_lorentz_;
      TLorentzVector hme_totjets_lorentz_;
      TVector2 hmemet_vec2_;

      TLorentzVector nu1_lorentz_true_;
      TLorentzVector nu2_lorentz_true_;
      TLorentzVector onshellW_lorentz_true_;
      TLorentzVector offshellW_lorentz_true_;
      TLorentzVector b1_lorentz_;
      TLorentzVector b2_lorentz_;
      TLorentzVector htoWW_lorentz_true_;
      TLorentzVector htoBB_lorentz_true_;
      TLorentzVector h2tohh_lorentz_true_;
      
      TLorentzVector mu_onshellW_lorentz_;
      TLorentzVector mu_offshellW_lorentz_;
      TLorentzVector jets_lorentz_;
      TVector2 met_vec2_;
      TLorentzVector nu_onshellW_lorentz_;
      TLorentzVector nu_offshellW_lorentz_;
      TLorentzVector offshellW_lorentz_;
      TLorentzVector onshellW_lorentz_;
      TLorentzVector htoWW_lorentz_;
      TLorentzVector htoBB_lorentz_;
      TLorentzVector h2tohh_lorentz_;


    public:
      void setlepton1kinematic(float px, float py, float pz, float E) {hme_lep1_lorentz_.SetPxPyPzE(px, py,pz, E);}
      void setlepton2kinematic(float px, float py, float pz, float E) {hme_lep2_lorentz_.SetPxPyPzE(px, py,pz, E);}
      void setb1jetkinematic(float px, float py, float pz, float E) {hme_b1jet_lorentz_.SetPxPyPzE(px, py,pz, E);}
      void setb2jetkinematic(float px, float py, float pz, float E) {hme_b2jet_lorentz_.SetPxPyPzE(px, py,pz, E);}
      void settotjetskinematic(float px, float py, float pz, float E) {hme_totjets_lorentz_.SetPxPyPzE(px, py,pz, E);}
      void setMET(float px,float py) {hmemet_vec2_.Set(px, py); }


    private:
	
      TLorentzVector ideal_met_lorentz_;
      TLorentzVector h2tohh_expect_lorentz_;
    private:
      //branches
      float eta_mean_;
      float eta_rms_;
      float eta_gen_; 
      float phi_gen_;
      float wmass_gen_;
      float hmass_gen_;
      float metpx_gen_;
      float metpy_gen_;
       
      
      int control_;
      float weight_;
      float weight1_;//extra weight
      float weight2_;//extra weight
      float weight3_;//extra weight
      float weight4_;//extra weight
 
      float mu_onshellW_Eta_;
      float mu_onshellW_Phi_;
      float mu_onshellW_Pt_;
      float mu_onshellW_E_;
      float mu_offshellW_Eta_;
      float mu_offshellW_Phi_;
      float mu_offshellW_Pt_;
      float mu_offshellW_E_;
      float nu_onshellW_Eta_;
      float nu_onshellW_Phi_;
      float nu_onshellW_Pt_;
      float nu_onshellW_E_;
      float nu_offshellW_Eta_;
      float nu_offshellW_Phi_;
      float nu_offshellW_Pt_;
      float nu_offshellW_E_;
      
      float onshellW_Eta_;
      float onshellW_Phi_;
      float onshellW_Pt_;
      float onshellW_E_;
      float onshellW_Mass_;
      float offshellW_Eta_;
      float offshellW_Phi_;
      float offshellW_Pt_;
      float offshellW_E_;
      float offshellW_Mass_;
    
      float b1_Pt_;
      float b1_Px_;
      float b1_Py_;
      float b1_E_;
      float b1_Eta_;
      float b1_Phi_;
      float b1_Mass_;
      float b2_Pt_;
      float b2_Px_;
      float b2_Py_;
      float b2_E_;
      float b2_Eta_;
      float b2_Phi_;
      float b2_Mass_;

      float b1jet_Pt_;
      float b1jet_Px_;
      float b1jet_Py_;
      float b1jet_Energy_;
      float b1jet_Eta_;
      float b1jet_Phi_;
      float b1jet_Mass_;
      float b2jet_Pt_;
      float b2jet_Px_;
      float b2jet_Py_;
      float b2jet_Energy_;
      float b2jet_Eta_;
      float b2jet_Phi_;
      float b2jet_Mass_;
      float b1jet_dR_;
      float b2jet_dR_;
      float b1rescalefactor_true_;
      float b2rescalefactor_true_;
      float rescalec1_true_;
      float rescalec2_true_;
       
      float htoBB_jets_Eta_;
      float htoBB_jets_Phi_;
      float htoBB_jets_Pt_;
      float htoBB_jets_E_;
      float htoBB_jets_Mass_;
      float htoBB_Eta_;
      float htoBB_Phi_;
      float htoBB_Pt_;
      float htoBB_E_;
      float htoBB_Mass_;
      float htoWW_Eta_;
      float htoWW_Phi_;
      float htoWW_Pt_;
      float htoWW_E_;
      float htoWW_Mass_;

      float heavyMassEstimatormet_E_;
      float heavyMassEstimatormet_Phi_;
      float heavyMassEstimatormet_Px_;
      float heavyMassEstimatormet_Py_;
      float ideal_met_E_;
      float ideal_met_Px_;
      float ideal_met_Py_;

      float h2tohh_Eta_;
      float h2tohh_Phi_;
      float h2tohh_Pt_;
      float h2tohh_E_;
      float h2tohh_Mass_;


      float met_;
      float met_phi_;
      float met_px_;
      float met_py_;

      float eta_nuoffshellW_true_;
      float phi_nuoffshellW_true_;
      float pt_nuoffshellW_true_;
      float px_nuoffshellW_true_;
      float py_nuoffshellW_true_;
      float eta_nuonshellW_true_;
      float phi_nuonshellW_true_;
      float pt_nuonshellW_true_;
      float px_nuonshellW_true_;
      float py_nuonshellW_true_;
      float mass_offshellW_true_;
      float mass_onshellW_true_;
      float mass_htoWW_true_;
      float pt_h2tohh_true_;
      float mass_h2tohh_true_;
      float mass_h2_expect_;

};
   






#endif
