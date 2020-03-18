// -*- C++ -*-
//

#ifndef heavyMassEstimator_cc
#define heavyMassEstimator_cc



//Simulation or not
#include <iostream>
#include "../interface/heavyMassEstimator.h"
#include <time.h>       /* time */

//constructor
heavyMassEstimator::heavyMassEstimator(bool PUsample, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
        int iterations, const std::string& RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose)
  : hmetree_(nullptr)
  , heavyMassEstimator_h2Mass_(nullptr)
  , heavyMassEstimator_h2Massweight1_(nullptr)
  , heavyMassEstimator_h2Massweight4_(nullptr)
  , file_(nullptr)
  , wmasshist_(nullptr)
  , onoffshellWmass_hist_(nullptr)
  , onshellnupt_hist_(nullptr)
  , bjetrescalec1_hist_(nullptr)
  , bjetrescalec2_hist_(nullptr)
  , rnd_(nullptr)
{

   heavyMassEstimatordebug_=false;
   if (heavyMassEstimatordebug_) std::cout <<" heavyMassEstimator::Debug::1 "<< std::endl;

   PUsample_ = PUsample;

   verbose_ = verbose; 

   weightfromonshellnupt_func_ = weightfromonshellnupt_func;
   weightfromonshellnupt_hist_ = weightfromonshellnupt_hist;
   weightfromonoffshellWmass_hist_ = weightfromonoffshellWmass_hist;
   iterations_ = iterations;
   RefPDFfile_ = RefPDFfile;
   useMET_ = useMET;
  ///0.no correction; 1. simple rescale, 2.elaborate rescale, -1.ideal cas, 4. simple rescale from bjet simple rescale. 5: elaborate rescale from bjet elaborate rescale 
   metcorrection_ = metcorrection;//0.no correction; 1. simple rescale, 2.elaborate rescale, -1.ideal case
   bjetrescale_ = bjetrescaleAlgo;//0.no rescale; 1.simple rescale; 2.elaborate rescale, -1.ideal case
   writehmetree_ = false;
   
   std::stringstream histss;
   std::stringstream histweight1ss;
   std::stringstream histweight4ss;
   histss << "heavyMassEstimator_h2Mass_" << iev_;
   histweight1ss << "heavyMassEstimator_h2Mass_weight1_" << iev_;
   histweight4ss << "heavyMassEstimator_h2Mass_weight4_" << iev_;
   const std::string histname(histss.str());
   const std::string histweight1name(histweight1ss.str());
   const std::string histweight4name(histweight4ss.str());
   heavyMassEstimator_h2Mass_ = new TH1F(histname.c_str(),histname.c_str(), 3800, 200, 4000);
   heavyMassEstimator_h2Massweight1_ = new TH1F(histweight1name.c_str(),histweight1name.c_str(), 3800, 200, 4000);
   heavyMassEstimator_h2Massweight4_ = new TH1F(histweight4name.c_str(),histweight4name.c_str(), 3800, 200, 4000);
   
   file_ = TFile::Open(RefPDFfile_.c_str(),"READ");
   wmasshist_ = readoutonshellWMassPDF(); 
   onoffshellWmass_hist_ = readoutonoffshellWMassPDF(); 
   onshellnupt_hist_ = readoutonshellnuptPDF(); 
   bjetrescalec1_hist_ = readoutbjetrescalec1PDF(); 
   bjetrescalec2_hist_ = readoutbjetrescalec2PDF(); 
   //std::cout <<" rescale priori distribution 1" << std::endl;
   //std::cout <<" onshellnupt max content " <<onshellnupt_hist_->GetBinContent(onshellnupt_hist_->GetMaximumBin()) << std::endl;
   (const_cast<TH1F*>(onshellnupt_hist_))->Scale(1.0/onshellnupt_hist_->GetBinContent(onshellnupt_hist_->GetMaximumBin()));
   (const_cast<TH2F*>(onoffshellWmass_hist_))->Scale(1.0/onoffshellWmass_hist_->GetBinContent(onoffshellWmass_hist_->GetMaximumBin()));
   (const_cast<TH1F*>(bjetrescalec1_hist_))->Scale(1.0/bjetrescalec1_hist_->GetBinContent(bjetrescalec1_hist_->GetMaximumBin()));
   (const_cast<TH1F*>(bjetrescalec2_hist_))->Scale(1.0/bjetrescalec2_hist_->GetBinContent(bjetrescalec2_hist_->GetMaximumBin()));
   //replace above by normalization
   //(const_cast<TH1F*>(onshellnupt_hist_))->Scale(1.0/onshellnupt_hist->GetBinContent(onshellnupt_hist_->GetMaximumBin()));
   //(const_cast<TH1F*>(onoffshellWmass_hist_))->Scale(1.0/onoffshellWmass_hist->GetBinContent(onoffshellWmass_hist_->GetMaximumBin()));
   //(const_cast<TH1F*>(bjetrescalec2_hist_))->Scale(1.0/bjetrescalec2_hist_->GetBinContent(bjetrescalec2_hist_->GetMaximumBin()));
   //std::cout <<" rescale priori distribution 2" << std::endl;
   //file_->Close();

   rnd_ = new TRandom3();

   if (heavyMassEstimatordebug_) std::cout <<" heavyMassEstimator::Debug::2 "<< std::endl;
}

heavyMassEstimator::heavyMassEstimator(){

  //std::cout <<" empty constructor " << std::endl;

}


//deconstructor
heavyMassEstimator::~heavyMassEstimator(){

  delete hmetree_;
  delete heavyMassEstimator_h2Mass_;
  delete heavyMassEstimator_h2Massweight1_;
  delete heavyMassEstimator_h2Massweight4_;
  delete wmasshist_;
  delete onoffshellWmass_hist_;
  delete onshellnupt_hist_;
  delete bjetrescalec1_hist_;
  delete bjetrescalec2_hist_;
  delete file_;
  delete rnd_;

}

void 
heavyMassEstimator::set_inputs(const TLorentzVector& lep1_lorentz, const TLorentzVector& lep2_lorentz, 
        const TLorentzVector& b1jet_lorentz, const TLorentzVector& b2jet_lorentz,
	const TLorentzVector& totjets_lorentz, const TLorentzVector& met_lorentz, 
        int ievent)
{

  set_inputs(lep1_lorentz, lep2_lorentz, 
    b1jet_lorentz, b2jet_lorentz, 
    totjets_lorentz, met_lorentz, 
    NULL, NULL, 
    NULL, NULL, 
    NULL, 
    -1, 
    false,
    ievent);

}

void 
heavyMassEstimator::set_inputs(const TLorentzVector& lep1_lorentz, const TLorentzVector& lep2_lorentz, 
        const TLorentzVector& b1jet_lorentz, const TLorentzVector& b2jet_lorentz,
	const TLorentzVector& totjets_lorentz, const TLorentzVector& met_lorentz, 
        const TLorentzVector* nu1_lorentz, const TLorentzVector* nu2_lorentz,
	const TLorentzVector* b_genp_lorentz, const TLorentzVector* bbar_genp_lorentz, 
        const TLorentzVector* h2_lorentz, 
        int onshellMarker, // simulation only
        bool simulation, 
        int ievent)
{
  hme_lep1_lorentz_ = lep1_lorentz;
  hme_lep2_lorentz_ = lep2_lorentz;
  hme_b1jet_lorentz_ = b1jet_lorentz;
  hme_b2jet_lorentz_ = b2jet_lorentz;
  hme_totjets_lorentz_ = totjets_lorentz;
  hme_bjets_lorentz_ = b1jet_lorentz + b2jet_lorentz;
  hmemet_vec2_ = TVector2(met_lorentz.Px(),met_lorentz.Py());

  //std::cout <<" lep1 lorentz "; lep1_lorentz.Print();
  //std::cout <<" b1 jet lorentz "; b1jet_lorentz.Print();
  //std::cout <<" b1 gen lorentz "; b_genp_lorentz.Print();
  simulation_ = simulation;

  onshellMarker_ = onshellMarker;
  if (simulation_){
    assert(nu1_lorentz && nu2_lorentz && b_genp_lorentz && bbar_genp_lorentz && h2_lorentz); 
    nu1_lorentz_true_ = *nu1_lorentz;
    nu2_lorentz_true_ = *nu2_lorentz;
        
    if (onshellMarker_ == 1){
      onshellW_lorentz_true_ = nu1_lorentz_true_ + hme_lep1_lorentz_;
      offshellW_lorentz_true_ = nu2_lorentz_true_ + hme_lep2_lorentz_;
    }
    else if (onshellMarker_ ==2 ){
      offshellW_lorentz_true_ = nu1_lorentz_true_ + hme_lep1_lorentz_;
      onshellW_lorentz_true_ = nu2_lorentz_true_ + hme_lep2_lorentz_;
    } 
    else std::cout <<" onshellMarker input error"  << std::endl;

    htoWW_lorentz_true_ = offshellW_lorentz_true_ + onshellW_lorentz_true_;
    h2tohh_lorentz_true_ = *h2_lorentz; 
    b1_lorentz_ = *b_genp_lorentz;
    b2_lorentz_ = *bbar_genp_lorentz;
    htoBB_lorentz_true_ = *b_genp_lorentz + *bbar_genp_lorentz;
    ideal_met_lorentz_.SetXYZM(nu1_lorentz_true_.Px()+nu2_lorentz_true_.Px(),nu1_lorentz_true_.Py()+nu2_lorentz_true_.Py(),0,0);
  }
   
  b1rescalefactor_ = 1;
  b2rescalefactor_ = 1;
  rescalec1_ = 1;
  rescalec2_ = 1;

  nu_onshellW_lorentz_ = TLorentzVector();
  nu_offshellW_lorentz_ = TLorentzVector();
  offshellW_lorentz_ = TLorentzVector();
  onshellW_lorentz_ = TLorentzVector();
  htoWW_lorentz_ = TLorentzVector();
  htoBB_lorentz_ = hme_bjets_lorentz_;
  h2tohh_lorentz_ = TLorentzVector();
  met_vec2_ = TVector2(hmemet_vec2_.Px(), hmemet_vec2_.Py());
  //printTrueLorentz(); 

  iev_ = ievent;

  heavyMassEstimator_h2Mass_->Reset();
  heavyMassEstimator_h2Massweight1_->Reset();
  heavyMassEstimator_h2Massweight4_->Reset();

  if ( writehmetree_ ) {
    delete hmetree_;
    std::stringstream ss;
    ss << "hmetree_" << iev_;
    const std::string name(ss.str());
    //hmetree_ = fs->make<TTree>(name.c_str(), name.c_str());
    hmetree_ = new TTree(name.c_str(), name.c_str());
    initTree(hmetree_);
  }

  runheavyMassEstimator();
}


//================================================================================================================
//
// runheavyMassEstimator algo
//
//
//================================================================================================================
//----------- method called to run heavyMassEstimator method for each case -------------------
// control 0 : take muon from onshellW as muon from onshell W and nu_offshellW_eta = some_eta+deltaeta
// control 1 : take muon from onshellW  as muon from onshell W and nu_offshellW_eta = some_eta-deltaeta
// control 2 : take muon from offshellW as muon from onshell W and nu_offshellW_eta = some_eta+deltaeta
// control 3 : take muon from offshellW as muon from onshell W and nu_offshellW_eta = some_eta-deltaeta
bool 
heavyMassEstimator::runheavyMassEstimator(){//should not include any gen level information here in final version

  // genetated (eta,phi) pair
  eta_gen_ = 0;
  phi_gen_ = 0;
  //PU0:14.8,  PU40:25.2
  float met_sigma = (PUsample_? 25.2:14.8);
  //std::cout <<(PUsample_?" PUsample ":"Not PUsample ")<< " met_sigma "<< met_sigma << std::endl;
  //if (writehmetree_) hmetree_->SetDebug(100,0,9999999);
  //int count = 100000;
  bool validrun = false;
  eta_mean_=0;
  eta_rms_=1.403;
  //std::cout <<" time null " << time(NULL) << std::endl;
  seed_ = time(NULL);
  rnd_->SetSeed(seed_+iev_);
  //TF1* wmasspdf = new TF1("wmasspdf","exp(x*7.87e-3+1.69)+603.47*exp(-0.5*((x-80.1)/2.0)**2)",50,90);

   if (heavyMassEstimatordebug_) std::cout <<" heavyMassEstimator::Debug::3  start runheavyMassEstimator() in heavyMassEstimator class "  << std::endl; 
   float nu_onshellW_pt =0;
   wmass_gen_ = 80.3;// initial value
   float step,random01;
  // printTrueLorentz();

  for (int i = 0; i < iterations_ ; i++){
    eta_gen_ = rnd_->Uniform(-6,6); 
    phi_gen_ = rnd_->Uniform(-3.1415926, 3.1415926);
    //wmass_gen_ = rnd_->Gaus(80.385,0.015);
    hmass_gen_ = rnd_->Gaus(125.03,0.3);
    if (metcorrection_>3){
    	metpx_gen_ = rnd_->Gaus(0.0,met_sigma);
   	metpy_gen_ = rnd_->Gaus(0.0,met_sigma);
    }else {
	metpx_gen_ = 0;
   	metpy_gen_ = 0;
	}
    TVector2 met_gen = TVector2(metpx_gen_, metpy_gen_);

    //generate onshell Wmass
    step = rnd_->Uniform(-4,4);
    //step = rnd_->Gaus(0,8);
    random01 = rnd_->Uniform(0,1);
    //wmass_gen_ = onshellWMassRandomWalk(wmass_gen_, step, random01);
    wmass_gen_ = onshellWMassRandomWalk(wmass_gen_, step, random01, wmasshist_);
    if (bjetrescale_ ==1){
	//type1 bjet correction
	b1rescalefactor_ = 125/hme_bjets_lorentz_.M();
	b2rescalefactor_ = 125/hme_bjets_lorentz_.M();
	rescalec1_ = b1rescalefactor_; 
        rescalec2_ = b2rescalefactor_;
    }
    if (bjetrescale_ ==2){
	//type2 bjet correction
	rescalec1_ = bjetrescalec1_hist_->GetRandom();
	//std::cout <<" rescale c1 " << rescalec1 << std::endl;
	bool hascorrection =  bjetsCorrection();
	if (not hascorrection) continue;
	//htoBB_lorentz_ = b1rescalefactor_*b1lorentz_ + b2rescalefactor_*b2lorentz_; 
    } 
    if (bjetrescale_>0){
	htoBB_lorentz_ = b1rescalefactor_*hme_b1jet_lorentz_ + b2rescalefactor_*hme_b2jet_lorentz_;
	/*if (b2rescalefactor_ > 4 or b2rescalefactor_ < 0.1){
		//continue;
        	std::cerr <<" heavyMassEstimator bjetrescale b1 "<< b1rescalefactor_ << " b2 "<< b2rescalefactor_ << std::endl;
		//std::cerr <<" htobb after correction mass "<< htoBB_lorentz_.M(); htoBB_lorentz_.Print();
	 }*/
        if (fabs(htoBB_lorentz_.M()-125)>1 && verbose_ > 0){
		continue;
		std::cerr <<" error the htobb mass is not close 125, htobb_mass "<< htoBB_lorentz_.M() << std::endl;
        	std::cerr <<" heavyMassEstimator bjetrescale b1 "<< b1rescalefactor_ << " b2 "<< b2rescalefactor_ << std::endl;
	 }
		
	}

    //metcorrection_ == 4 or 5 and metcorrection_ == bjetrescale_+3, do bjet correction and then propagate correction to met
    if ((metcorrection_-3)>0 and  ((metcorrection_ -3)== bjetrescale_ or metcorrection_==bjetrescale_)) metCorrection();
    else if ((metcorrection_-3) ==1 or metcorrection_==1){
	//no bjet correction but does correct MET according to type1 bjet correction
	b1rescalefactor_ = 125/hme_bjets_lorentz_.M();
	b2rescalefactor_ = 125/hme_bjets_lorentz_.M();
	metCorrection(); 
    }
    else if ((metcorrection_-3) ==2 or metcorrection_==2){
	//no bjet correction but does correct MET according to type2 bjet correction
	rescalec1_ = bjetrescalec1_hist_->GetRandom();
	bool hascorrection = bjetsCorrection();//calculate b1rescalefactor_ b2rescalefactor_
	if (not hascorrection) continue;
	metCorrection(); 
    }
    else if (metcorrection_ > 5){
	//no bjet correction. only smearing MET according to MET resolution
        b1rescalefactor_ = 1.0;
        b2rescalefactor_ = 1.0;
	//bjetsCorrection();
	metCorrection(); 
    }

    //std::cout <<" Met input px "<< hmemet_vec2_.Px() << " py "<< hmemet_vec2_.Py() <<" pt "<< hmemet_vec2_.Mod() <<std::endl;
    //std::cout <<" met before smearing metpx_gen " << met_vec2_.Px() <<" metpy_gen " << met_vec2_.Py() <<std::endl;
    //std::cout <<" metpx_gen " << metpx_gen_ <<" metpy_gen " << metpy_gen_ <<std::endl;
    if (metcorrection_>3 and useMET_) met_vec2_ = met_vec2_+met_gen;
    if (heavyMassEstimatordebug_) std::cout <<" heavyMassEstimator::Debug::4 in heavyMassEstimator loop met for heavyMassEstimator px " << met_vec2_.Px() <<" py " << met_vec2_.Py() <<std::endl;
    if (bjetrescale_ == -1 && simulation_)
	htoBB_lorentz_ = htoBB_lorentz_true_;
    if (metcorrection_ ==-1 && simulation_)
	met_vec2_ = TVector2(nu1_lorentz_true_.Px()+nu2_lorentz_true_.Px(),nu1_lorentz_true_.Py()+nu2_lorentz_true_.Py());
    /*
    std::cout <<" heavyMassEstimator metCorrection b1 "<< b1rescalefactor_ << " b2 "<< b2rescalefactor_ << std::endl;
    std::cout <<"b1jet Px "<<hme_b1jet_lorentz_.Px() <<" Py "<<hme_b1jet_lorentz_.Py() 
	<<" b2jet Px "<< hme_b2jet_lorentz_.Px() <<" Py "<<hme_b2jet_lorentz_.Py() << std::endl;
    std::cout <<"sum of two nu  px "<<nu1_lorentz_true_.Px()+nu2_lorentz_true_.Px() <<" py "<<nu1_lorentz_true_.Py()+nu2_lorentz_true_.Py()<<std::endl;
    std::cout <<" heavyMassEstimator input met  px "<< hmemet_vec2_.Px() << " py "<<hmemet_vec2_.Py() <<" pt "<< hmemet_vec2_.Mod() <<std::endl;
    std::cout <<" met after correction px "<< met_vec2_.Px() << " py "<< met_vec2_.Py() <<" pt "<< met_vec2_.Mod() <<std::endl;
    // std::cout <<" bjets input M_h= "<< htoBB_lorentz_.M(); htoBB_lorentz_.Print();
    //wmass_gen_ = wmasspdf->GetRandom(50.0,90.0);
    //test
    //eta_gen_ = eta_nuonshellW_true_;
    //phi_gen_ = phi_nuonshellW_true_;
    //wmass_gen_ = mass_onshellW_true_; 
    std::cout << "true eta phi of nuonshell ("<<eta_nuonshellW_true_ <<","<<phi_nuonshellW_true_<<"), pt " <<pt_nuonshellW_true_ 
	<<" mass of onshellW " << mass_onshellW_true_ <<" wmass_gen "<< wmass_gen_  << std::endl;
	std::cout << "true eta phi of nuoffshell ("<<eta_nuoffshellW_true_ <<","<<phi_nuoffshellW_true_<<"), pt " <<pt_nuoffshellW_true_ 
	<<" mass of offshellW " << mass_offshellW_true_ <<  std::endl;
     */
    //std::cout <<" eta_gen "<< eta_gen_ << " phi_gen " << phi_gen_ << " wmass "<<wmass_gen_ << std::endl;
    int solutions = 0;//count num of soluble case
    bool solution[4]={false, false, false, false}; //record whether the case is soluble or not
    for (int j = 0; j < 4; j++){
	assignMuLorentzVec(j/2);
	nu_onshellW_pt = nu1pt_onshellW(std::make_pair(eta_gen_, phi_gen_), mu_onshellW_lorentz_, wmass_gen_); 
	nu_onshellW_lorentz_.SetPtEtaPhiM(nu_onshellW_pt, eta_gen_, phi_gen_,0);
	//solution[j] = nulorentz_offshellW(jets_lorentz_, mu_onshellW_lorentz_,
	//std::cout << " calculate nu1_pt " << nu_onshellW_pt << " px " << nu_onshellW_lorentz_.Px() <<" py "<< nu_onshellW_lorentz_.Py() << std::endl;
	if (useMET_)
	  solution[j] = nulorentz_offshellW(met_vec2_, mu_onshellW_lorentz_,
		mu_offshellW_lorentz_, nu_onshellW_lorentz_,
		nu_offshellW_lorentz_, j%2, hmass_gen_);
	else
	  solution[j] = nulorentz_offshellW(hme_totjets_lorentz_, mu_onshellW_lorentz_,
		mu_offshellW_lorentz_, nu_onshellW_lorentz_,
		nu_offshellW_lorentz_, j%2, hmass_gen_);

	//std::cout << j << " nu_offshellW_eta " << nu_offshellW_lorentz_.Eta()<<" phi " << nu_offshellW_lorentz_.Phi() << std::endl; 
	if (solution[j]) solutions++;
    }
    //       nu_offshellW_lorentz_= NULL; 
    for (int j = 0; j < 4; j++){
	//if ( writehmetree_ ) hmetree_->Fill();
	if (!solution[j])  continue;
	// reassign muons LorentzVector
	if (simulation_){
	  TLorentzVector tmpnu12;
	  tmpnu12.SetXYZM(met_vec2_.Px(),met_vec2_.Py(),nu1_lorentz_true_.Pz()+nu2_lorentz_true_.Pz(),0);
	  h2tohh_expect_lorentz_ = mu_onshellW_lorentz_ + mu_offshellW_lorentz_ + htoBB_lorentz_ + tmpnu12;
	  mass_h2_expect_ = h2tohh_expect_lorentz_.M();
	}

	control_ = j;
	assignMuLorentzVec(j/2);
	nu_onshellW_pt = nu1pt_onshellW(std::make_pair(eta_gen_, phi_gen_), mu_onshellW_lorentz_, wmass_gen_); 
	nu_onshellW_lorentz_.SetPtEtaPhiM(nu_onshellW_pt, eta_gen_, phi_gen_,0);
	//nulorentz_offshellW(jets_lorentz_, mu_onshellW_lorentz_,
	if (useMET_)
	  nulorentz_offshellW(met_vec2_, mu_onshellW_lorentz_,
		mu_offshellW_lorentz_, nu_onshellW_lorentz_,
		nu_offshellW_lorentz_, j%2, hmass_gen_);
	else
	  nulorentz_offshellW(hme_totjets_lorentz_, mu_onshellW_lorentz_,
		mu_offshellW_lorentz_, nu_onshellW_lorentz_,
		nu_offshellW_lorentz_, j%2, hmass_gen_);

	weight_ = 1.0/solutions;// change weight if we consider possibility factor  like matrix elements
	mu_onshellW_Eta_ = mu_onshellW_lorentz_.Eta();
	mu_onshellW_Phi_ = mu_onshellW_lorentz_.Phi();
	mu_onshellW_Pt_ = mu_onshellW_lorentz_.Pt();
	mu_onshellW_E_ = mu_onshellW_lorentz_.E();

	mu_offshellW_Eta_ = mu_offshellW_lorentz_.Eta();
	mu_offshellW_Phi_ = mu_offshellW_lorentz_.Phi();
	mu_offshellW_Pt_ = mu_offshellW_lorentz_.Pt();
	mu_offshellW_E_ = mu_offshellW_lorentz_.E();

	nu_onshellW_Eta_ = nu_onshellW_lorentz_.Eta();
	nu_onshellW_Phi_ = nu_onshellW_lorentz_.Phi();
	nu_onshellW_Pt_ = nu_onshellW_lorentz_.Pt();
	nu_onshellW_E_ = nu_onshellW_lorentz_.E();

	nu_offshellW_Eta_ = nu_offshellW_lorentz_.Eta();
	nu_offshellW_Phi_ = nu_offshellW_lorentz_.Phi();
	nu_offshellW_Pt_ = nu_offshellW_lorentz_.Pt();
	nu_offshellW_E_ = nu_offshellW_lorentz_.E();

	onshellW_lorentz_ = mu_onshellW_lorentz_ + nu_onshellW_lorentz_;
	offshellW_lorentz_ = mu_offshellW_lorentz_ + nu_offshellW_lorentz_;
	htoWW_lorentz_ = onshellW_lorentz_ + offshellW_lorentz_;
	h2tohh_lorentz_ = htoWW_lorentz_ + htoBB_lorentz_;
	if (h2tohh_lorentz_.M()<245 or h2tohh_lorentz_.M()>3800) {
                 if (verbose_ > 0) {
			std::cerr <<" heavyMassEstimator h2 mass is too small, or too large,  M_h " <<h2tohh_lorentz_.M() << std::endl;
			std::cerr <<" gen nu eta "<< eta_gen_ <<" nu phi "<< phi_gen_ << std::endl;
			std::cerr <<" from heavyMassEstimator mu_onshell (px,py,pz, E)= ("<< mu_onshellW_lorentz_.Px()<<", "<<  mu_onshellW_lorentz_.Py()<<", "<< mu_onshellW_lorentz_.Pz()<<", "<< mu_onshellW_lorentz_.E() <<")"<< std::endl;
			std::cerr <<" from heavyMassEstimator mu_offshell (px,py,pz, E)= ("<< mu_offshellW_lorentz_.Px()<<", "<<  mu_offshellW_lorentz_.Py()<<", "<< mu_offshellW_lorentz_.Pz()<<", "<< mu_offshellW_lorentz_.E() <<")"<< std::endl;
			std::cerr <<" from heavyMassEstimator nu_onshell (px,py,pz, E)= ("<< nu_onshellW_lorentz_.Px()<<", "<<  nu_onshellW_lorentz_.Py()<<", "<< nu_onshellW_lorentz_.Pz()<<", "<< nu_onshellW_lorentz_.E() <<")"<< std::endl;
			std::cerr <<" from heavyMassEstimator nu_offshell (px,py,pz, E)= ("<< nu_offshellW_lorentz_.Px()<<", "<<  nu_offshellW_lorentz_.Py()<<", "<< nu_offshellW_lorentz_.Pz()<<", "<< nu_offshellW_lorentz_.E() <<")"<< std::endl;
			std::cerr <<" from heavyMassEstimator htoBB, mass "<< htoBB_lorentz_.M()<<"(px,py,pz, E)= ("<<htoBB_lorentz_.Px()<<", "<< htoBB_lorentz_.Py() <<", "<< htoBB_lorentz_.Pz() <<", "<< htoBB_lorentz_.E()<<")" <<std::endl;
                 }
                if (simulation_ && verbose_ > 0) {
    			std::cerr <<"following is pure gen level infromation " << std::endl;
    			std::cerr <<" nu1 px "<<nu1_lorentz_true_.Px() << " py " <<nu1_lorentz_true_.Py() << " pt "<< nu1_lorentz_true_.Pt() 
			<< " eta "<<nu1_lorentz_true_.Eta() << " phi "<< nu1_lorentz_true_.Phi() << std::endl;
    			std::cerr <<" nu2 px "<<nu2_lorentz_true_.Px() << " py " <<nu2_lorentz_true_.Py() << " pt "<< nu2_lorentz_true_.Pt() 
			<< " eta "<<nu2_lorentz_true_.Eta() << " phi "<< nu2_lorentz_true_.Phi() << std::endl;
    			std::cerr <<" onshellW mass "<< onshellW_lorentz_true_.M(); onshellW_lorentz_true_.Print();  
    			std::cerr <<"offshellW mass " <<offshellW_lorentz_true_.M(); offshellW_lorentz_true_.Print();  
    			std::cerr <<" htoWW mass "<< htoWW_lorentz_true_.M(); htoWW_lorentz_true_.Print();
    			std::cerr <<" htoBB mass "<< htoBB_lorentz_true_.M(); htoBB_lorentz_true_.Print();
    			std::cerr <<" h2tohh, pz " <<h2tohh_lorentz_true_.Pz() << " Mass " << h2tohh_lorentz_true_.M() << std::endl;
   		}	
		
		continue;
	}
	//*h2tohh_lorentz_ = *htoWW_lorentz_+*htoBB_lorentz_true_;

	//*met_vec2_ = TVector2(nu_onshellW_lorentz_.Px()+nu_offshellW_lorentz_.Px(),
	//				nu_onshellW_lorentz_.Py()+nu_offshellW_lorentz_.Py());
        if (fabs(hmass_gen_-htoWW_lorentz_.M()) > 2 && verbose_ > 0) {
	  std::cout << "  hmass_gen " << hmass_gen_ << " Higgs mass from heavyMassEstimator " << htoWW_lorentz_.M() <<std::endl;
          //verbose_ = 4;
	}

	if (verbose_ > 3){
	  std::cout << " onshell W mass "<< onshellW_lorentz_.M();   onshellW_lorentz_.Print();
	  std::cout << " offshell W mass "<< offshellW_lorentz_.M(); offshellW_lorentz_.Print();
	  std::cout << " htoWW mass "<< htoWW_lorentz_.M(); htoWW_lorentz_.Print();
	  //std::cout << " htoBB mass "<< htoBB_lorentz_.M(); htoBB_lorentz_.Print();
	  std::cout << " htoBB mass "<< htoBB_lorentz_true_.M(); htoBB_lorentz_true_.Print();
          //verbose_ = 0;
	}
	if (verbose_ > 3 && (h2tohh_lorentz_.Pt()/h2tohh_lorentz_.E())>0.0000001) {
	  std::cout << " h2tohh mass "<< h2tohh_lorentz_.M() <<" pt " << h2tohh_lorentz_.Pt();
	  h2tohh_lorentz_.Print();
	}
	onshellW_Eta_ = onshellW_lorentz_.Eta();
	onshellW_Phi_ = onshellW_lorentz_.Phi();
	onshellW_Pt_ = onshellW_lorentz_.Pt();
	onshellW_E_ = onshellW_lorentz_.E();
	onshellW_Mass_ = onshellW_lorentz_.M();
	offshellW_Eta_ = offshellW_lorentz_.Eta();
	offshellW_Phi_ = offshellW_lorentz_.Phi();
	offshellW_Pt_ = offshellW_lorentz_.Pt();
	offshellW_E_ = offshellW_lorentz_.E();
	offshellW_Mass_ = offshellW_lorentz_.M();
	htoWW_Eta_ = htoWW_lorentz_.Eta();
	htoWW_Phi_ = htoWW_lorentz_.Phi();
	htoWW_Pt_ = htoWW_lorentz_.Pt();
	htoWW_E_ = htoWW_lorentz_.E();
	htoWW_Mass_ = htoWW_lorentz_.M();
	htoBB_jets_Eta_ = htoBB_lorentz_.Eta();
	htoBB_jets_Phi_ = htoBB_lorentz_.Phi();
	htoBB_jets_Pt_ = htoBB_lorentz_.Pt();
	htoBB_jets_E_ = htoBB_lorentz_.E();
	htoBB_jets_Mass_ = htoBB_lorentz_.M();
	h2tohh_Pt_ = h2tohh_lorentz_.Pt();
	h2tohh_E_ = h2tohh_lorentz_.E();
	h2tohh_Mass_ = h2tohh_lorentz_.M();
	heavyMassEstimatormet_Px_ = met_vec2_.Px();
	heavyMassEstimatormet_Py_ = met_vec2_.Py();
	heavyMassEstimatormet_E_ = met_vec2_.Mod();
	heavyMassEstimatormet_Phi_ = met_vec2_.Phi();

	if (weightfromonshellnupt_func_) weight1_ = weightfromonshellnupt(nu_onshellW_pt); 
	if (weightfromonshellnupt_hist_) weight1_ = weightfromhist(onshellnupt_hist_, nu_onshellW_pt); 
	if (weightfromonoffshellWmass_hist_) weight2_ = weightfromhist(onoffshellWmass_hist_, wmass_gen_, offshellW_lorentz_.M()); 
	if (weightfromonoffshellWmass_hist_) weight3_ = weightfromhist(onoffshellWmass_hist_, wmass_gen_, offshellW_lorentz_.M(), false);
	if (weightfrombjetrescalec1c2_hist_) weight4_ = weightfromhist(bjetrescalec2_hist_, rescalec2_);
	weight1_ = weight1_*weight_;
	weight2_ = weight2_*weight1_;
	weight3_ = weight1_*weight3_;
	weight4_ = weight4_*weight1_;
	if ((h2tohh_lorentz_.Pt()/h2tohh_lorentz_.E())>0.0000001){
	  h2tohh_Eta_ = h2tohh_lorentz_.Eta();
	  h2tohh_Phi_ = h2tohh_lorentz_.Phi();
	}else {//pt =0, strange case here
	  h2tohh_Eta_ = 1000000;
	  h2tohh_Phi_ = 0;
	}

	//printheavyMassEstimatorresult();
	//if ( writehmetree_ ) hmetree_->Fill();
	if (weight1_<=0.0) continue;
	heavyMassEstimator_h2Mass_->Fill(h2tohh_Mass_, weight_);
	heavyMassEstimator_h2Massweight1_->Fill(h2tohh_Mass_, weight1_);
	heavyMassEstimator_h2Massweight4_->Fill(h2tohh_Mass_, weight4_);
	if (heavyMassEstimatordebug_)  std::cout <<" h2tohh mass "<< h2tohh_Mass_ <<" weight " << weight_ <<" weight1 "<< weight1_ <<" nu_onshellW_pt "<< nu_onshellW_pt << std::endl;
        validrun =true;
    }//end controls loop,(0,1,2,3)
    //hmetree_->Fill();
  }//end of iteration
  if (heavyMassEstimatordebug_){
      std::cout <<"initial heavyMassEstimator input met  px "<<hmemet_vec2_.Px() << " py "<<hmemet_vec2_.Py() <<" pt "<< hmemet_vec2_.Mod() <<std::endl;
      std::cout <<"last iteration heavyMassEstimator input met  px "<<met_vec2_.Px() << " py "<<met_vec2_.Py() <<" pt "<< met_vec2_.Mod() <<std::endl;
      if (validrun) std::cout <<" true validrun "<< std::endl;
      else std::cout <<" false validrun "<< std::endl;
      //std::cout <<"last iteration bjets input M_h= "<< htoBB_lorentz_.M(); htoBB_lorentz_.Print();
      //std::cout <<"num of solutions " << heavyMassEstimator_h2Mass_->GetEntries() << std::endl;
      }

  //std::cout <<"gFile get name "<< gFile->GetName() <<" gFile get options " << gFile->GetOption() << std::endl;
  //file_->Close();
  return validrun;
}

//------------ method called to initialize a tree for heavyMassEstimator for this event ------------
void
heavyMassEstimator::initTree(TTree* hmetree){

  //std::cout <<" init tree "<< hmetree->GetTitle() << std::endl; 
  //initial branch value if necessary
  //
  //
  weight1_ = 1.0;
  weight2_ = 1.0;
  weight3_ = 1.0;
  weight4_ = 1.0;

  if (simulation_ && onshellMarker_ == 1){
    eta_nuoffshellW_true_ = nu2_lorentz_true_.Eta();
    phi_nuoffshellW_true_ = nu2_lorentz_true_.Phi();
    pt_nuoffshellW_true_ = nu2_lorentz_true_.Pt();
    px_nuoffshellW_true_ = nu2_lorentz_true_.Px();
    py_nuoffshellW_true_ = nu2_lorentz_true_.Py();
    eta_nuonshellW_true_ = nu1_lorentz_true_.Eta();
    phi_nuonshellW_true_ = nu1_lorentz_true_.Phi();
    pt_nuonshellW_true_ = nu1_lorentz_true_.Pt();
    px_nuonshellW_true_ = nu1_lorentz_true_.Px();
    py_nuonshellW_true_ = nu1_lorentz_true_.Py();
  }
  else if (simulation_ && onshellMarker_ == 2){
    eta_nuoffshellW_true_ = nu1_lorentz_true_.Eta();
    phi_nuoffshellW_true_ = nu1_lorentz_true_.Phi();
    pt_nuoffshellW_true_ = nu1_lorentz_true_.Pt();
    px_nuoffshellW_true_ = nu1_lorentz_true_.Px();
    py_nuoffshellW_true_ = nu1_lorentz_true_.Py();
    eta_nuonshellW_true_ = nu2_lorentz_true_.Eta();
    phi_nuonshellW_true_ = nu2_lorentz_true_.Phi();
    pt_nuonshellW_true_ = nu2_lorentz_true_.Pt();
    px_nuonshellW_true_ = nu2_lorentz_true_.Px();
    py_nuonshellW_true_ = nu2_lorentz_true_.Py();
  }

  if (simulation_){ 
    htoBB_Eta_ = htoBB_lorentz_true_.Eta();
    htoBB_Phi_ = htoBB_lorentz_true_.Phi();
    htoBB_Pt_ = htoBB_lorentz_true_.Pt();
    htoBB_E_ = htoBB_lorentz_true_.E();
    htoBB_Mass_ = htoBB_lorentz_true_.M();
    b1_Eta_ = b1_lorentz_.Eta();
    b1_Phi_ = b1_lorentz_.Phi();
    b1_Pt_ = b1_lorentz_.Pt();
    b1_Px_ = b1_lorentz_.Px();
    b1_Py_ = b1_lorentz_.Py();
    b1_E_ = b1_lorentz_.E();
    b1_Mass_ = b1_lorentz_.M();
    b2_Eta_ = b2_lorentz_.Eta();
    b2_Phi_ = b2_lorentz_.Phi();
    b2_Pt_ = b2_lorentz_.Pt();
    b2_Px_ = b2_lorentz_.Px();
    b2_Py_ = b2_lorentz_.Py();
    b2_E_ = b2_lorentz_.E();
    b2_Mass_ = b2_lorentz_.M();

    mass_offshellW_true_ = offshellW_lorentz_true_.M();
    mass_onshellW_true_ = onshellW_lorentz_true_.M();
    mass_htoWW_true_ = htoWW_lorentz_true_.M();
    pt_h2tohh_true_ = h2tohh_lorentz_true_.Pt();
    mass_h2tohh_true_ = h2tohh_lorentz_true_.M();
    ideal_met_Px_ = ideal_met_lorentz_.Px();
    ideal_met_Py_ = ideal_met_lorentz_.Py();
    ideal_met_E_ = ideal_met_lorentz_.Energy();
  }
  else {
    eta_nuoffshellW_true_ = -1;
    phi_nuoffshellW_true_ = -1;
    pt_nuoffshellW_true_ = -1;
    eta_nuonshellW_true_ = -1;
    phi_nuonshellW_true_ = -1;
    pt_nuonshellW_true_ = -1;

    mass_offshellW_true_ = -1;
    mass_onshellW_true_ = -1;
    mass_htoWW_true_ = -1;
    pt_h2tohh_true_ = -1;
    mass_h2tohh_true_ = -1;



  }
  b1jet_Eta_ = hme_b1jet_lorentz_.Eta();
  b1jet_Phi_ = hme_b1jet_lorentz_.Phi();
  b1jet_Pt_ = hme_b1jet_lorentz_.Pt();
  b1jet_Px_ = hme_b1jet_lorentz_.Px();
  b1jet_Py_ = hme_b1jet_lorentz_.Py();
  b1jet_Energy_ = hme_b1jet_lorentz_.Energy();
  b1jet_Mass_ = hme_b1jet_lorentz_.M();
  b2jet_Eta_ = hme_b2jet_lorentz_.Eta();
  b2jet_Phi_ = hme_b2jet_lorentz_.Phi();
  b2jet_Pt_ = hme_b2jet_lorentz_.Pt();
  b2jet_Px_ = hme_b2jet_lorentz_.Px();
  b2jet_Py_ = hme_b2jet_lorentz_.Py();
  b2jet_Energy_ = hme_b2jet_lorentz_.Energy();
  b2jet_Mass_ = hme_b2jet_lorentz_.M();

  if (simulation_){
    b1jet_dR_ = hme_b1jet_lorentz_.DeltaR(b1_lorentz_);
    b2jet_dR_ = hme_b2jet_lorentz_.DeltaR(b2_lorentz_);
    //b1jet_dR_ = deltaR(b1jet_Eta_, b1jet_Phi_, b1_lorentz_.Eta(), b1_lorentz_.Phi());
    //b2jet_dR_ = deltaR(b2jet_Eta_, b2jet_Phi_, b2_lorentz_.Eta(), b2_lorentz_.Phi());
    b1rescalefactor_true_ = b1_Pt_/b1jet_Pt_;
    b2rescalefactor_true_ = b2_Pt_/b2jet_Pt_;
    if (b1jet_Pt_>b2jet_Pt_){
	rescalec1_true_ = b1rescalefactor_true_;
	rescalec2_true_ = b2rescalefactor_true_;
    }
    else {
	rescalec1_true_ = b2rescalefactor_true_;
	rescalec2_true_ = b1rescalefactor_true_;
    }
  }
  met_ = hmemet_vec2_.Mod(); 
  met_px_ = hmemet_vec2_.Px();
  met_py_ = hmemet_vec2_.Py();
  met_phi_ = hmemet_vec2_.Phi();

  hmetree->Branch("ievent", &iev_);
  hmetree->Branch("eta_mean", &eta_mean_);
  hmetree->Branch("eta_rms", &eta_rms_);
  hmetree->Branch("eta_gen",&eta_gen_);
  hmetree->Branch("phi_gen",&phi_gen_);
  hmetree->Branch("wmass_gen",&wmass_gen_);
  hmetree->Branch("hmass_gen",&hmass_gen_);

  hmetree->Branch("mu_onshellW_eta", &mu_onshellW_Eta_);
  hmetree->Branch("mu_onshellW_phi", &mu_onshellW_Phi_);
  hmetree->Branch("mu_onshellW_pt", &mu_onshellW_Pt_);
  hmetree->Branch("mu_onshellW_E", &mu_onshellW_E_);
  hmetree->Branch("mu_offshellW_eta", &mu_offshellW_Eta_);
  hmetree->Branch("mu_offshellW_phi", &mu_offshellW_Phi_);
  hmetree->Branch("mu_offshellW_pt", &mu_offshellW_Pt_);
  hmetree->Branch("mu_offshellW_E", &mu_offshellW_E_);
  hmetree->Branch("nu_onshellW_eta", &nu_onshellW_Eta_);
  hmetree->Branch("nu_onshellW_phi", &nu_onshellW_Phi_);
  hmetree->Branch("nu_onshellW_pt", &nu_onshellW_Pt_);
  hmetree->Branch("nu_onshellW_E", &nu_onshellW_E_);
  hmetree->Branch("nu_offshellW_eta", &nu_offshellW_Eta_);
  hmetree->Branch("nu_offshellW_phi", &nu_offshellW_Phi_);
  hmetree->Branch("nu_offshellW_pt", &nu_offshellW_Pt_);
  hmetree->Branch("nu_offshellW_E", &nu_offshellW_E_);
  hmetree->Branch("onshellW_eta", &onshellW_Eta_);
  hmetree->Branch("onshellW_phi", &onshellW_Phi_);
  hmetree->Branch("onshellW_pt", &onshellW_Pt_);
  hmetree->Branch("onshellW_E", &onshellW_E_);
  hmetree->Branch("onshellW_Mass", &onshellW_Mass_);
  hmetree->Branch("offshellW_eta", &offshellW_Eta_);
  hmetree->Branch("offshellW_phi", &offshellW_Phi_);
  hmetree->Branch("offshellW_pt", &offshellW_Pt_);
  hmetree->Branch("offshellW_E", &offshellW_E_);
  hmetree->Branch("offshellW_Mass", &offshellW_Mass_);
  hmetree->Branch("htoWW_Eta", &htoWW_Eta_);
  hmetree->Branch("htoWW_Phi", &htoWW_Phi_);
  hmetree->Branch("htoWW_Pt", &htoWW_Pt_);
  hmetree->Branch("htoWW_E", &htoWW_E_);
  hmetree->Branch("htoWW_Mass", &htoWW_Mass_);

  hmetree->Branch("b1_Eta", &b1_Eta_);
  hmetree->Branch("b1_Phi", &b1_Phi_);
  hmetree->Branch("b1_Pt", &b1_Pt_);
  hmetree->Branch("b1_Px", &b1_Px_);
  hmetree->Branch("b1_Py", &b1_Py_);
  hmetree->Branch("b1_E", &b1_E_);
  hmetree->Branch("b1_Mass", &b1_Mass_);
  hmetree->Branch("b2_Eta", &b2_Eta_);
  hmetree->Branch("b2_Phi", &b2_Phi_);
  hmetree->Branch("b2_Pt", &b2_Pt_);
  hmetree->Branch("b2_Px", &b2_Px_);
  hmetree->Branch("b2_Py", &b2_Py_);
  hmetree->Branch("b2_E", &b2_E_);
  hmetree->Branch("b2_Mass", &b2_Mass_);

  hmetree->Branch("b1jet_Eta", &b1jet_Eta_);
  hmetree->Branch("b1jet_Phi", &b1jet_Phi_);
  hmetree->Branch("b1jet_Pt", &b1jet_Pt_);
  hmetree->Branch("b1jet_Px", &b1jet_Px_);
  hmetree->Branch("b1jet_Py", &b1jet_Py_);
  hmetree->Branch("b1jet_E", &b1jet_Energy_);
  hmetree->Branch("b1jet_Mass", &b1jet_Mass_);
  hmetree->Branch("b2jet_Eta", &b2jet_Eta_);
  hmetree->Branch("b2jet_Phi", &b2jet_Phi_);
  hmetree->Branch("b2jet_Pt", &b2jet_Pt_);
  hmetree->Branch("b2jet_Px", &b2jet_Px_);
  hmetree->Branch("b2jet_Py", &b2jet_Py_);
  hmetree->Branch("b2jet_E", &b2jet_Energy_);
  hmetree->Branch("b2jet_Mass", &b2jet_Mass_);
  hmetree->Branch("b2jet_dR", &b2jet_dR_);
  hmetree->Branch("b1jet_dR", &b1jet_dR_);

  hmetree->Branch("htoBB_Eta", &htoBB_Eta_);
  hmetree->Branch("htoBB_Phi", &htoBB_Phi_);
  hmetree->Branch("htoBB_Pt", &htoBB_Pt_);
  hmetree->Branch("htoBB_E", &htoBB_E_);
  hmetree->Branch("htoBB_Mass", &htoBB_Mass_);
  hmetree->Branch("b1rescalefactor",&b1rescalefactor_);
  hmetree->Branch("b2rescalefactor",&b2rescalefactor_);
  hmetree->Branch("b1rescalefactor_true_",&b1rescalefactor_true_);
  hmetree->Branch("b2rescalefactor_true_",&b2rescalefactor_true_);
  hmetree->Branch("rescalec1",&rescalec1_);
  hmetree->Branch("rescalec2",&rescalec2_);
  hmetree->Branch("rescalec1_true_",&rescalec1_true_);
  hmetree->Branch("rescalec2_true_",&rescalec2_true_);

  hmetree->Branch("htoBB_jets_Eta", &htoBB_jets_Eta_);
  hmetree->Branch("htoBB_jets_Phi", &htoBB_jets_Phi_);
  hmetree->Branch("htoBB_jets_Pt", &htoBB_jets_Pt_);
  hmetree->Branch("htoBB_jets_E", &htoBB_jets_E_);
  hmetree->Branch("htoBB_jets_Mass", &htoBB_jets_Mass_);

  hmetree->Branch("heavyMassEstimatormet_E",&heavyMassEstimatormet_E_);
  hmetree->Branch("heavyMassEstimatormet_Phi",&heavyMassEstimatormet_Phi_);
  hmetree->Branch("heavyMassEstimatormet_Px",&heavyMassEstimatormet_Px_);
  hmetree->Branch("heavyMassEstimatormet_Py",&heavyMassEstimatormet_Py_);
  hmetree->Branch("ideal_met_Px", &ideal_met_Px_);
  hmetree->Branch("ideal_met_Py", &ideal_met_Py_);
  hmetree->Branch("ideal_met_E", &ideal_met_E_);


  hmetree->Branch("h2tohh_Eta", &h2tohh_Eta_);
  hmetree->Branch("h2tohh_Phi", &h2tohh_Phi_);
  hmetree->Branch("h2tohh_Pt", &h2tohh_Pt_);
  hmetree->Branch("h2tohh_E", &h2tohh_E_);
  hmetree->Branch("h2tohh_Mass", &h2tohh_Mass_);


  hmetree->Branch("met_true_",&met_);
  hmetree->Branch("met_phi_true_",&met_phi_);
  hmetree->Branch("met_px_true_",&met_px_);
  hmetree->Branch("met_py_true_",&met_py_);

  hmetree->Branch("eta_nuoffshellW_true_", &eta_nuoffshellW_true_);
  hmetree->Branch("phi_nuoffshellW_true_", &phi_nuoffshellW_true_);
  hmetree->Branch("eta_nuonshellW_true_", &eta_nuonshellW_true_);
  hmetree->Branch("phi_nuonshellW_true_", &phi_nuonshellW_true_);
  hmetree->Branch("pt_nuoffshellW_true_", &pt_nuoffshellW_true_);
  hmetree->Branch("px_nuoffshellW_true_", &px_nuoffshellW_true_);
  hmetree->Branch("py_nuoffshellW_true_", &py_nuoffshellW_true_);
  hmetree->Branch("pt_nuonshellW_true_", &pt_nuonshellW_true_);
  hmetree->Branch("px_nuonshellW_true_", &px_nuonshellW_true_);
  hmetree->Branch("py_nuonshellW_true_", &py_nuonshellW_true_);
  hmetree->Branch("mass_offshellW_true_", &mass_offshellW_true_);
  hmetree->Branch("mass_onshellW_true_", &mass_onshellW_true_);
  hmetree->Branch("mass_h2_true_", &mass_h2tohh_true_);
  hmetree->Branch("pt_h2_true_", &pt_h2tohh_true_);
  hmetree->Branch("mass_htoWW_true_", &mass_htoWW_true_);
  hmetree->Branch("mass_h2_expect", &mass_h2_expect_);

  hmetree->Branch("weight", &weight_);
  hmetree->Branch("weight1", &weight1_);
  hmetree->Branch("weight2", &weight2_);
  hmetree->Branch("weight3", &weight3_);
  hmetree->Branch("weight4", &weight4_);
  hmetree->Branch("control", &control_);

  //also init tree


}


TLorentzVector 
heavyMassEstimator::calculateMET(){

  TLorentzVector METlorentz = TLorentzVector();
  //TVector2 met_pxpy(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py());
  //METlorentz.SetPxPyPzE(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(),0,met_pxpy.Mod());

  return METlorentz;
}


//------------ method called to assign muons lorenz vector --------------
void 
heavyMassEstimator::assignMuLorentzVec(int control){

  //  runheavyMassEstimator() control/2 == 0, namely control =0 here, we have correct muon lorentz Vector pair
  //
  //std::cout <<" beign muon assignment " << std::endl;
  if (simulation_){
    if (onshellMarker_ == 1 && control == 0){
	mu_onshellW_lorentz_ = hme_lep1_lorentz_;
	mu_offshellW_lorentz_ = hme_lep2_lorentz_; }
    else if (onshellMarker_ == 1 && control == 1){
	mu_onshellW_lorentz_ = hme_lep2_lorentz_;
	mu_offshellW_lorentz_ = hme_lep1_lorentz_;}
    else if (onshellMarker_ == 2 && control == 0){
	mu_onshellW_lorentz_ = hme_lep2_lorentz_;
	mu_offshellW_lorentz_ = hme_lep1_lorentz_;}
    else if (onshellMarker_ == 2 && control == 1){
	mu_onshellW_lorentz_ = hme_lep1_lorentz_;
	mu_offshellW_lorentz_ = hme_lep2_lorentz_;}
  }//simulation case
  else {
    if (control == 0){
	mu_onshellW_lorentz_ = hme_lep1_lorentz_;
	mu_offshellW_lorentz_ = hme_lep2_lorentz_;}
    else if (control == 1){
	mu_onshellW_lorentz_ = hme_lep2_lorentz_;
	mu_offshellW_lorentz_ = hme_lep1_lorentz_;}

  }//real case, assign them randomly
  //std::cout <<" end muon assignment " << std::endl;


}

// ------------ method called to generate a pair (eta,phi) for nuetrino1  ------------
EtaPhi 
heavyMassEstimator::generatenu1_etaphi(){

  float eta=0.0;
  float phi=0.0;

  float mean=0;
  float rms=1.403;
  eta = genEtaGuass(mean, rms);
  phi = genPhiFlat();

  return std::make_pair(eta, phi);
}

// ------------ method called to generate eta from Gauss distribution  ------------
float 
heavyMassEstimator::genEtaGuass(float mean, float rms){

  float eta = rnd_->Gaus(mean, rms);

  return eta;

}

// ------------ method called to generate phi from Flat distribution  ------------
float 
heavyMassEstimator::genPhiFlat(){

  float pi = 3.1415926;
  float phi = rnd_->Uniform(-pi, pi);

  return phi;
}

//------------ method to load TH1F object from ROOT file -----------------------------
const TH1F* loadHistogram1d(const TFile* file, const std::string& histName)
{
  const TH1F* hist = dynamic_cast<TH1F*>((const_cast<TFile*>(file))->Get(histName.data()));
  assert(hist);
  std::string histName_cloned = std::string(histName).append("_cloned");
  const TH1F* hist_cloned = (TH1F*)hist->Clone(histName_cloned.data());
  return hist_cloned;
}

//------------ method to load TH1F object from ROOT file -----------------------------
const TH2F* loadHistogram2d(const TFile* file, const std::string& histName)
{
  const TH2F* hist = dynamic_cast<TH2F*>((const_cast<TFile*>(file))->Get(histName.data()));
  assert(hist);
  std::string histName_cloned = std::string(histName).append("_cloned");
  const TH2F* hist_cloned = (TH2F*)hist->Clone(histName_cloned.data());
  return hist_cloned;
}

//------------ method called to readout TH1F onshellWmasspdf from root file -----------------------------
//
const TH1F*
heavyMassEstimator::readoutonshellWMassPDF(){

  const TH1F* onshellWmasspdf = loadHistogram1d(file_, "onshellWmasspdf");
  return onshellWmasspdf;

}

//------------ method called to readout TH1F offshellWmasspdf from root file -----------------------------
//
const TH1F*
heavyMassEstimator::readoutoffshellWMassPDF(){

  const TH1F* offshellWmasspdf = loadHistogram1d(file_, "offshellWmasspdf");
  return offshellWmasspdf;

}


//------------ method called to readout TH2F onoffshellWmasspdf from root file -----------------------------
//
const TH2F*
heavyMassEstimator::readoutonoffshellWMassPDF(){

  const TH2F* onoffshellWmasspdf = loadHistogram2d(file_, "onoffshellWmasspdf");
  return onoffshellWmasspdf;

}


//------------ method called to readout TH1F onshellnuptpdf from root file -----------------------------
//
const TH1F*
heavyMassEstimator::readoutonshellnuptPDF(){

  const TH1F* onshellnuptpdf = loadHistogram1d(file_, "onshellnuptpdf");
  return onshellnuptpdf;

}

//------------ method called to readout TH1F bjetrescalec1pdf from root file -----------------------------
//
const TH1F*
heavyMassEstimator::readoutbjetrescalec1PDF(){

  std::string histName = ( PUsample_ ) ? "recobjetrescalec1pdfPU40" : "bjetrescalec1dR4pdf"; 
  const TH1F* bjetrescalec1pdf = loadHistogram1d(file_, histName);
  return bjetrescalec1pdf;

}

//------------ method called to readout TH1F bjetrescalec2pdf from root file -----------------------------
//
const TH1F*
heavyMassEstimator::readoutbjetrescalec2PDF(){

  std::string histName = ( PUsample_ ) ? "recobjetrescalec2pdfPU40" : "bjetrescalec2dR4pdf";
  const TH1F* bjetrescalec2pdf = loadHistogram1d(file_, histName);
  return bjetrescalec2pdf;

}

//------------ method called to readout TH2F bjetrescalec1c2pdf from root file -----------------------------
//
const TH2F*
heavyMassEstimator::readoutbjetrescalec1c2PDF(){

  const TH2F* bjetrescalec1c2pdf = loadHistogram2d(file_, "bjetrescalec1c2pdf");
  return bjetrescalec1c2pdf;

}

//------------ method to describe onshellW mass Probability density function ------------------------------
//
float 
heavyMassEstimator::onshellWMassPDF(float mass){

  // float sigma = 1.75;
  // float mean = 80.1;
  float p0 =7.87161e-03;
  float p1 =1.69085;
  float p2 =603.474 ;
  float p = 0;
  p = exp(mass*p0+p1)+p2*exp(-0.5*((mass-80.1)/2.00)*((mass-80.1)/2.00));
  return p;
}

//------------ use random walk to generate random onshellW mass accroding to wmass pdf --------------
//
float
heavyMassEstimator::onshellWMassRandomWalk(float x0, float step, float random){

  float xmin = 50;
  float xmax = 90;
  float x1 = x0+step;
  while (x1 > xmax || x1 < xmin){
    if (x1 > xmax) x1 = x1-xmax+xmin;
    if (x1 < xmin) x1 = xmax-(xmin-x1);
  }
  //transition probability
  float w = onshellWMassPDF(x1)/onshellWMassPDF(x0);
  //std::cout <<" initial " <<x0 <<" step " << step << " x1 "<< x1 << " transition probability " << w << " random " << random << std::endl;
  if (w >= 1.00) return x1;
  if (w < 1.00 && random < w) return x1;
  else return x0;

}  


//------------ use random walk to generate random onshellW mass accroding to wmass pdf --------------
//
float
heavyMassEstimator::onshellWMassRandomWalk(float x0, float step, float random, const TH1F* hist){
  float xmin = 50;
  float xmax = 90;
  //periodic boundary codition
  while (x0 > xmax || x0 < xmin){
    if (x0 > xmax) x0 = x0-xmax+xmin;
    if (x0 < xmin) x0 = xmax-(xmin-x0);
  }

  float x1 = x0+step;
  while (x1 > xmax || x1 < xmin){
    if (x1 > xmax) x1 = x1-xmax+xmin;
    if (x1 < xmin) x1 = xmax-(xmin-x1);
  }
  //find
  int binx0_1,binx0_2;
  int binx1_1,binx1_2;
  double bincent0_1,bincont0_1;// center and content
  double bincent1_1,bincont1_1;
  binx0_1 = (const_cast<TH1F*>(hist))->FindBin(x0);
  binx1_1 = (const_cast<TH1F*>(hist))->FindBin(x1);
  if ((float)hist->GetBinCenter(binx0_1) < x0){
    binx0_2 = binx0_1+1;
  }
  else {
    binx0_2 = binx0_1;
    binx0_1 = binx0_1-1;
  }

  if ((float)hist->GetBinCenter(binx1_1) < x1){
    binx1_2 = binx1_1+1;
  }
  else {
    binx1_2 = binx1_1;
    binx1_1 = binx1_1-1;
  }
  bincent0_1 = hist->GetBinCenter(binx0_1);
  bincont0_1 = hist->GetBinContent(binx0_1);
  bincent1_1 = hist->GetBinCenter(binx1_1);
  bincont1_1 = hist->GetBinContent(binx1_1);
  double w0 = (x0-bincent0_1)*(bincont0_1-hist->GetBinContent(binx0_2))/(bincent0_1-hist->GetBinCenter(binx0_2))+bincont0_1;
  double w1 = (x1-bincent1_1)*(bincont1_1-hist->GetBinContent(binx1_2))/(bincent1_1-hist->GetBinCenter(binx1_2))+bincont1_1;
  //transition probability
  double w = w1/w0;

  //std::cout <<" initial " <<x0 <<" step " << step << " x1 "<< x1 << " transition probability " << w << " random " << random << std::endl;
  if (w >= 1.00) return x1;
  if (w < 1.00 && random < (float)w) return x1;
  else return x0;

}  


//---------- weight solution by a histogram --------------------------------------------------------
//
float
heavyMassEstimator::weightfromhist(const TH1F* hist, float x){
  //hist should be scaled

  float weight = 0.0;
  int bin1 = (const_cast<TH1F*>(hist))->FindBin(x);
  //first make sure that x is within range
  if (bin1 == 0 || bin1 == hist->GetNbinsX()+1) return weight=0;

  weight = (const_cast<TH1F*>(hist))->Interpolate(x);
  return weight;
}

//---------- weight solution by a 2d histogram --------------------------------------------------------
//
float
heavyMassEstimator::weightfromhist(const TH2F* hist, float x, float y, bool whole){
  //hist should be scaled

  float weight = 0.0;
  int bin1 = hist->GetXaxis()->FindBin(x);
  int bin2 = hist->GetYaxis()->FindBin(y);
  //first make sure that x is within range
  if (bin1 == 0 || bin1 == hist->GetNbinsX()+1) return weight=0;
  if (bin2 == 0 || bin2 == hist->GetNbinsY()+1) return weight=0;
  weight = hist->GetBinContent(bin1, bin2);
  if (whole){
    return weight;
  }else {
    if (hist->GetBinContent(bin1, bin2) < .1) return 0;
    float integral = hist->Integral(bin1,bin1, 0, hist->GetNbinsY()+1);
    return weight/integral;
  }

}


//---------- weight solution by nupt --------------------------------------------------------
//
float
heavyMassEstimator::weightfromonshellnupt(float nupt){

  float weight = 0.0;
  float max = 170;
  if (nupt<0 || nupt>125) return 0.0;

  weight = -16.925+12.4066*nupt-0.2884*std::pow(nupt,2)+0.00203*std::pow(nupt,3)+7.695e-7*std::pow(nupt,4)
    -7.2191e-8*std::pow(nupt,5)+2.499e-10*std::pow(nupt,6);
  if (weight < 0 && nupt<5) return 0.0;
  if (weight < 0 && verbose_ > 0) std::cout << " error! nupt " << nupt << " weight " << weight << std::endl;
  weight = weight/max;
  return weight;
}


//------------- method called to calculate pt of nuetrinos from on-shell W decay ------------
float 
heavyMassEstimator::nu1pt_onshellW(EtaPhi nu1_etaphi, const TLorentzVector& lep1lorentz, float wMass){
   
  float nu1_pt=0.0;
  //   TVector2 numu_phi = TVector2(nu_etaphi.first(),lep1lorentz.Eta());
  float deltaeta = nu1_etaphi.first - lep1lorentz.Eta();
  float deltaphi = nu1_etaphi.second - lep1lorentz.Phi();
  nu1_pt = wMass*wMass/(2*lep1lorentz.Pt()*(cosh(deltaeta)-cos(deltaphi)));
  return nu1_pt;

}

//------------ method called to check whether the solution in this case exist or not -------------
// not use now, may be helpful later 
bool
heavyMassEstimator::checkSolution(const TLorentzVector& jetslorentz,
    const TLorentzVector& lep1lorentz,
    const TLorentzVector& lep2lorentz,
    const TLorentzVector& nu1lorentz, int control, float hMass){

  TLorentzVector tmplorentz(lep1lorentz.Px()+lep2lorentz.Px()+nu1lorentz.Px(),
	lep1lorentz.Py()+lep2lorentz.Py()+nu1lorentz.Py(),
	lep1lorentz.Pz()+lep2lorentz.Pz()+nu1lorentz.Pz(),
	lep1lorentz.Energy()+lep2lorentz.Energy()+nu1lorentz.Energy());

  float nu_tmp_px;
  float nu_tmp_py;
  float nu_tmp_pt;

  nu_tmp_px = -jetslorentz.Px()-lep1lorentz.Px()-lep2lorentz.Px()-nu1lorentz.Px();
  nu_tmp_py = -jetslorentz.Py()-lep1lorentz.Py()-lep2lorentz.Py()-nu1lorentz.Py();
  TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);

  nu_tmp_pt = nu_pxpy.Mod();

  float chdeltaeta;//cosh(nu2_eta-tmp2lorenz_eta)
  TLorentzVector tmp2lorentz(sqrt(pow(tmplorentz.Pt(),2)+pow(tmplorentz.M(),2)),0,tmplorentz.Pz(),tmplorentz.Energy());// construct massless lorentzvector with same pz and E as tmplorentzvector

  chdeltaeta = (pow(hMass,2)+pow(jetslorentz.Pt(),2)-pow(tmplorentz.M(),2)-pow(tmplorentz.Pt(),2)-pow(nu_tmp_pt,2))/(2*tmp2lorentz.Pt()*nu_tmp_pt);

  // place the cuts we may need 
  //
  //at present if (|chdeltaeta|>1) return true; 
  return (fabs(chdeltaeta)>1);
}



//------------- method called to calculate lorentzvector of second nuetrinos, which is from offshell W -----------
// return true if we can get nu_offshellW_lorentz_
bool 
heavyMassEstimator::nulorentz_offshellW(const TLorentzVector& jetslorentz, 
    const TLorentzVector& lep1lorentz, 
    const TLorentzVector& lep2lorentz, 
    TLorentzVector& nu1lorentz, 
    TLorentzVector& nu2lorentz, int control, float hMass){

  TLorentzVector tmplorentz(lep1lorentz.Px()+lep2lorentz.Px()+nu1lorentz.Px(),
	lep1lorentz.Py()+lep2lorentz.Py()+nu1lorentz.Py(),
	lep1lorentz.Pz()+lep2lorentz.Pz()+nu1lorentz.Pz(),
	lep1lorentz.Energy()+lep2lorentz.Energy()+nu1lorentz.Energy());
  float nu_tmp_px;
  float nu_tmp_py;
  float nu_tmp_pt;

  nu_tmp_px = -jetslorentz.Px()-lep1lorentz.Px()-lep2lorentz.Px()-nu1lorentz.Px();
  nu_tmp_py = -jetslorentz.Py()-lep1lorentz.Py()-lep2lorentz.Py()-nu1lorentz.Py();
  TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);

  nu_tmp_pt = nu_pxpy.Mod();

  float chdeltaeta;//cosh(nu_offshellW_eta-tmp2lorentz_eta)
  TLorentzVector tmp2lorentz(sqrt(pow(tmplorentz.Pt(),2)+pow(tmplorentz.M(),2)),0,tmplorentz.Pz(),tmplorentz.Energy());//fake one massless lorentzvector with same pz and E

  chdeltaeta = (pow(hMass,2)+pow(jetslorentz.Pt(),2)-pow(tmplorentz.M(),2)-pow(tmplorentz.Pt(),2)-pow(nu_tmp_pt,2))/(2*tmp2lorentz.Pt()*nu_tmp_pt);

  if (verbose_ >0 ){

    std::cout << "From jetLorentz nu2 px: " << nu_tmp_px << " py: "<< nu_tmp_py << " chdeltaeta: " << chdeltaeta << std::endl;
    float chdeltaeta_tmp = (pow(hMass,2)+2*(hmemet_vec2_.Px()*tmp2lorentz.Px()+hmemet_vec2_.Py()*tmp2lorentz.Py())-pow(nu_tmp_pt,2))/(2*tmp2lorentz.Pt()*nu_tmp_pt);
    std::cout << "From hmemet nu2 px: "<< hmemet_vec2_.Px()-nu1lorentz.Px() <<" py: "<< hmemet_vec2_.Py()-nu1lorentz.Py()
	<<" chdeltaeta: " << chdeltaeta_tmp << std::endl; 
  }
  if (chdeltaeta < 1.0) {
    nu2lorentz.SetPtEtaPhiM(0, 0, 0, 0);
    return false;
  }
  float nu_tmp_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi());
  float deltaeta = acosh(chdeltaeta);
  float nu_tmp_eta = (control_ == 1) ? (tmp2lorentz.Eta()-deltaeta) : (tmp2lorentz.Eta()+deltaeta);//control_ = j%2 
  // should check whether deltaeta > 1
  // std::cout <<"control "<< control_ <<" nu_tmp_px " << nu_tmp_px << "  nu_tmp_py " << nu_tmp_py << " nu_tmp_pt " << nu_tmp_pt 
  //         << " cosh(deltaeta2) " << chdeltaeta << " nu_tmp_eta " << nu_tmp_eta << " nu_tmp_phi " << nu_tmp_phi << std::endl; 
  if (fabs(nu_tmp_eta) > 7) {
    nu2lorentz.SetPtEtaPhiM(0, 0, 0, 0);
    return false;  //from simulation, |nu_offshellW_Eta|<6
  }
  nu2lorentz.SetPtEtaPhiM(nu_tmp_pt, nu_tmp_eta, nu_tmp_phi, 0);
  TLorentzVector htoww_tmp = tmplorentz + nu2lorentz;
  if (abs(htoww_tmp.M()-hMass) >2 && verbose_ > 0){
    std::cout <<" set Higgs Mass" << hMass << " heavyMassEstimator higgs mass" << htoww_tmp.M() << std::endl;
    htoww_tmp.Print();
    //verbose_ = 1;
  }
  if (verbose_ > 0){
    std::cout << "tmplorentz mass " << tmplorentz.M(); tmplorentz.Print();
    std::cout << "tmp2lorentz mass " << tmp2lorentz.M(); tmp2lorentz.Print();
    std::cout << " jets lorentz"; jetslorentz.Print(); 
    std::cout << " lep1 lorentz "; lep1lorentz.Print();
    std::cout << " lep2 lorentz "; lep2lorentz.Print();
    std::cout << " nu1 lorentz "; nu1lorentz.Print();
    std::cout << " tmp lorentz "; tmplorentz.Print();
    std::cout << " nu2 lorentz "; nu2lorentz.Print();
  }
  // std::cout << " nu_offshellW lorentz "; nu2lorentz.Print();

  return true; 
}



//------------- method called to calculate lorentzvector of second nuetrinos, which is from offshell W -----------
// return true if we can get nu_offshellW_lorentz_
bool 
heavyMassEstimator::nulorentz_offshellW(const TVector2& met, 
    const TLorentzVector& lep1lorentz, 
    const TLorentzVector& lep2lorentz, 
    TLorentzVector& nu1lorentz, 
    TLorentzVector& nu2lorentz, int control, float hMass){

  TLorentzVector tmplorentz(lep1lorentz.Px()+lep2lorentz.Px()+nu1lorentz.Px(),
	lep1lorentz.Py()+lep2lorentz.Py()+nu1lorentz.Py(),
	lep1lorentz.Pz()+lep2lorentz.Pz()+nu1lorentz.Pz(),
	lep1lorentz.Energy()+lep2lorentz.Energy()+nu1lorentz.Energy());
  float nu_tmp_px;
  float nu_tmp_py;
  float nu_tmp_pt;

  nu_tmp_px = met.Px()-nu1lorentz.Px();
  nu_tmp_py = met.Py()-nu1lorentz.Py();
  TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);

  nu_tmp_pt = nu_pxpy.Mod();

  float chdeltaeta;//cosh(nu_offshellW_eta-tmp2lorentz_eta)
  TLorentzVector tmp2lorentz(sqrt(pow(tmplorentz.Pt(),2)+pow(tmplorentz.M(),2)),0,tmplorentz.Pz(),tmplorentz.Energy());//fake one massless lorentzvector with same pz and E

  chdeltaeta = (pow(hMass,2)+2*(nu_pxpy.Px()*tmplorentz.Px()+nu_pxpy.Py()*tmplorentz.Py())-pow(tmplorentz.M(),2))/(2*tmp2lorentz.Pt()*nu_tmp_pt);
  if (verbose_ >0 ){
    std::cout << "nu2 px: " << nu_tmp_px << " py: "<< nu_tmp_py << std::endl;
    std::cout << "chdeltaeta " << chdeltaeta << std::endl;
    std::cout << "tmp2lorentz "; tmp2lorentz.Print();
  }
  if (chdeltaeta < 1.0) {
    nu2lorentz.SetPtEtaPhiM(0, 0, 0, 0);
    return false;
  }
  float nu_tmp_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi());
  float deltaeta = acosh(chdeltaeta);
  float nu_tmp_eta = (control_ == 1) ? (tmp2lorentz.Eta()-deltaeta) : (tmp2lorentz.Eta()+deltaeta);//control = j%2 
  // should check whether deltaeta > 1
  // std::cout <<"control "<< control_ <<" nu_tmp_px " << nu_tmp_px << "  nu_tmp_py " << nu_tmp_py << " nu_tmp_pt " << nu_tmp_pt 
  //         << " cosh(deltaeta2) " << chdeltaeta << " nu_tmp_eta " << nu_tmp_eta << " nu_tmp_phi " << nu_tmp_phi << std::endl; 
  if (fabs(nu_tmp_eta) > 7) {
    nu2lorentz.SetPtEtaPhiM(0, 0, 0, 0);
    return false;  //from simulation, |nu_offshellW_Eta|<6
  }
  nu2lorentz.SetPtEtaPhiM(nu_tmp_pt, nu_tmp_eta, nu_tmp_phi, 0);
  TLorentzVector htoww_tmp = tmplorentz + nu2lorentz;
  if (abs(htoww_tmp.M()-hMass) >2 && verbose_ > 0){
    std::cout <<" set Higgs Mass" << hMass << " heavyMassEstimator higgs mass" << htoww_tmp.M() << std::endl;
    htoww_tmp.Print();
    //verbose_ = 1;
  }
  if (verbose_ > 0){
    std::cout << "tmplorentz mass " << tmplorentz.M(); tmplorentz.Print();
    std::cout << "tmp2lorentz mass " << tmp2lorentz.M(); tmp2lorentz.Print();
    std::cout << " met Tvector2 "; met.Print(); 
    std::cout << " lep1 lorentz "; lep1lorentz.Print();
    std::cout << " lep2 lorentz "; lep2lorentz.Print();
    std::cout << " nu1 lorentz "; nu1lorentz.Print();
    std::cout << " tmp lorentz "; tmplorentz.Print();
    std::cout << " nu2 lorentz "; nu2lorentz.Print();
  }
  // std::cout << " nu_offshellW lorentz "; nu2lorentz.Print();

  return true; 
}

//--------------------------- bjets correction, based on c1, calculate c2 here  ----------------------------------------------------------
//use rescalec1, rescalec2 to correct bjets
bool 
heavyMassEstimator::bjetsCorrection(){
  //c1rescale taken from pdf
  TLorentzVector b1lorentz;
  TLorentzVector b2lorentz;
  if (hme_b1jet_lorentz_.Pt()> hme_b2jet_lorentz_.Pt()){
    b1lorentz = hme_b1jet_lorentz_;
    b2lorentz = hme_b2jet_lorentz_;
  }
  else {
    if(verbose_ > 0) {
        std::cout <<"wired b1jet is not jet with larger pt "<< std::endl;
    }
    b1lorentz = hme_b2jet_lorentz_;
    b2lorentz = hme_b1jet_lorentz_;
  }
  //x1*c2*c2+x2*c2+x3=0, slove for c2
  float x1 = b2lorentz.M2();
  float x2 = 2*rescalec1_*(b1lorentz*b2lorentz);
  float x3 = rescalec1_*rescalec1_*b1lorentz.M2()-125*125;
  if (x2<0 && verbose_ > 0) std::cerr <<"error bjets lorentzvector dot productor less than 0 " << std::endl;
  /*std::cout <<" b2lorentz mass "<<  b2lorentz.M2();  b2lorentz.Print();
  std::cout <<" rescale1  "<< rescalec1_ <<" x1 "<< x1 <<" x2 " << x2 << " x3 "<< x3 << std::endl;
  std::cout <<"c2 solution1 " << (-x2+std::sqrt(x2*x2-4*x1*x3))/(2*x1)<<" solution2 "<<(-x2-std::sqrt(x2*x2-4*x1*x3))/(2*x1)<<std::endl;
  std::cout <<"b1jet pt "<< hme_b1jet_lorentz_.Pt() <<" b2jet pt " << hme_b2jet_lorentz_.Pt() << std::endl;
   */
  if ((x2*x2-4*x1*x3) <0 or x1==0){
	return false;
	std::cout <<" error ! there is no soluations for bjetsCorrection "<< std::endl;
  }
  rescalec2_ = (-x2+std::sqrt(x2*x2-4*x1*x3))/(2*x1);
  /*if (rescalec2_<=0){
	std::cout <<" b1jet mass "<< b1lorentz.M(); b1lorentz.Print(); 
	std::cout <<" b2jet mass "<< b2lorentz.M(); b2lorentz.Print(); 
	std::cout <<" x1 "<< x1 <<" x2 "<< x2 <<" x3 "<< x3 << std::endl;
  	std::cout <<"c2 solution1 " << (-x2+std::sqrt(x2*x2-4*x1*x3))/(2*x1)<<" solution2 "<<(-x2-std::sqrt(x2*x2-4*x1*x3))/(2*x1)<<std::endl;
	}*/
  if (hme_b1jet_lorentz_.Pt()> hme_b2jet_lorentz_.Pt()){
    b1rescalefactor_ = rescalec1_;
    b2rescalefactor_ = rescalec2_;
  }else{
    if(verbose_ > 0) {
        std::cout <<"wired b1jet is not jet with larger pt "<< std::endl;
    }
    b2rescalefactor_ = rescalec1_;
    b1rescalefactor_ = rescalec2_;
  }
  //TLorentzVector htobb_corr = b1rescalefactor_*b1lorentz+ b2rescalefactor_*b2lorentz;
  //finally b1rescalefactor_->b1jet b2rescalefactor_->b2jet;
  //std::cout <<" htoBB m_h after correction "<< htobb_corr.M(); htobb_corr.Print();
  //std::cout <<" htoBB hme m_h"<< htoBB_lorentz_.M(); htoBB_lorentz_.Print();
  return true;

}

//--------------------------- MET correction  ----------------------------------------------------------
void 
heavyMassEstimator::metCorrection(){
  //std::cout <<" b1rescalefactor " << b1rescalefactor << " b2rescalefactor " << b2rescalefactor << std::endl;
  float metpx_correction = hmemet_vec2_.Px()-(b1rescalefactor_-1)*hme_b1jet_lorentz_.Px()-(b2rescalefactor_-1)*hme_b2jet_lorentz_.Px();
  float metpy_correction = hmemet_vec2_.Py()-(b1rescalefactor_-1)*hme_b1jet_lorentz_.Py()-(b2rescalefactor_-1)*hme_b2jet_lorentz_.Py();
  met_vec2_ = TVector2(metpx_correction, metpy_correction);

  //std::cout <<" metCorrection metpx_correction " <<  metpx_correction <<" metpy_correction " << metpy_correction << std::endl;
}

//--------------------------- retrun heavyMassEstimator result ----------------------------------------------------------
const TH1F&
heavyMassEstimator::getheavyMassEstimatorh2(){

  //std::cout <<" RMS "<< heavyMassEstimator_h2Mass_->GetRMS() << " entries " << heavyMassEstimator_h2Mass_->GetEntries() << std::endl;
  return *heavyMassEstimator_h2Mass_;

}

//--------------------------- retrun heavyMassEstimator result ----------------------------------------------------------
const TH1F&
heavyMassEstimator::getheavyMassEstimatorh2weight1(){

  //std::cout <<" RMS "<< heavyMassEstimator_h2Mass_->GetRMS() << " entries " << heavyMassEstimator_h2Mass_->GetEntries() << std::endl;
  return *heavyMassEstimator_h2Massweight1_;

}

//--------------------------- retrun heavyMassEstimator result ----------------------------------------------------------
const TH1F&
heavyMassEstimator::getheavyMassEstimatorh2weight4(){

  //std::cout <<" RMS "<< heavyMassEstimator_h2Mass_->GetRMS() << " entries " << heavyMassEstimator_h2Mass_->GetEntries() << std::endl;
  return *heavyMassEstimator_h2Massweight4_;

}

//--------------------------- retrun heavyMassEstimator result ----------------------------------------------------------
const TTree*
heavyMassEstimator::getheavyMassEstimatorTree(){

  return hmetree_;

}



//----------------------------- print lorentz vectors from analyzer --------------------------------------------
//
void 
heavyMassEstimator::printTrueLorentz(){

  std::cout <<"  print out lorentz vector pass to heavyMassEstimator " << std::endl;
  if (simulation_) std::cout <<" onshell channel is " << onshellMarker_ << std::endl;
  std::cout <<" lep1 " ; hme_lep1_lorentz_.Print();
  std::cout <<" lep2 " ; hme_lep2_lorentz_.Print();
  std::cout <<"bjets,  M_h= " << htoBB_lorentz_.M(); htoBB_lorentz_.Print();
  std::cout <<"met px " << hmemet_vec2_.Px() <<" py" <<hmemet_vec2_.Py() << std::endl;
  if (simulation_) {
    std::cout <<"following is pure gen level infromation " << std::endl;
    std::cout <<" nu1 px "<<nu1_lorentz_true_.Px() << " py " <<nu1_lorentz_true_.Py() << " pt "<< nu1_lorentz_true_.Pt() 
	<< " eta "<<nu1_lorentz_true_.Eta() << " phi "<< nu1_lorentz_true_.Phi() << std::endl;
    std::cout <<" nu2 px "<<nu2_lorentz_true_.Px() << " py " <<nu2_lorentz_true_.Py() << " pt "<< nu2_lorentz_true_.Pt() 
	<< " eta "<<nu2_lorentz_true_.Eta() << " phi "<< nu2_lorentz_true_.Phi() << std::endl;
    std::cout <<" onshellW mass "<< onshellW_lorentz_true_.M(); onshellW_lorentz_true_.Print();  
    std::cout <<"offshellW mass " <<offshellW_lorentz_true_.M(); offshellW_lorentz_true_.Print();  
    std::cout <<" htoWW mass "<< htoWW_lorentz_true_.M(); htoWW_lorentz_true_.Print();
    std::cout <<" htoBB mass "<< htoBB_lorentz_true_.M(); htoBB_lorentz_true_.Print();
    std::cout <<" h2tohh, pz " <<h2tohh_lorentz_true_.Pz() << " Mass " << h2tohh_lorentz_true_.M() << std::endl;
  }
}


//----------------------------- print lorentz vectors from heavyMassEstimator results --------------------------------------------
void 
heavyMassEstimator::printheavyMassEstimatorresult(){

  std::cout <<" print out results from heavyMassEstimator, namely survival soltions " << std::endl;
  std::cout <<" onshell mu ";mu_onshellW_lorentz_.Print();
  std::cout <<" onshell nu ";nu_onshellW_lorentz_.Print();
  std::cout <<"offshell mu ";mu_offshellW_lorentz_.Print();
  std::cout <<"offshell nu ";nu_offshellW_lorentz_.Print();
  std::cout <<" onshell W  "; onshellW_lorentz_.Print();
  std::cout <<"offshell W  "; offshellW_lorentz_.Print();
  std::cout <<" htoBB      "; htoBB_lorentz_.Print();
  std::cout <<" h2tohh  pz"<< h2tohh_lorentz_.Pz() << " mass "<< h2tohh_lorentz_.M() << std::endl; 

}

#endif
