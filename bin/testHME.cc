using namespace std;
#include "../interface/heavyMassEstimator.h"

#include <iostream>
//ROOT

#include "TLorentzVector.h"

struct event {
    
   TLorentzVector lep1_p4;
   TLorentzVector lep2_p4;
   TLorentzVector b1jet_p4;
   TLorentzVector b2jet_p4;
   TLorentzVector met_p4;
   TLorentzVector totjets_p4;
   void printparticle(TLorentzVector tmp, string message){ std::cout << message <<" ( "<< tmp.Px()<<", "<< tmp.Py() <<", "<< tmp.Pz() <<", "<< tmp.Energy() <<") "<< std::endl;}
   void print(){
       printparticle(lep1_p4, "lepton1");
       printparticle(lep2_p4, "lepton2");
       printparticle(b1jet_p4, "b1jet");
       printparticle(b2jet_p4, "b2jet");
       printparticle(met_p4, "MET");
   }

};

int main(int argc, char *argv[])
{                 

    //event from h2tohh signal, M=400 
    event evlist[10];

    //first event: how HME input variables are constructed
    TLorentzVector nu1_p4, nu2_p4, nu_p4;
    nu1_p4.SetPxPyPzE(-6.67451, -55.3138, 47.3741, 73.1332);
    nu2_p4.SetPxPyPzE(-15.8995, -20.3445, 68.0647, 72.7976);
    nu_p4 = nu1_p4+nu2_p4;
    evlist[0].lep1_p4.SetPxPyPzE(-41.852, -16.4778, 1.89502, 45.019);//first lep
    evlist[0].lep2_p4.SetPxPyPzE(-1.81875, -18.7215, -6.85567, 20.0204);//second lep
    evlist[0].b1jet_p4.SetPxPyPzE(-62.1751, 82.7153, 190.54, 216.876); //first bjet
    evlist[0].b2jet_p4.SetPxPyPzE(-34.9725, -0.875287, 325.264, 327.173);//second bjet
    evlist[0].totjets_p4.SetPxPyPzE(-21.004410,17.735535,429.205327,557.353233);//sum of all jets with pt>20, |eta|<2.5, but usully not used in HME
    evlist[0].met_p4.SetPxPyPzE(nu_p4.Px(), nu_p4.Py(), 0, nu_p4.Pt());//met

    //real events from radion signal with M=400, narrow width
    evlist[1].lep1_p4.SetPxPyPzE(8.003993,47.391953,-6.067338,48.444656);
    evlist[1].lep2_p4.SetPxPyPzE(-10.217114,22.575941,-14.320704,28.620905);
    evlist[1].b1jet_p4.SetPxPyPzE(-50.568604,-73.598099,-107.297661,139.793747);
    evlist[1].b2jet_p4.SetPxPyPzE(-54.522942,-41.693054,24.696209,73.266808);
    evlist[1].totjets_p4.SetPxPyPzE(-72.348271,-37.953197,-143.769315,349.266855);
    evlist[1].met_p4.SetXYZM(57.0932, 58.3219, 0.0, 0.0);

    evlist[2].lep1_p4.SetPxPyPzE(-90.618065,-18.197742,-106.434578,140.964966);
    evlist[2].lep2_p4.SetPxPyPzE(-48.817795,1.260729,-61.836735,78.794411);
    evlist[2].b1jet_p4.SetPxPyPzE(109.058380,-43.815208,-134.934738,179.951065);
    evlist[2].b2jet_p4.SetPxPyPzE(55.421257,-17.291174,28.242508,65.646118);
    evlist[2].totjets_p4.SetPxPyPzE(25.986875,-56.534228,-241.287724,613.837237);
    evlist[2].met_p4.SetXYZM(-51.647, 56.8002, 0.0,0.0);
    //evlist[2].met_p4.SetXYZM(0.0,0.0);

    evlist[3].lep1_p4.SetPxPyPzE(-20.479839,61.715431,-7.073586,65.408447);
    evlist[3].lep2_p4.SetPxPyPzE(-15.327445,11.297248,-13.903660,23.577135);
    evlist[3].b1jet_p4.SetPxPyPzE(17.658623,-71.915688,11.239569,76.712524);
    evlist[3].b2jet_p4.SetPxPyPzE(55.094479,19.823299,-1.950312,58.932232);
    evlist[3].totjets_p4.SetPxPyPzE(31.018055,27.664730,-15.734976,235.484109);
    evlist[3].met_p4.SetXYZM(-61.773, 67.0091, 0.0,0.0);

 
    evlist[4].lep1_p4.SetPxPyPzE(-40.740326,-34.756702,41.160843,67.542793);
    evlist[4].lep2_p4.SetPxPyPzE(-46.830639,-14.090547,22.967823,54.029476);
    evlist[4].b1jet_p4.SetPxPyPzE(67.761467,56.386566,-12.797042,90.260597);
    evlist[4].b2jet_p4.SetPxPyPzE(11.292420,-27.073809,-39.151920,49.180706);
    evlist[4].totjets_p4.SetPxPyPzE(-21.354155,-26.302453,20.757114,278.642616);
    evlist[4].met_p4.SetXYZM(27.0983, 55.9654,0.0,0.0);

    evlist[5].lep1_p4.SetPxPyPzE(-25.073132,-48.923740,62.819412,83.477448);
    evlist[5].lep2_p4.SetPxPyPzE(-22.622330,6.422928,-27.284948,36.020874);
    evlist[5].b1jet_p4.SetPxPyPzE(72.667816,43.315197,-23.022663,88.966263);
    evlist[5].b2jet_p4.SetPxPyPzE(23.668390,-19.554901,-65.903366,72.972054);
    evlist[5].totjets_p4.SetPxPyPzE(1.658958,59.296977,386.106201,909.979552);
    evlist[5].met_p4.SetXYZM(1.14413, -53.0413,0.0,0.0);

    evlist[6].lep1_p4.SetPxPyPzE(-53.952442,133.644470,5.644073,144.234451);
    evlist[6].lep2_p4.SetPxPyPzE(-9.095098,24.572216,-6.853642,27.083172);
    evlist[6].b1jet_p4.SetPxPyPzE(-48.201416,-31.002796,-30.228876,65.041664);
    evlist[6].b2jet_p4.SetPxPyPzE(29.605610,-33.769974,-96.033966,106.385162);
    evlist[6].totjets_p4.SetPxPyPzE(-86.736130,108.754678,-127.794368,359.246207);
    evlist[6].met_p4.SetXYZM( -73.6003, 5.65941, 0.0,0.0);

 
    evlist[7].lep1_p4.SetPxPyPzE(-22.381014,50.897514,149.366302,159.379333);
    evlist[7].lep2_p4.SetPxPyPzE(-26.758770,27.116207,24.360731,45.219208);
    evlist[7].b1jet_p4.SetPxPyPzE(30.404705,-97.673676,206.070602,230.575378);
    evlist[7].b2jet_p4.SetPxPyPzE(65.745499,-27.916668,134.110291,152.165344);
    evlist[7].totjets_p4.SetPxPyPzE(-36.051720,5.614648,526.256547,934.900260);
    evlist[7].met_p4.SetXYZM(11.9203, 9.41891, 0.0,0.0);

 
    evlist[8].lep1_p4.SetPxPyPzE(-43.754608,109.016647,-136.843246,180.347397);
    evlist[8].lep2_p4.SetPxPyPzE(-11.171001,43.853645,-38.688629,59.537842);
    evlist[8].b1jet_p4.SetPxPyPzE(73.284187,-65.016823,-47.734806,109.983322);
    evlist[8].b2jet_p4.SetPxPyPzE(18.789080,-66.921143,-220.692169,231.570847);
    evlist[8].totjets_p4.SetPxPyPzE(33.074733,33.899972,-458.715957,601.624808);
    evlist[8].met_p4.SetXYZM(-27.3289, -12.2805,0.0,0.0);

 
    evlist[9].lep1_p4.SetPxPyPzE(-149.203964,44.106693,239.248840,285.389618);
    evlist[9].lep2_p4.SetPxPyPzE(-13.539427,-16.102177,56.490829,60.281185);
    evlist[9].b1jet_p4.SetPxPyPzE(31.546118,-54.460289,81.689758,103.500099);
    evlist[9].b2jet_p4.SetPxPyPzE(55.698837,6.437282,-5.329213,57.301113);
    evlist[9].totjets_p4.SetPxPyPzE(-31.302394,11.471296,459.979051,704.240572);
    //evlist[9].met_p4.SetPxPyPzE();
    evlist[9].met_p4.SetXYZM(12.5409, 1.68792,0.0,0.0);

 
    //evlist[0].lep1_p4.SetPxPyPzE();
    //evlist[0].lep2_p4.SetPxPyPzE();
    //evlist[0].b1jet_p4.SetPxPyPzE();
    //evlist[0].b2jet_p4.SetPxPyPzE();
    //evlist[0].totjets_p4.SetPxPyPzE();
    //evlist[0].met_p4.SetPxPyPzE();
 
 
    int nevent = 10;
    //HME configuration
    bool PUSample_ = true;//whether event is from PU sample or not.
    bool weightfromonshellnupt_func_ = false;
    bool weightfromonshellnupt_hist_ = true;
    bool weightfromonoffshellWmass_hist_ = true;
    string RefPDFfile_ = "../data/REFPDFPU40.root";//the root file contains histogram for weighting 
    
    int iterations_ = 10000;
    bool useMET_ = true;//use MET or totjets_p4 to estimate kinematic sum of two nuetrino
    int bjetrescaleAlgo_ = 2;//jet correction
    int metcorrection_ = 5;//met correction

    float h2tohh_mass = 400.0;//signal benmark M=400, narrow width. in other word, HME output of all above events should be close to 400.0
 
    for (int ievent=0; ievent <nevent; ievent++){
	evlist[ievent].print();
	heavyMassEstimator *thishme = new heavyMassEstimator(&evlist[ievent].lep1_p4, &evlist[ievent].lep2_p4, &evlist[ievent].b1jet_p4, &evlist[ievent].b2jet_p4, &evlist[ievent].totjets_p4, &evlist[ievent].met_p4, 
	    PUSample_, ievent, weightfromonshellnupt_func_, weightfromonshellnupt_hist_, weightfromonoffshellWmass_hist_,
	    iterations_, RefPDFfile_, useMET_, bjetrescaleAlgo_, metcorrection_);
	bool runheavyMassEstimatorok = thishme->runheavyMassEstimator();
	if (runheavyMassEstimatorok) {
	  //heavyMassEstimatortree =  (thishme->getheavyMassEstimatorTree())->CloneTree();
	    std::stringstream ss;
	  ss <<"heavyMassEstimator_h2mass_event"<< ievent;
	  const std::string histname(ss.str());
	  std::stringstream ss1;
	  ss1 <<"heavyMassEstimator_h2massweight1_event"<< ievent;
	  const std::string histname1(ss1.str());
	  std::stringstream ss4;
	  ss4 <<"heavyMassEstimator_h2massweight4_event"<< ievent;
	  const std::string histname4(ss4.str());
	  TH1F* heavyMassEstimator_h2mass =(TH1F*)(thishme->getheavyMassEstimatorh2()).Clone(histname.c_str());
	  TH1F* heavyMassEstimator_h2mass_weight1 =(TH1F*)(thishme->getheavyMassEstimatorh2weight1()).Clone(histname1.c_str());
	  TH1F* heavyMassEstimator_h2mass_weight4 =(TH1F*)(thishme->getheavyMassEstimatorh2weight4()).Clone(histname4.c_str());
	  //std::cout <<" Mass_h2mass in Analyzer " << std::endl;
	  float heavyMassEstimator_h2mass_prob = (heavyMassEstimator_h2mass->GetXaxis())->GetBinCenter(heavyMassEstimator_h2mass->GetMaximumBin());
	  float heavyMassEstimator_h2massweight1_prob = (heavyMassEstimator_h2mass_weight1->GetXaxis())->GetBinCenter(heavyMassEstimator_h2mass_weight1->GetMaximumBin());
	  float heavyMassEstimator_h2massweight4_prob = (heavyMassEstimator_h2mass_weight4->GetXaxis())->GetBinCenter(heavyMassEstimator_h2mass_weight4->GetMaximumBin());
	  float heavyMassEstimator_h2mass_RMS = heavyMassEstimator_h2mass->GetRMS();
	  float heavyMassEstimator_h2massweight1_RMS = heavyMassEstimator_h2mass_weight1->GetRMS();
	  float heavyMassEstimator_h2massweight4_RMS = heavyMassEstimator_h2mass_weight4->GetRMS();
	  //float heavyMassEstimator_h2mass_Entries = heavyMassEstimator_h2mass->GetEntries();
	  //float heavyMassEstimator_h2mass_Mean = heavyMassEstimator_h2mass->GetMean();
	  //int nbin=(heavyMassEstimator_h2mass->GetXaxis())->GetNbins();
	  //float heavyMassEstimator_h2mass_overflow = heavyMassEstimator_h2mass->GetBinContent(nbin+1);
	  //float heavyMassEstimator_h2mass_underflow = heavyMassEstimator_h2mass->GetBinContent(-1);
	  //float heavyMassEstimator_h2mass_MaxBin = heavyMassEstimator_h2mass->GetBinContent( heavyMassEstimator_h2mass->GetMaximumBin() );
	  //float heavyMassEstimator_h2mass_weight1_MaxBin = heavyMassEstimator_h2mass_weight1->GetBinContent( heavyMassEstimator_h2mass_weight1->GetMaximumBin() );
	  //float heavyMassEstimator_h2mass_weight4_MaxBin = heavyMassEstimator_h2mass_weight4->GetBinContent( heavyMassEstimator_h2mass_weight4->GetMaximumBin() );
	  std::cout <<"True HH mass "<< h2tohh_mass <<"; Reconstructed HH mass " << heavyMassEstimator_h2mass_prob <<" +/- "<< heavyMassEstimator_h2mass_RMS << "; reconstructed HH mass with type1 weight "<< heavyMassEstimator_h2massweight1_prob <<" +/- " << heavyMassEstimator_h2massweight1_RMS <<"; Reconstructed HH mass with type2 weight "<< heavyMassEstimator_h2massweight4_prob <<" +/- "<< heavyMassEstimator_h2massweight4_RMS <<std::endl;
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
	      //if (keepheavyMassEstimatorhist_)
	      //	file->WriteObject(heavyMassEstimator_h2mass, histname.c_str());
	}//runheavyMassEstimatorok
	delete thishme;

    }
    std::cout << std::endl;
    std::cout << "*****************************************************************************************************************************************" << std::endl;
    std::cout << "if you want to use this code, please cite:                                                                                               " << std::endl;
    std::cout << "T. Huang, J. M. No, L. Pernié, M. Ramsey-Musolf, A. Safonov, M. Spannowsky, and P. Winslow                                               " << std::endl;
    std::cout << "\" Resonant di-Higgs boson production in the bbWW channel: Probing the electroweak phase transition at the LHC\"                         " << std::endl;
    std::cout << "Phys. Rev. D 96, 035007 – Published 11 August 2017                                                                                       " << std::endl;
    std::cout << "*****************************************************************************************************************************************" << std::endl;
    std:: cout << std::endl;

}
