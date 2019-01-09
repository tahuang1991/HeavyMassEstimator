heavy mass estimator used ROOT package

how to use heavy mass estimator(HME) package
C++ version:
the HME class is defined in interface/heavyMassEstimator.h and src/heavyMassEstimator.cc
the constructor of HME is 

    heavyMassEstimator(TLorentzVector* lep1_lorentz, TLorentzVector* lep2_lorentz, TLorentzVector* b1jet_lorentz, TLorentzVector* b2jet_lorentz, 
	TLorentzVector* totjets_lorentz,TLorentzVector* met_lorentz, TLorentzVector* nu1_lorentz, TLorentzVector* nu2_lorentz,
	TLorentzVector* b_genp_lorentz, TLorentzVector* bbar_genp_lorentz, TLorentzVector* h2tohh_lorentz, int onshellMarker, bool simulation,	       bool PUsample_,
	int ievent, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
        int iterations, std::string RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose_=0
	);

in this constructor:  lep1_lorentz, lep2_lorentz are the lorentz vector for two leptons
b1jet_lorentz, b2jet_lorentz are the lorentz vector for two jets,
totjets_lorentz is the lorentz vector sum of all jets with pT > 20 and abs(eta)<2.5
nu1_lorentz, nu2_lorentz, b_genp_lorentz, bbar_genp_lorentz, h2tohh_lorentz, simulation contains simulation information and are used for debugging.
and the other constructor can ignore those inputs. 

PUsample_ is a boolean whether the input sample is PU or not. MET in HME is smeared according to MET resolution and MET resolution depends on PU
ievent is event id and used for debug purpose 
weightfromonshellnupt_func, weightfromonshellnupt_hist, weightfromonoffshellWmass_hist are booleans to define how to assign the weights to iteration results
and RefPDFfile is a root file containing the histograms used for weighing
iterations is the number of eta-phi pair random generation
useMET is whether to use MET from met collection or calculate MET from leptons and jets 
bjetrescaleAlgo and metcorrection are used to define how to correct bjets and met

the example code to use heavyMassEstimator is shown in testHME.cc and 10 events are hard coded in testHME.cc
use make to compile the package here and then run ./testHME.exe


in python version:
scripts/HeavyMassEstimator.py define the a class called HeavyMassEstimator
very similar to C++ version
HeavyMassEstimator used two leptons and two jets, MET as inputs to get HME results.
PUsample_ is by defualt is true and MET is taken from met collect in event (for useMET in C++ version)
bjetrescaleAlgo and metcorrection, by default, are taken the values that give the best performance in past



