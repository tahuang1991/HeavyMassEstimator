- [Heavy mass estimator(HME) in general](#heavy-mass-estimator)
- [Check HME package](#check-out-package)
- [HME in C++](#cplusplus-version)
- [HME in python](#python-version)
- [References](#references)

# Heavy mass estimator
Designed for estimating the mass of heavy resonance in X->HH->bbWW->bblvlv with two neutrinos in final states.  Both C++ and python version code are provided in this pacakge 

The HME version for X->HH->ZZbb->llvvbb is only implemented in python in python/HeavyMassEstimatorHHZZBB.py 
# Check out package
```
git clone https://github.com/tahuang1991/HeavyMassEstimator
```

# Cplusplus version
ROOT must be installed ahead as heavy mass estimator package depends on ROOT package 

## How to use heavy mass estimator c++ version package

the HME class is defined in interface/heavyMassEstimator.h and src/heavyMassEstimator.cc
the constructor of HME is 
>```
>    heavyMassEstimator(TLorentzVector* lep1_lorentz, TLorentzVector* lep2_lorentz, TLorentzVector* b1jet_lorentz, TLorentzVector* b2jet_lorentz, 
>	TLorentzVector* totjets_lorentz,TLorentzVector* met_lorentz, TLorentzVector* nu1_lorentz, TLorentzVector* nu2_lorentz,
>	TLorentzVector* b_genp_lorentz, TLorentzVector* bbar_genp_lorentz, TLorentzVector* h2tohh_lorentz, int onshellMarker, bool simulation,	       bool PUsample_,
>	int ievent, bool weightfromonshellnupt_func, bool weightfromonshellnupt_hist, bool weightfromonoffshellWmass_hist,
>        int iterations, std::string RefPDFfile, bool useMET, int bjetrescaleAlgo, int metcorrection, int verbose_=0
>	);
>```
in this constructor:  

  - lep1_lorentz, lep2_lorentz are the [lorentz vector](https://root.cern.ch/doc/v608/LorentzVectorPage.html) in ROOT for two leptons
     
  - b1jet_lorentz, b2jet_lorentz are the lorentz vectors for two jets,
  
  - totjets_lorentz is the lorentz vector sum of all jets with pT > 20 and abs(eta)<2.5
     
  - nu1_lorentz, nu2_lorentz, b_genp_lorentz, bbar_genp_lorentz, h2tohh_lorentz are from pure simulation and are used for debugging. The constructor can ignore those simluation inputs. 
  
  - PUsample_ is a boolean to determine whether the pileup (PU) of input sample is 40 or 0. MET in HME is smeared according to MET resolution and MET resolution depends on PU. If PUsample_=true, then MET resolution is 25.2 GeV, otherwise the MET resolution is 14.8 GeV
  
  - ievent is event id and used for debug purpose 
 
  - weightfromonshellnupt_func,weightfromonshellnupt_hist, weightfromonoffshellWmass_hist are booleans to define how to assign the weights to iteration results
 
  - RefPDFfile is a root file containing the histograms used for weighting and jet correction. The histograms are derived from delphes sample and used in Pheno study(https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.035007). These profiled PDFs could be generally used with CMS Run2 settings. 
 
  - iterations is the number of eta-phi pair random generations for one events.  The recommanded value is 10k.  With larger iterations, the HME could be slightly improved but it also comsumes more much computing resouce
 
  - useMET is whether to use MET from met collection or calculate MET from leptons and jets with momentum and energy conservation in transverse plane. useMET=True is recommended for Reconstuction events 
 
  - bjetrescaleAlgo and metcorrection are used to define how to correct bjets and met. bjetrescaleAlgo=2 and metcorrection = 5 are recommended from general studies.

In addition, it is also possible to set the MET covariance correction for each event and enable MET covriance matrix smearing
  >```
  >hme.setMETCovMatrix(met_covxx, met_covyy, met_covxy, True)
  >```




## How to compile and run C++ version heavy mass estimator
prerequisites to compile the package:
  -ROOT

To compile the code:
```
cd HeavyMassEstimator
make
```

To run heavy mass estimator with testHME:
```
cd exec/
./testHME
```
The above executable command is to run heavyMassEstimator with 10 events hard coded in testHME.cc.  The event with ievent=1 also enabled MET covariance matrix smearing in this example. 

To run heavy mass estimator with the provided Radion MC sample with mass = 700 GeV as input file
```
cd exec/
./runHME_ntuple ../data/Radion_M700_1kevents.root
```
The above command would read in the events from ROOT TTree "Events" in data/Radion_M700_1kevents.root and compute the HME for each event. The events with HME are stored in one output file, which by default is named "out_runHME_ntuple.root"




# Python version
python version HME code is under python/

python/HeavyMassEstimator.py defined a class called HeavyMassEstimator, which is very similar to C++ version.

## How to use python version HME class

First step is to initialze HME class
>```
>hme = HeavyMassEstimator()
>```

Second step is to customize the HME with kinemamtic inputs and other controling parameters

  - Resolved case: HeavyMassEstimator used two leptons and two jets, MET as inputs.
  >```
  >hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
  >```
  - Boosted case: HeavyMassEstimator used two leptons, the merged jet (H->bb) and MET as inputs 
  >```
  >hme.setKinematic_boosted(lep1_p4, lep2_p4, dijet_p4, met_vec2, 0.0)
  >```
  - Set the number of iteration for one event
  >```
  > hme.setIterations(iterations)
  >```
  - Set the MET resolution that is used to smeared MET. By default HME uses 25.2GeV as MET resolution 
  >```
  >hme.setMETResolution(metRes)
  >```
  - Set the MET covariance correction and enable MET covriance matrix smearing
  >```
  >hme.setMETCovMatrix(met_covxx, met_covyy, met_covxy, True)
  >```

Third step is to run HME
>```
>hme.runHME()
>```

And the last step is to get the most probable mass from HME
>```
>hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
>```

bjetrescaleAlgo and metcorrection in python version HME, by default, are using the following  values  compared to C++ version HME class:
  - PUSample = true, which means MET resolution is 25.2 GeV
  - useMET = true
  - bjetrescaleAlgo = 2
  - metcorrection = 5

The profiled PDFs are hard coded in python/HardcodeREFPDF.py which upgraded to include boosted case in python/HardcodeREFPDF_boosted.py

## Examples to use python version HME class

  - test/testHME.py runs the HME class with  one hard code event 
  - test/runHME_HHbbWW_general.py and test/runHME_HHbbWW_boosted.py read in the events from ROOT TTree and then compute the HME for each event. The following is example code from 
test/runHME_HHbbWW_general.py  
>```
>chain.GetEntry(nEv)
>lep1_p4 = ROOT.TLorentzVector()
>lep2_p4 = ROOT.TLorentzVector()
>jet1_p4 = ROOT.TLorentzVector()
>jet2_p4 = ROOT.TLorentzVector()
>dijet_p4 = ROOT.TLorentzVector()
>met_vec2 = ROOT.TVector2()
>
>lep1_p4.SetXYZT(chain.l1_px, chain.l1_py, chain.l1_pz, chain.l1_E)
>lep2_p4.SetXYZT(chain.l2_px, chain.l2_py, chain.l2_pz, chain.l2_E)
>jet1_p4.SetXYZT(chain.j1_px, chain.j1_py, chain.j1_pz, chain.j1_E)
>jet2_p4.SetXYZT(chain.j2_px, chain.j2_py, chain.j2_pz, chain.j2_E)
>met_vec2.SetX(chain.met_px)
>met_vec2.SetY(chain.met_py)
>
>hme = HeavyMassEstimator()
>hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
>hme.setIterations(iterations)
>hme.runHME()
>hme_mass_peak[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
>```
and here is the example about how to run runHME_HHbbWW_boosted.py with Radion_M700_1kevents.root as input
```
cd test/
python runHME_HHbbWW_boosted.py -i ../data/Radion_M700_1kevents.root -o Radion_M700_HME.root -it 10000
```
The output file to store events with HME is Radion_M700_HME.root and also the iterations for one event is 10000.
# References

T. Huang, J. M. No, L. Pernié, M. Ramsey-Musolf, A. Safonov, M. Spannowsky, and P. Winslow                                               
" Resonant di-Higgs boson production in the bbWW channel: Probing the electroweak phase transition at the LHC"                         
Phys. Rev. D 96, 035007 – Published 11 August 2017  
https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.035007


