import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
from math import *
from HeavyMassEstimator import *
import argparse
import numpy as np
ROOT.gROOT.SetBatch(1)
#sys.argv.append( '-b' ) 
from array import array
import mt2_bisect

import Samplelist_2016Run as Samplelist
import warnings
sys.argv.append( '-b' )
sys.argv.append( '-q' )

warnings.filterwarnings("ignore", category=RuntimeWarning)
#warnings.filterwarnings("ignore", category=RuntimeWarning)
#with warnings.catch_warnings():
        #warnings.simplefilter('ignore', RuntimeWarning)
execfile("start.py")
execfile("functions.py")


sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD')
from localSamplelist import *
#Creating folders and parameters
def AddFiles(chain, inputDir):
    theInputFiles = []
    for d in range(len(inputDir)):
        my_dir = inputDir[d]
        if not os.path.isdir(my_dir):
            print "ERROR: This is not a valid directory: ", my_dir
            if d==len(inputDir)-1:
                print "ERROR: No input files were selected"
                exit()
            continue
        print "Proceed to next directory"
        ls = os.listdir(my_dir)
        chain.Add(my_dir[:] + x for x in ls if x.endswith('root'))

doTest = False
doHME = True

#tree_name="DiHiggsWWBBAna/evtree"
#tree_name = "t"
tree_name = "Friends"
benchmarks = ["ttV", "Wjet", "sT", "VV", "DY", "TT", "Data"]

#whichSample = "Graviton_0419"
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
	return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
	return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-n", "--njobs", dest="njobs", type=int, default=100, help="total splitted jobs. [Default: 100]")
parser.add_argument("-i", "--ijob", dest="ijob", type=int, default=0, help="ith job in splitting. [Default: 0]")
parser.add_argument("-jt", "--jobtype", dest="jobtype", default="RadionM300", help="sample type to run. [Default: radion]")
parser.add_argument("-d", "--dataname", dest="dataname", default=None, help="sample type to run. [Default: None]")
parser.add_argument("-it", "--iterations", dest="iterations", type=int, default=10000, help="iteration number used HME [Default: 100000]")
parser.add_argument("-o", "--outputdir",dest="output", default="./",help="output root file directory, [Defualt:PWD] ")
parser.add_argument("-DY", "--DYestimation",dest="DYestimation",type=str2bool, default=False,help="do DYestimation or not, [Defualt:False] ")
#args = parser.parse_args()
args, unknown = parser.parse_known_args()
if args.output[-1] != '/':
   args.output = args.output+'/'
print " output ",args.output
whichSample = args.jobtype
if whichSample not in benchmarks and "Radion" not in whichSample and "Graviton" not in whichSample and "untagged" not in whichSample:
    print "error!!! ",whichSample," is not in ",benchmarks
    sys.exit()

if len(full_local_samplelist[args.jobtype].keys()) == 1:
    args.dataname = full_local_samplelist[args.jobtype].keys()[0]
elif args.dataname == None:
     print "error!!! no dataname is found "
     sys.exit()

maxn = 10 
#DY_BDT_flat                 = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
#dy_nobtag_to_btagM_weight   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
passCC                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
passpre                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
mt2                         = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
mt2_bb                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
mt2_ll			    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_gen              = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_reco             = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_reco               = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_reco             = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entries_reco            = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entry_peak_reco         = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_weight2_gen      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_weight2_reco     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_weight2_reco       = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_weight2_reco     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entries_weight2_reco    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entry_peak_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
#hme_offshellWmass_reco      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mostprob_offshellWmass_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_offshellWmass_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_offshellWmass_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mostprob_offshellWmass_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_offshellWmass_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_offshellWmass_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hmecputime                  = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
cross_section               = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
event_weight_sum            = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_bins = []
for i in range(100):
    hme_bins.append(array( 'f', maxn*[ 0. ] ))



thisfile = full_local_samplelist[args.jobtype][args.dataname]["path"]
if args.DYestimation:
    #thisfile = untagged_samplelist[args.jobtype][args.dataname]["path"]
    thisfile = full_local_samplelist[args.jobtype][args.dataname]["path"]


if not os.path.isfile(thisfile):
    print "error!!! not valid file ",thisfile
    exit()

if args.jobtype != "Data" and args.jobtype != "Data_untagged":
    tfiletmp = ROOT.TFile(thisfile, "READ") 
    h_cutflow = tfiletmp.Get("h_cutflow")
    event_weight_sum[0] = h_cutflow.GetBinContent(1)
    cross_section[0] = full_local_samplelist[args.jobtype][args.dataname]["cross_section"]##FIXME
    tfiletmp.Close()
    print "sample ",args.jobtype," dataset ",args.dataname ," cross_section[0] ",cross_section[0]," event_weight_sum[0]  ",event_weight_sum[0]," file ",thisfile

TCha = ROOT.TChain(tree_name)
TCha.Add(thisfile)
#TCha.Add("/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root")

TotalEv = TCha.GetEntries()
nStart = (TotalEv/args.njobs + 1)*args.ijob
nEnd = (TotalEv/args.njobs + 1)*(args.ijob+1)
if args.ijob == args.njobs-1 or nEnd>TCha.GetEntries():
    nEnd = TCha.GetEntries()
print whichSample, "TChain has", TCha.GetEntries(), "entries."," job star from ",nStart," end at ",nEnd, " events to run ",nEnd-nStart
#f = ROOT.TFile('/fdata/hepx/store/user/%s/Hhh_For_Test/'%user + whichSample + '_ijob%d.root'%args.ijob,'recreate'); f.cd()
f = ROOT.TFile(args.output + args.dataname + '_ijob%d.root'%args.ijob, 'recreate'); f.cd()
TCha2 = TCha.CloneTree(0)


#TCha2.Branch("DY_BDT_flat",                   DY_BDT_flat,                  "DY_BDT_flat/F")
#TCha2.Branch("dy_nobtag_to_btagM_weight",     dy_nobtag_to_btagM_weight,    "dy_nobtag_to_btagM_weight/F")
TCha2.Branch("passCC",                        passCC,                       "passCC/F")
TCha2.Branch("passpre",                       passpre,                      "passpre/F")
TCha2.Branch("mt2",                           mt2,                          "mt2/F")
TCha2.Branch("mt2_bb",                        mt2_bb,                       "mt2_bb/F")
TCha2.Branch("mt2_ll",                        mt2_ll,                       "mt2_ll/F")
TCha2.Branch("hme_h2mass_gen",                hme_h2mass_gen,               "hme_h2mass_gen/F")
TCha2.Branch("hme_h2mass_reco",               hme_h2mass_reco,              "hme_h2mass_reco/F")
TCha2.Branch("hme_mean_reco",                 hme_mean_reco,                "hme_mean_reco/F")
TCha2.Branch("hme_stddev_reco",               hme_stddev_reco,              "hme_stddev_reco/F")
TCha2.Branch("hme_entries_reco",              hme_entries_reco,             "hme_entries_reco/F")
TCha2.Branch("hme_entry_peak_reco",           hme_entry_peak_reco,          "hme_entry_peak_reco/F")
TCha2.Branch("hme_h2mass_weight2_gen",        hme_h2mass_weight2_gen,       "hme_h2mass_weight2_gen/F")
TCha2.Branch("hme_h2mass_weight2_reco",       hme_h2mass_weight2_reco,      "hme_h2mass_weight2_reco/F")
TCha2.Branch("hme_mean_weight2_reco",         hme_mean_weight2_reco,        "hme_mean_weight2_reco/F")
TCha2.Branch("hme_stddev_weight2_reco",        hme_stddev_weight2_reco,      "hme_stddev_weight2_reco/F")
TCha2.Branch("hme_entries_weight2_reco",      hme_entries_weight2_reco,     "hme_entries_weight2_reco/F")
TCha2.Branch("hme_entry_peak_weight2_reco",   hme_entry_peak_weight2_reco,  "hme_entry_peak_weight2_reco/F")
TCha2.Branch("hme_mostprob_offshellWmass_reco",   hme_mostprob_offshellWmass_reco,  "hme_mostprob_offshellWmass_reco/F")
TCha2.Branch("hme_mean_offshellWmass_reco",   hme_mean_offshellWmass_reco,  "hme_mean_offshellWmass_reco/F")
TCha2.Branch("hme_stddev_offshellWmass_reco", hme_stddev_offshellWmass_reco,"hme_stddev_offshellWmass_reco/F")
TCha2.Branch("hme_mean_offshellWmass_weight2_reco",   hme_mean_offshellWmass_weight2_reco,  "hme_mean_offshellWmass_weight2_reco/F")
TCha2.Branch("hme_mostprob_offshellWmass_weight2_reco",   hme_mostprob_offshellWmass_weight2_reco,  "hme_mostprob_offshellWmass_weight2_reco/F")
TCha2.Branch("hme_stddev_offshellWmass_weight2_reco", hme_stddev_offshellWmass_weight2_reco,"hme_stddev_offshellWmass_weight2_reco/F")
TCha2.Branch("hmecputime",                    hmecputime,                   "hmecputime/F")
TCha2.Branch("cross_section",                 cross_section,                "cross_section/F")
TCha2.Branch("event_weight_sum",              event_weight_sum,             "event_weight_sum/F")

#ignore HME bins now 
#for i in range(100):
#    TCha2.Branch("hme_bin%d"%i, hme_bins[i], "hme_bin%d/F"%(i))

h_h2mass_weight_gen      = ROOT.TH1F("h_h2mass_weight_gen","",100,200.,1200.);  h_h2mass_weight_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_gen_sum  = ROOT.TH1F("h_h2mass_weight_gen_sum","",1000,200.,1200.);  h_h2mass_weight_gen_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_gen     = ROOT.TH1F("h_h2mass_weight1_gen","",100,200.,1200.);  h_h2mass_weight1_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight2_gen     = ROOT.TH1F("h_h2mass_weight2_gen","",100,200.,1200.);  h_h2mass_weight2_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco     = ROOT.TH1F("h_h2mass_weight_reco","",100,200.,1200.);  h_h2mass_weight_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco_sum = ROOT.TH1F("h_h2mass_weight_reco_sum","",1000,200.,1200.);  h_h2mass_weight_reco_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_reco    = ROOT.TH1F("h_h2mass_weight1_reco","",100,200.,1200.);  h_h2mass_weight1_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight2_reco     = ROOT.TH1F("h_h2mass_weight2_reco","",100,200.,1200.);  h_h2mass_weight2_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
offshellWmass_gen	 = ROOT.TH1F("hme_offshellWmass_gen","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_gen.GetXaxis().SetTitle("HME offshell W mass [GeV]");
offshellWmass_reco	 = ROOT.TH1F("hme_offshellWmass_reco","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_reco.GetXaxis().SetTitle("HME offshell W mass [GeV]");
# adding more branches to TTree
def initbr():
    passCC[0] = 0
    passpre[0] = 0
    hme_h2mass_gen[0] = -1
    hme_h2mass_reco[0] = -1
    hme_mean_reco[0] = -1
    hme_stddev_reco[0] = -1
    hme_entry_peak_reco[0] = -1
    hme_h2mass_weight2_reco[0] = -1
    hme_mean_weight2_reco[0] = -1
    hme_stddev_weight2_reco[0] = -1
    hme_entry_peak_weight2_reco[0] = -1
    #print "hme_bins ",hme_bins
    hme_entries_reco[0] = -1
    for i in range(len(hme_bins)):
	hme_bins[i][0] = 0.0

stop_watch = ROOT.TStopwatch()
stop_watch2 = ROOT.TStopwatch()

stop_watch2.Start()
for nEv in range(nStart, nEnd):
  if (doTest and nEv%10 == 0 ):
      print "nEv ",nEv
  elif (nEv%10000 == 0):
      print "nEv ",nEv
  if (doTest and nEv-nStart>=100):
      break

  print "nEv ",nEv
  if nEv == nStart:
      print "First event to run: ",nEv, " last event ", nEnd," HME iterations ",args.iterations
  TCha.GetEntry(nEv)
  initbr()

  #if whichSample.startswith("DY"):
  #    DY_BDT_flat[0]    =  TCha.DY_BDT_flat
  #    dy_nobtag_to_btagM_weight[0] = TCha.dy_nobtag_to_btagM_weight
  #else:
  #   dy_nobtag_to_btagM_weight[0] = 1.0 

  MET_cut               = (TCha.met_pt>20)
  MuMu_cut              = (TCha.isMuMu and TCha.lep1_pt>20 and fabs(TCha.lep1_eta)<2.4 and TCha.lep2_pt>10 and fabs(TCha.lep2_eta)<2.4 and TCha.ll_M>12)
  MuEl_cut              = (TCha.isMuEl and TCha.lep1_pt>20 and fabs(TCha.lep1_eta)<2.4 and TCha.lep2_pt>15 and fabs(TCha.lep2_eta)<2.5 and TCha.ll_M>12)
  ElMu_cut              = (TCha.isElMu and TCha.lep1_pt>25 and fabs(TCha.lep1_eta)<2.5 and TCha.lep2_pt>10 and fabs(TCha.lep2_eta)<2.4 and TCha.ll_M>12)
  ElEl_cut              = (TCha.isElEl and TCha.lep1_pt>25 and fabs(TCha.lep1_eta)<2.5 and TCha.lep2_pt>15 and fabs(TCha.lep2_eta)<2.5 and TCha.ll_M>12)
  #B1B2_cut              = (TCha.jet1_pt>20 and fabs(TCha.jet1_eta)<2.4 and TCha.jet2_pt>20 and fabs(TCha.jet2_eta)<2.4)
  B1B2_cut              = (TCha.jet1_pt>20 and TCha.jet2_pt>20)
  #preselection          = (MET_cut and MuMu_cut and B1B2_cut)
  preselection          = (MET_cut)### Louvain sample, MuMuCut is okay, B1B2Cut is okay 
  mt_cleancut           = (TCha.llmetjj_MTformula>10)
  NoZ_cleancut          = (TCha.ll_M<(91-15))
  if not NoZ_cleancut:
      TCha2.Fill()
      continue
  mJJ_cleancut          = (TCha.jj_M>30)
  #dR_lljj_cleancut      = (TCha.dR_l1l2b1b2>0.2 and TCha.dR_l1l2b1b2<4)
  dR_lljj_cleancut      = 1
  dR_ll_cleancut        = (TCha.ll_DR_l_l<3.8)
  dR_jj_cleancut        = (TCha.jj_DR_j_j<3.8)
  cleaning_cuts         = (mt_cleancut and NoZ_cleancut and mJJ_cleancut and dR_lljj_cleancut and dR_ll_cleancut and dR_jj_cleancut)
  #if (preselection and cleaning_cuts ):
  #    print "event ",nEv," pass CC cut"
  if (preselection):
      passpre[0] = 1
  if (cleaning_cuts):
      passCC[0] = 1
  

  #if (cleaning_cuts and TCha.findAllGenParticles and doHME and False):
  if (False):
      
      mu1_p4 	      = ROOT.TLorentzVector(TCha.mu1_px, TCha.mu1_py, TCha.mu1_pz, TCha.mu1_energy)
      mu2_p4 	      = ROOT.TLorentzVector(TCha.mu2_px, TCha.mu2_py, TCha.mu2_pz, TCha.mu2_energy)
      jet1_p4 	      = ROOT.TLorentzVector(TCha.b1_px, TCha.b1_py, TCha.b1_pz, TCha.b1_energy)
      jet2_p4 	      = ROOT.TLorentzVector(TCha.b2_px, TCha.b2_py, TCha.b2_pz, TCha.b2_energy)
      met_vec2        = ROOT.TVector2(TCha.nu1_px + TCha.nu2_px, TCha.nu1_py + TCha.nu2_py)
      hme_gen = HeavyMassEstimator()
      onshellW_mu  = 1 
      if TCha.w1_mass < TCha.w2_mass:
      	  onshellW_mu  = 2
      hme_gen.setKinematic(mu1_p4, mu2_p4, jet1_p4, jet2_p4, met_vec2, onshellW_mu)
      #hme_gen.setonshellWmasspdf(onshellWmasspdf)
      #hme_gen.setoffshellWmasspdf(offshellWmasspdf)
      #hme_gen.setonshellnuptpdf(onshellnuptpdf)
      #hme.showKinematic()
      hme_gen.setIterations(args.iterations)
      hme_gen.runHME()
      if hme_gen.hme_h2Mass.GetEntries()>0 and hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())>200:
	  print "Gen level most probable reco mass ",hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())," entries ",hme_gen.hme_h2Mass.GetEntries()
	  hme_gen.hme_h2Mass.SetName("hme_h2Mass_ev%d_genlevel"%nEv)
          hme_gen.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_genlevel"%nEv)
	  hme_gen.hmetree.SetName("hmetree_ev%d_genlevel"%nEv)
	  hme_gen.hme_h2Mass_correctmunupair.SetName("hme_h2Mass_correctmunupair_ev%d_genlevel"%nEv)
	  hme_gen.hme_h2Mass_incorrectmunupair.SetName("hme_h2Mass_incorrectmunupair_ev%d_genlevel"%nEv)
	  hme_gen.hme_h2Mass_correctmunupair.SetLineColor(ROOT.kRed)
	  hme_gen.hme_h2Mass_incorrectmunupair.SetLineColor(ROOT.kBlue)
    	  hme_gen.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_genlevel"%nEv)
	  h_h2mass_weight_gen_sum.Add(hme_gen.hme_h2Mass)
          offshellWmass_gen.Add(hme_gen.hme_offshellWmass)
	  if nEv%100 == 0:
	      hme_gen.hme_h2Mass.Write()
	      hme_gen.hme_h2MassWeight1.Write()
    	      hme_gen.hme_h2Mass_correctmunupair.Write()
    	      hme_gen.hme_h2Mass_incorrectmunupair.Write()
    	      hme_gen.hme_offshellWmass.Write()
              hme_gen.hmetree.Write()
          h_h2mass_weight_gen.Fill(hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin()))
          h_h2mass_weight1_gen.Fill(hme_gen.hme_h2MassWeight1.GetXaxis().GetBinCenter(hme_gen.hme_h2MassWeight1.GetMaximumBin()))
          hme_h2mass_gen[0] = hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())
  
  #if (cleaning_cuts and doHME):
  stop_watch.Start()
  if (doHME):
      #lep1_p4 	  = ROOT.TLorentzVector(TCha.lep1_px, TCha.lep1_py, TCha.lep1_pz, TCha.lep1_energy)
      #lep2_p4 	  = ROOT.TLorentzVector(TCha.lep2_px, TCha.lep2_py, TCha.lep2_pz, TCha.lep2_energy)
      #jet1_p4 	  = ROOT.TLorentzVector(TCha.jet1_px, TCha.jet1_py, TCha.jet1_pz, TCha.jet1_energy)
      #jet2_p4 	  = ROOT.TLorentzVector(TCha.jet2_px, TCha.jet2_py, TCha.jet2_pz, TCha.jet2_energy)
      #met_vec2    = ROOT.TVector2(TCha.met_px, TCha.met_py)
      lep1_p4 	  = ROOT.TLorentzVector(); lep1_p4.SetPtEtaPhiM(TCha.lep1_pt, TCha.lep1_eta, TCha.lep1_phi, 0)
      lep2_p4 	  = ROOT.TLorentzVector(); lep2_p4.SetPtEtaPhiM(TCha.lep2_pt, TCha.lep2_eta, TCha.lep2_phi, 0) 
      #jet1_p4 	  = ROOT.TLorentzVector(); jet1_p4.SetPtEtaPhiM(TCha.jet1_pt, TCha.jet1_eta, TCha.jet1_phi, 4.18)
      #jet2_p4 	  = ROOT.TLorentzVector(); jet2_p4.SetPtEtaPhiM(TCha.jet2_pt, TCha.jet2_eta, TCha.jet2_phi, 4.18)
      jet1_p4 	  = ROOT.TLorentzVector(); jet1_p4.SetPtEtaPhiE(TCha.jet1_pt, TCha.jet1_eta, TCha.jet1_phi, TCha.jet1_E)
      jet2_p4 	  = ROOT.TLorentzVector(); jet2_p4.SetPtEtaPhiE(TCha.jet2_pt, TCha.jet2_eta, TCha.jet2_phi, TCha.jet2_E)
      met_vec2    = ROOT.TVector2();       met_vec2.SetMagPhi(TCha.met_pt, TCha.met_phi)
      
      mt2_event1 =  mt2_bisect.mt2()
      sumesBl11_p4 = lep1_p4+jet1_p4
      sumesBl22_p4 = lep2_p4+jet2_p4
      mt2_event1.set_momenta(sumesBl11_p4.M(), sumesBl11_p4.Px(), sumesBl11_p4.Py(), sumesBl22_p4.M(), sumesBl22_p4.Px(), sumesBl22_p4.Py(), 0.0, met_vec2.Px(), met_vec2.Py())
      mt2_event1.set_mn(0.0)
      MT2_1 = mt2_event1.get_mt2()

      mt2_event2 =  mt2_bisect.mt2()
      sumesBl12_p4 = lep2_p4+jet1_p4
      sumesBl21_p4 = lep1_p4+jet2_p4
      mt2_event2.set_momenta(sumesBl12_p4.M(), sumesBl12_p4.Px(), sumesBl12_p4.Py(), sumesBl21_p4.M(), sumesBl21_p4.Px(), sumesBl21_p4.Py(), 0.0, met_vec2.Px(), met_vec2.Py())
      mt2_event2.set_mn(0.0)
      MT2_2 = mt2_event2.get_mt2()
      
      mt2[0] = min(MT2_1, MT2_2)

      #jj_p4 = jet1_p4+jet2_p4
      #ll_p4 = lep1_p4+lep2_p4
      mt2_eventbb = mt2_bisect.mt2()
      mt2_eventbb.nevt = nEv
      mt2_eventbb.set_momenta(jet1_p4.M(), jet1_p4.Px(), jet1_p4.Py(), jet2_p4.M(), jet2_p4.Px(), jet2_p4.Py(), 0.0, met_vec2.Px()+lep1_p4.Px()+lep2_p4.Px(), met_vec2.Py()+lep1_p4.Py()+lep1_p4.Py())
      mt2_eventbb.set_mn(80.4)
      #mt2_eventbb._print()
      mt2_bb[0] = mt2_eventbb.get_mt2()

      mt2_eventll = mt2_bisect.mt2()
      mt2_eventll.nevt = nEv
      mt2_eventll.set_momenta(lep1_p4.M(), lep1_p4.Px(), lep1_p4.Py(), lep2_p4.M(), lep2_p4.Px(), lep2_p4.Py(), 0.0, met_vec2.Px(), met_vec2.Py())
      mt2_eventll.set_mn(0.0)
      mt2_ll[0] = mt2_eventll.get_mt2()
      #mt2_eventll._print()

      #print "MT2_1 ", MT2_1," MT2_2 ",MT2_2, " MT2 final ", mt2[0]," mt2_ll ",mt2_ll[0]," mt2_jj ",mt2_bb[0]
      hme = HeavyMassEstimator()
      hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
      #hme.setonshellWmasspdf(onshellWmasspdf)
      #hme.setoffshellWmasspdf(offshellWmasspdf)
      #hme.setonshellnuptpdf(onshellnuptpdf)
      #hme.setrecobjetrescalec1pdf(recobjetrescalec1pdfPU40)
      #hme.showKinematic()
      hme.setIterations(args.iterations)
      #hme.setDebug(True)
      hme.runHME()
      #hme.hme_offshellWmass.SetName("hme_offshellWmass_TCha.d_genlTCha.e"%nEv)
      if hme.hme_h2Mass.GetEntries() <= 0:
          print "NO solution found!!!!! "
      elif hme.hme_h2Mass.GetEntries() >0 and  hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin()) < 249.0 :
          print "Num solutions ",hme.hme_h2Mass.GetEntries()," BUT the maximum is ",hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
	  #hme.hme_h2Mass.Print("ALL")

      if hme.hme_h2Mass.GetEntries()>0 and hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())>=250.0 :
	  #print "Reco Level most probable reco mass ",hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())," entries ",hme.hme_h2Mass.GetEntries()," stddev ",hme.hme_h2Mass.GetStdDev(1)
	  hme.hme_h2Mass.SetName("hme_h2Mass_ev%d_recolevel"%nEv)
          hme.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_recolevel"%nEv)
          hme.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_recolevel"%nEv)
	  h_h2mass_weight_reco_sum.Add(hme.hme_h2Mass)
	  #if nEv%100 == 0:
	  #    hme.hme_h2Mass.Write()
	  #    hme.hme_h2MassWeight1.Write()
    	  #    hme.hme_offshellWmass.Write()
          h_h2mass_weight_reco.Fill(hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin()))
          h_h2mass_weight1_reco.Fill(hme.hme_h2MassWeight1.GetXaxis().GetBinCenter(hme.hme_h2MassWeight1.GetMaximumBin()))
          h_h2mass_weight2_reco.Fill(hme.hme_h2MassWeight2.GetXaxis().GetBinCenter(hme.hme_h2MassWeight2.GetMaximumBin()))

    	  hme_h2mass_reco[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
          hme_mean_reco[0] = hme.hme_h2Mass.GetMean()
          #hme_stddev_reco[0] = hme.hme_h2Mass.GetStdDev(1)
          hme_entries_reco[0] = float(hme.hme_h2Mass.GetEntries())/args.iterations
	  hme_entry_peak_reco[0] = hme.hme_h2Mass.Integral(hme.hme_h2Mass.GetMaximumBin()-5, hme.hme_h2Mass.GetMaximumBin()+5)
    	  hme_h2mass_weight2_reco[0] = hme.hme_h2MassWeight2.GetXaxis().GetBinCenter(hme.hme_h2MassWeight2.GetMaximumBin())
          hme_mean_weight2_reco[0] = hme.hme_h2MassWeight2.GetMean()
          #hme_stddev_weight2_reco[0] = hme.hme_h2MassWeight2.GetStdDev(1)
          hme_entries_weight2_reco[0] = float(hme.hme_h2MassWeight2.GetEntries())/args.iterations
	  hme_entry_peak_weight2_reco[0] = hme.hme_h2MassWeight2.Integral(hme.hme_h2MassWeight2.GetMaximumBin()-5, hme.hme_h2MassWeight2.GetMaximumBin()+5)
	  print "hme_h2mass_reco[0] ",hme_h2mass_reco[0]," hme_h2mass_weight2_reco[0] ",hme_h2mass_weight2_reco[0]
	  #offshell Wmass
	  h_offshellWmass_recoh2mass = hme.hme_h2MassAndoffshellWmass.ProjectionY("h_offshellWmass_recoh2mass",hme.hme_h2Mass.GetMaximumBin()-5, hme.hme_h2Mass.GetMaximumBin()+5)
	  h_offshellWmass_recoh2mass_weight2 = hme.hme_h2MassAndoffshellWmass_weight2.ProjectionY("h_offshellWmass_recoh2mass_weight2",hme.hme_h2MassWeight2.GetMaximumBin()-5, hme.hme_h2MassWeight2.GetMaximumBin()+5)
	  hme_mostprob_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetXaxis().GetBinCenter(h_offshellWmass_recoh2mass.GetMaximumBin())
	  hme_mean_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetMean()
	  #hme_stddev_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetStdDev(1)
	  hme_mostprob_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetXaxis().GetBinCenter(h_offshellWmass_recoh2mass_weight2.GetMaximumBin())
	  hme_mean_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetMean()
	  #hme_stddev_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetStdDev(1)

          hme.hme_h2Mass.Rebin(10)
          hme.hme_h2MassWeight2.Rebin(10)
          nbins = hme.hme_h2Mass.GetXaxis().GetNbins()
          #print "nbins ",nbins
	  if nbins != len(hme_bins):
	      print "error!!! nbins of hme_h2mass is not the same as the one written to TTree "
	  for i in range(1, len(hme_bins)+1):
	      hme_bins[i-1][0] =  hme.hme_h2Mass.GetBinContent(i)
        
  # WEIGHTS
  hmecputime[0] =  stop_watch.CpuTime()
  #print "Cputime ",stop_watch.CpuTime()," realtime ",stop_watch.RealTime()
  stop_watch.Stop()
  TCha2.Fill()

TCha2.Write()
f.Close()
print "Done. %s ijob %d "%( args.jobtype,  args.ijob)
stop_watch2.Stop()
print "stop_watch2 Cputime ",stop_watch2.CpuTime()," realtime ",stop_watch2.RealTime()
