import ROOT
from array import array
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.argv.append( '-b' )
sys.argv.append( '-q' )

sys.path.append(os.path.join(os.path.abspath('..'), 'python'))
from HeavyMassEstimator import *

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-i", "--inputFile", dest="infile", type=str, default=None, help="input file name. [Default: None]")
parser.add_argument("-o", "--outputFile", dest="outfile", type=str, default="out.root", help="output file name. [Default: out.root]")
parser.add_argument("-it", "--iterations", dest="iterations", type=int, default=10000, help="iteration number used HME [Default: 100000]")
parser.add_argument("-nrebins", "--nrebins", dest="nrebins", type=int, default=10, help="rebining the histogram from HME [Default: 1]")
parser.add_argument("-mR", "--metRes", dest="metRes", type=float, default=None, help="MET resolution [Default: None, would value hard coded in HME]")
parser.add_argument("-gp", "--gravitypercent", dest="gravitypercent", type=float, default=0.18, help="gravitypercent [Default: 0.18]")
parser.add_argument("-nStart", "--nStart", dest="nStart", type=int, default=0, help="first event index to process. [Default: 0]")
parser.add_argument("-nEnd", "--nEnd", dest="nEnd", type=int, default=-1, help="last event index to process. -1 means last event in TTree [Default: -1]")
#args = parser.parse_args()
args, unknown = parser.parse_known_args()

if args.infile == None or not(os.path.isfile(args.infile)):
  print("No input file given ")
  sys.exit("Use 'python runHME_HHbbWW_boosted.py -i <inputFile> -o <outputFile>'")


nbins_rebin = args.nrebins
metRes = args.metRes
per_weight_gravitycenter = args.gravitypercent

def get_gravity_center(th1, percent):
    ##get the gravity center(aka weight average) around peak
    binwidth = th1.GetBinWidth(1)
    bin_max = th1.GetMaximumBin()
    bin_max_center = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin())
    if th1.GetEntries() == 0:
        return bin_max_center
    nbins_gravity = int(bin_max_center*percent/binwidth/2.0)
    total_weightedEntries = th1.GetBinContent(bin_max) * bin_max_center
    total_entries = th1.GetBinContent(bin_max) 
    for i in range(1, nbins_gravity+1):
      if bin_max - i >0:
        total_weightedEntries += th1.GetBinContent(bin_max - i) * th1.GetBinCenter(bin_max - i)
        total_entries += th1.GetBinContent(bin_max - i)
      if bin_max + i <= th1.GetNbinsX():
        total_weightedEntries += th1.GetBinContent(bin_max + i) * th1.GetBinCenter(bin_max + i)
        total_entries += th1.GetBinContent(bin_max + i)
    if total_entries == 0:
        return bin_max_center
    weighted_average = total_weightedEntries/total_entries

    if (bin_max_center - 250)*(weighted_average - 250) < 0.0:
      print("percent ", percent, " nbins for gravity on each side ", nbins_gravity, " max bin ",bin_max," most probable mass ", bin_max_center, " total_entries ", total_entries, " weighted sum ", total_weightedEntries, " New weighted probable mass ", weighted_average)
    return weighted_average

chain = ROOT.TChain("Events")
chain.Add(args.infile)
TotalEv = chain.GetEntries()

out = ROOT.TFile(args.outfile, 'recreate')
out.cd()
tree = chain.CloneTree(0)

iterations = args.iterations
nStart = args.nStart
nEnd = TotalEv
if args.nEnd > 0 and args.nEnd <= TotalEv:
  nEnd = args.nEnd
if nStart >= nEnd:
  sys.exit("Error! first event index to process is larger than last one: nStart ", nStart, " nEnd ", nEnd)
print("Total events = ", TotalEv, " nStart ", nStart, " nEnd ", nEnd, " iteration for HME ", iterations," gravitypercent ", args.gravitypercent)

hme_mass_peak                = array('f', [0.])
hme_mass_peak_gravity        = array('f', [0.])
hme_mass_peak_divSol         = array('f', [0.])
hme_mass_peak_boosted        = array('f', [0.])
hme_mass_peak_boosted_divSol = array('f', [0.])
nsolution0 = array('f', [0.0])
nsolution1 = array('f', [0.0])
nsolution2 = array('f', [0.0])
nsolution3 = array('f', [0.0])
nsolution4 = array('f', [0.0])

massBranch                = tree.Branch("hme_mass_peak", hme_mass_peak, "hme_mass_peak/F")
massBranch_gravity        = tree.Branch("hme_mass_peak_gravity", hme_mass_peak_gravity, "hme_mass_peak_gravity/F")
massBranch_divSol         = tree.Branch("hme_mass_peak_divSol", hme_mass_peak_divSol, "hme_mass_peak_divSol/F")
massBranch_boosted        = tree.Branch("hme_mass_peak_boosted", hme_mass_peak_boosted, "hme_mass_peak_boosted/F")
massBranch_boosted_divSol = tree.Branch("hme_mass_peak_boosted_divSol", hme_mass_peak_boosted_divSol, "hme_mass_peak_boosted_divSol/F")
br_nsolution0             = tree.Branch("hme_nsolution0", nsolution0, "nsolution0/F")
br_nsolution1             = tree.Branch("hme_nsolution1", nsolution1, "nsolution1/F")
br_nsolution2             = tree.Branch("hme_nsolution2", nsolution2, "nsolution2/F")
br_nsolution3             = tree.Branch("hme_nsolution3", nsolution3, "nsolution3/F")
br_nsolution4             = tree.Branch("hme_nsolution4", nsolution4, "nsolution4/F")

def initBranches():
  hme_mass_peak[0] = 0.0
  hme_mass_peak_divSol[0] = 0.0
  hme_mass_peak_gravity[0] = 0.0
  hme_mass_peak_boosted[0] = 0.0
  hme_mass_peak_boosted_divSol[0] = 0.0
  nsolution0[0] = 0.0
  nsolution1[0] = 0.0
  nsolution2[0] = 0.0
  nsolution3[0] = 0.0
  nsolution4[0] = 0.0

for nEv in range(nStart, nEnd):
  initBranches()
  chain.GetEntry(nEv)
  lep1_p4 = ROOT.TLorentzVector()
  lep2_p4 = ROOT.TLorentzVector()
  jet1_p4 = ROOT.TLorentzVector()
  jet2_p4 = ROOT.TLorentzVector()
  dijet_p4 = ROOT.TLorentzVector()
  met_vec2 = ROOT.TVector2()

  lep1_p4.SetPxPyPzE(chain.l1_px,   chain.l1_py,    chain.l1_pz,    chain.l1_E)
  lep2_p4.SetPxPyPzE(chain.l2_px,   chain.l2_py,    chain.l2_pz,    chain.l2_E)
  jet1_p4.SetPxPyPzE(chain.j1_px,   chain.j1_py,    chain.j1_pz,    chain.j1_E)
  jet2_p4.SetPxPyPzE(chain.j2_px,   chain.j2_py,    chain.j2_pz,    chain.j2_E)
  #jet1_p4.SetPxPyPzE(chain.fatbjet_subjet1_Px,   chain.fatbjet_subjet1_Py,    chain.fatbjet_subjet1_Pz,    chain.fatbjet_subjet1_E)
  #jet2_p4.SetPxPyPzE(chain.fatbjet_subjet2_Px,   chain.fatbjet_subjet2_Py,    chain.fatbjet_subjet2_Pz,    chain.fatbjet_subjet2_E)
  met_vec2.Set(chain.met_px, chain.met_py)

  ##HME uses two leptons, two jets, MET as inputs
  hme = HeavyMassEstimator()
  hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
  #hme.showKinematic()
  hme.setIterations(iterations)
  ##hme.setMETCovMatrix(50.0, 40.0, 30.0, True) ## example to use MET cov smearing. met_covxx/covyy/covxy taken from events
  if metRes != None:
      hme.setMETResolution(metRes)
  hme.runHME()

  if hme.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      hme_mass_peak[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
      hme_mass_peak_divSol[0] = hme.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme.hme_h2Mass_divSolutions.GetMaximumBin())
      hme.hme_h2Mass.Rebin(nbins_rebin) 
      hme_mass_peak_gravity[0] = get_gravity_center(hme.hme_h2Mass, per_weight_gravitycenter)
      nsolution0[0] = hme.hme_nsolutions.GetBinContent(1)
      nsolution1[0] = hme.hme_nsolutions.GetBinContent(2)
      nsolution2[0] = hme.hme_nsolutions.GetBinContent(3)
      nsolution3[0] = hme.hme_nsolutions.GetBinContent(4)
      nsolution4[0] = hme.hme_nsolutions.GetBinContent(5)

  binwidth = hme.hme_h2Mass.GetBinWidth(1)
  print("ievent ", nEv," hme gave max mass = ", hme_mass_peak[0]," nsolution2[0] ",nsolution2[0]," weighted center  ", hme_mass_peak_gravity[0]," binwidth ",binwidth)

  ##HME uses two leptons, dijet, MET as inputs, later this is proven to be less efficient
  #dijet_p4.SetPxPyPzE(chain.fatbjet_Px,   chain.fatbjet_Py,    chain.fatbjet_Pz,    chain.fatbjet_E)
  #hme_boosted = HeavyMassEstimator()
  #hme_boosted.hme_h2Mass.Rebin(nbins_rebin)
  #hme_boosted.setKinematic_boosted(lep1_p4, lep2_p4, dijet_p4, met_vec2, 0.0)
  #hme_boosted.showKinematic()
  #hme_boosted.setIterations(iterations)
  #if metRes != None:
  #    hme_boosted.setMETResolution(metRes)
  #hme_boosted.runHME()
  #hme_mass_peak_boosted[0] = hme_boosted.hme_h2Mass.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass.GetMaximumBin())
  #hme_mass_peak_boosted_divSol[0] = hme_boosted.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass_divSolutions.GetMaximumBin())
  #print("hme boosted gave max mass = ", hme_mass_peak_boosted[0])

  tree.Fill()

tree.Write()
out.Close()
print("========== HME is Done! ==========")
