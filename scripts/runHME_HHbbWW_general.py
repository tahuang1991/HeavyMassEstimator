import string
import ROOT
from array import array
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.argv.append( '-b' )
sys.argv.append( '-q' )

from HeavyMassEstimator import *

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-i", "--inputFile", dest="infile", type=str, default=None, help="input file name. [Default: None]")
parser.add_argument("-o", "--outputFile", dest="outfile", type=str, default="out.root", help="output file name. [Default: out.root]")
parser.add_argument("-it", "--iterations", dest="iterations", type=int, default=10000, help="iteration number used HME [Default: 100000]")
parser.add_argument("-nStart", "--nStart", dest="nStart", type=int, default=0, help="first event index to process. [Default: 0]")
parser.add_argument("-nEnd", "--nEnd", dest="nEnd", type=int, default=-1, help="last event index to process. -1 means last event in TTree [Default: -1]")
#args = parser.parse_args()
args, unknown = parser.parse_known_args()

if args.infile == None:
  print("No input file given")
  sys.exit("Use 'python runHME_HHbbWW_v2.py -i <inputFile> -o <outputFile>'")

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
print("Total events = ", TotalEv, " nStart ", nStart, " nEnd ", nEnd, " iteration for HME ", iterations)

hme_mass_peak = array('f', [0.])
hme_mass_peak_divSol = array('f', [0.])
hme_mass_peak_boosted = array('f', [0.])
hme_mass_peak_boosted_divSol = array('f', [0.])

def initBranches():
  hme_mass_peak[0] = 0.0
  hme_mass_peak_divSol[0] = 0.0
  hme_mass_peak_boosted[0] = 0.0
  hme_mass_peak_boosted_divSol[0] = 0.0

massBranch = tree.Branch("hme_mass_peak", hme_mass_peak, "hme_mass_peak/F")
massBranch_divSol = tree.Branch("hme_mass_peak_divSol", hme_mass_peak_divSol, "hme_mass_peak_divSol/F")
#massBranch_boosted = tree.Branch("hme_mass_peak_boosted", hme_mass_peak_boosted, "hme_mass_peak_boosted/F")
#massBranch_boosted_divSol = tree.Branch("hme_mass_peak_boosted_divSol", hme_mass_peak_boosted_divSol, "hme_mass_peak_boosted_divSol/F")

for nEv in range(nStart, nEnd):
  initBranches()
  chain.GetEntry(nEv)
  lep1_p4 = ROOT.TLorentzVector()
  lep2_p4 = ROOT.TLorentzVector()
  jet1_p4 = ROOT.TLorentzVector()
  jet2_p4 = ROOT.TLorentzVector()
  dijet_p4 = ROOT.TLorentzVector()
  met_vec2 = ROOT.TVector2()

  #lep1_p4.SetPtEtaPhiM(chain.lep1_pt, chain.lep1_eta, chain.lep1_phi, 0)
  #lep2_p4.SetPtEtaPhiM(chain.lep2_pt, chain.lep2_eta, chain.lep2_phi, 0)
  #jet1_p4.SetPtEtaPhiE(chain.jet1_pt, chain.jet1_eta, chain.jet1_phi, chain.jet1_E)
  #jet2_p4.SetPtEtaPhiE(chain.jet2_pt, chain.jet2_eta, chain.jet2_phi, chain.jet2_E)
  #met_vec2.SetMagPhi(chain.MET_pt_nom, chain.MET_phi_nom)
  lep1_p4.SetXYZT(chain.l1_px, chain.l1_py, chain.l1_pz, chain.l1_E)
  lep2_p4.SetXYZT(chain.l2_px, chain.l2_py, chain.l2_pz, chain.l2_E)
  jet1_p4.SetXYZT(chain.j1_px, chain.j1_py, chain.j1_pz, chain.j1_E)
  jet2_p4.SetXYZT(chain.j2_px, chain.j2_py, chain.j2_pz, chain.j2_E)
  met_vec2.SetX(chain.met_px)
  met_vec2.SetY(chain.met_py)

  hme = HeavyMassEstimator()
  hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
  hme.showKinematic()
  hme.setIterations(iterations)
  hme.runHME()
  hme_mass_peak[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
  hme_mass_peak_divSol[0] = hme.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme.hme_h2Mass_divSolutions.GetMaximumBin())
  print("ievent ", nEv," hme gave max mass = ", hme_mass_peak[0])

  ## HME inputs with two leptons, dijets, met
  #dijet_p4 = jet1_p4 + jet2_p4
  #hme_boosted = HeavyMassEstimator()
  #hme_boosted.setKinematic_boosted(lep1_p4, lep2_p4, dijet_p4, met_vec2, 0.0)
  #hme_boosted.showKinematic()
  #hme_boosted.setIterations(iterations)
  #hme_boosted.runHME()
  #hme_mass_peak_boosted[0] = hme_boosted.hme_h2Mass.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass.GetMaximumBin())
  #hme_mass_peak_boosted_divSol[0] = hme_boosted.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass_divSolutions.GetMaximumBin())
  #print("hme boosted(with dijet as input) gave max mass = ", hme_boosted_maxMass)

  #massBranch.Fill()
  #massBranch_divSol.Fill()
  #massBranch_boosted.Fill()
  #massBranch_boosted_divSol.Fill()
  tree.Fill() 

tree.Write()
out.Close()
print("========== HME is Done! ==========")


