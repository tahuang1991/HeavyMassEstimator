import ROOT
from array import array
import argparse
import os
import sys

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
tree = chain.CloneTree(10)


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
nsolution0 = array('f', [0.0])
nsolution1 = array('f', [0.0])
nsolution2 = array('f', [0.0])
nsolution3 = array('f', [0.0])
nsolution4 = array('f', [0.0])

massBranch = tree.Branch("hme_mass_peak", hme_mass_peak, "hme_mass_peak/F")
massBranch_divSol = tree.Branch("hme_mass_peak_divSol", hme_mass_peak_divSol, "hme_mass_peak_divSol/F")
massBranch_boosted = tree.Branch("hme_mass_peak_boosted", hme_mass_peak_boosted, "hme_mass_peak_boosted/F")
massBranch_boosted_divSol = tree.Branch("hme_mass_peak_boosted_divSol", hme_mass_peak_boosted_divSol, "hme_mass_peak_boosted_divSol/F")
br_nsolution0 = tree.Branch("hme_nsolution0", nsolution0, "nsolution0/F")
br_nsolution1 = tree.Branch("hme_nsolution1", nsolution1, "nsolution1/F")
br_nsolution2 = tree.Branch("hme_nsolution2", nsolution2, "nsolution2/F")
br_nsolution3 = tree.Branch("hme_nsolution3", nsolution3, "nsolution3/F")
br_nsolution4 = tree.Branch("hme_nsolution4", nsolution4, "nsolution4/F")

def initBranches():
  hme_mass_peak[0] = 0.0
  hme_mass_peak_divSol[0] = 0.0
  hme_mass_peak_boosted[0] = 0.0
  hme_mass_peak_boosted_divSol[0] = 0.0
  nsolution0[0] = 0.0
  nsolution1[0] = 0.0
  nsolution2[0] = 0.0
  nsolution3[0] = 0.0
  nsolution4[0] = 0.0

nbins_rebin = 100
for nEv in range(nStart, nEnd):
  initBranches()
  chain.GetEntry(nEv)
  #print ("nEv ", nEv, " l1px ",chain.l1_Px, " py ",chain.l1_Py," pz ",chain.l1_Pz, " E ",chain.l1_E )
  lep1_p4 = ROOT.TLorentzVector()
  lep2_p4 = ROOT.TLorentzVector()
  jet1_p4 = ROOT.TLorentzVector()
  jet2_p4 = ROOT.TLorentzVector()
  dijet_p4 = ROOT.TLorentzVector()
  met_vec2 = ROOT.TVector2()

  lep1_p4.SetPxPyPzE(chain.l1_Px,   chain.l1_Py,    chain.l1_Pz,    chain.l1_E)
  lep2_p4.SetPxPyPzE(chain.l2_Px,   chain.l2_Py,    chain.l2_Pz,    chain.l2_E)
  jet1_p4.SetPxPyPzE(chain.fatbjet_subjet1_Px,   chain.fatbjet_subjet1_Py,    chain.fatbjet_subjet1_Pz,    chain.fatbjet_subjet1_E)
  jet2_p4.SetPxPyPzE(chain.fatbjet_subjet2_Px,   chain.fatbjet_subjet2_Py,    chain.fatbjet_subjet2_Pz,    chain.fatbjet_subjet2_E)
  dijet_p4.SetPxPyPzE(chain.fatbjet_Px,   chain.fatbjet_Py,    chain.fatbjet_Pz,    chain.fatbjet_E)
  #dijet_p4 = jet1_p4 + jet2_p4
  #met_vec2.SetMagPhi(chain.MET_pt_nom, chain.MET_phi_nom)
  met_vec2.Set(chain.met_Px, chain.met_Py)


  hme = HeavyMassEstimator()
  hme_boosted = HeavyMassEstimator()
  hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
  hme.showKinematic()
  hme.setIterations(iterations)
  hme.runHME()
  hme.hme_h2Mass.Rebin(nbins_rebin)
  hme_mass_peak[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
  hme_mass_peak_divSol[0] = hme.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme.hme_h2Mass_divSolutions.GetMaximumBin())
  nsolution0[0] = hme.hme_nsolutions.GetBinContent(1)
  nsolution1[0] = hme.hme_nsolutions.GetBinContent(2)
  nsolution2[0] = hme.hme_nsolutions.GetBinContent(3)
  nsolution3[0] = hme.hme_nsolutions.GetBinContent(4)
  nsolution4[0] = hme.hme_nsolutions.GetBinContent(5)
  print("hme gave max mass = ", hme_mass_peak[0]," bincontent 3 ",hme.hme_nsolutions.GetBinContent(3)," nsolution2[0] ",nsolution2[0])
  #hme.hme_nsolutions.Print("ALL")

  hme_boosted = HeavyMassEstimator()
  hme_boosted.setKinematic_boosted(lep1_p4, lep2_p4, dijet_p4, met_vec2, 0.0)
  hme_boosted.showKinematic()
  hme_boosted.setIterations(iterations)
  hme_boosted.runHME()
  hme_boosted.hme_h2Mass.Rebin(nbins_rebin)
  hme_mass_peak_boosted[0] = hme_boosted.hme_h2Mass.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass.GetMaximumBin())
  hme_mass_peak_boosted_divSol[0] = hme_boosted.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass_divSolutions.GetMaximumBin())
  print("hme boosted gave max mass = ", hme_mass_peak_boosted[0])

  #massBranch.Fill()
  #massBranch_divSol.Fill()
  #massBranch_boosted.Fill()
  #massBranch_boosted_divSol.Fill()
  #br_nsolution0.Fill()
  #br_nsolution1.Fill()
  #br_nsolution2.Fill()
  #br_nsolution3.Fill()
  #br_nsolution4.Fill()
  tree.Fill()

tree.Write()
out.Close()
print("========== HME is Done! ==========")