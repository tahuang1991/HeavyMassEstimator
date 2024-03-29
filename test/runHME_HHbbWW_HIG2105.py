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

sys.path.append(os.path.join(os.path.abspath('..'), 'python'))
from HeavyMassEstimator import *

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-i", "--inputFile", dest="infile", type=str, default=None, help="input file name. [Default: None]")
parser.add_argument("-o", "--outputFile", dest="outfile", type=str, default="out.root", help="output file name. [Default: out.root]")
parser.add_argument("-it", "--iterations", dest="iterations", type=int, default=10000, help="iteration number used HME [Default: 100000]")
parser.add_argument("-nStart", "--nStart", dest="nStart", type=int, default=0, help="first event index to process. [Default: 0]")
parser.add_argument("-nEnd", "--nEnd", dest="nEnd", type=int, default=-1, help="last event index to process. -1 means last event in TTree [Default: -1]")
#args = parser.parse_args()
args, unknown = parser.parse_known_args()

if args.infile == None or not(os.path.isfile(args.infile)):
  print("No input file given")
  sys.exit("Use 'python runHME_HHbbWW_v2.py -i <inputFile> -o <outputFile>'")

#chain = ROOT.TChain("Events")
chain = ROOT.TChain("Double_Tree")
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

hme_mass_peak_metunclusterup = array('f', [0.])
hme_mass_peak_metunclusterup_divSol = array('f', [0.])
hme_mass_peak_metunclusterdown = array('f', [0.])
hme_mass_peak_metunclusterdown_divSol = array('f', [0.])
hme_mass_peak_JER_up = array('f', [0.])
hme_mass_peak_JER_up_divSol = array('f', [0.])
hme_mass_peak_JER_down = array('f', [0.])
hme_mass_peak_JER_down_divSol = array('f', [0.])
hme_mass_peak_JES_up = array('f', [0.])
hme_mass_peak_JES_up_divSol = array('f', [0.])
hme_mass_peak_JES_down = array('f', [0.])
hme_mass_peak_JES_down_divSol = array('f', [0.])

def initBranches():
  hme_mass_peak[0] = 0.0
  hme_mass_peak_divSol[0] = 0.0
  hme_mass_peak_metunclusterup[0] = 0.0
  hme_mass_peak_metunclusterup_divSol[0] = 0.0
  hme_mass_peak_metunclusterdown[0] = 0.0
  hme_mass_peak_metunclusterdown_divSol[0] = 0.0
  hme_mass_peak_JER_up[0] = 0.0
  hme_mass_peak_JER_up_divSol[0] = 0.0
  hme_mass_peak_JER_down[0] = 0.0
  hme_mass_peak_JER_down_divSol[0] = 0.0
  hme_mass_peak_JES_up[0] = 0.0
  hme_mass_peak_JES_up_divSol[0] = 0.0
  hme_mass_peak_JES_down[0] = 0.0
  hme_mass_peak_JES_down_divSol[0] = 0.0

massBranch                = tree.Branch("hme_mass_peak", hme_mass_peak, "hme_mass_peak/F")
massBranch_divSol         = tree.Branch("hme_mass_peak_divSol", hme_mass_peak_divSol, "hme_mass_peak_divSol/F")
mass_JER_up_Branch           = tree.Branch("hme_mass_peak_JER_up",          hme_mass_peak_JER_up,          "hme_mass_peak_JER_up/F")
mass_JER_up_Branch_divSol    = tree.Branch("hme_mass_peak_JER_up_divSol",   hme_mass_peak_JER_up_divSol,   "hme_mass_peak_JER_up_divSol/F")
mass_JER_down_Branch         = tree.Branch("hme_mass_peak_JER_down",        hme_mass_peak_JER_down,        "hme_mass_peak_JER_down/F")
mass_JER_down_Branch_divSol  = tree.Branch("hme_mass_peak_JER_down_divSol", hme_mass_peak_JER_down_divSol, "hme_mass_peak_JER_down_divSol/F")
mass_JES_up_Branch           = tree.Branch("hme_mass_peak_JES_up",          hme_mass_peak_JES_up,          "hme_mass_peak_JES_up/F")
mass_JES_up_Branch_divSol    = tree.Branch("hme_mass_peak_JES_up_divSol",   hme_mass_peak_JES_up_divSol,   "hme_mass_peak_JES_up_divSol/F")
mass_JES_down_Branch         = tree.Branch("hme_mass_peak_JES_down",        hme_mass_peak_JES_down,        "hme_mass_peak_JES_down/F")
mass_JES_down_Branch_divSol  = tree.Branch("hme_mass_peak_JES_down_divSol", hme_mass_peak_JES_down_divSol, "hme_mass_peak_JES_down_divSol/F")

for nEv in range(nStart, nEnd):
  initBranches()
  chain.GetEntry(nEv)
  lep1_p4 = ROOT.TLorentzVector()
  lep2_p4 = ROOT.TLorentzVector()
  jet1_p4 = ROOT.TLorentzVector()
  jet2_p4 = ROOT.TLorentzVector()
  jet1_p4_JER_up = ROOT.TLorentzVector()
  jet2_p4_JER_up = ROOT.TLorentzVector()
  jet1_p4_JER_down = ROOT.TLorentzVector()
  jet2_p4_JER_down = ROOT.TLorentzVector()
  jet1_p4_JES_up = ROOT.TLorentzVector()
  jet2_p4_JES_up = ROOT.TLorentzVector()
  jet1_p4_JES_down = ROOT.TLorentzVector()
  jet2_p4_JES_down = ROOT.TLorentzVector()
  dijet_p4 = ROOT.TLorentzVector()
  met_vec2 = ROOT.TVector2()
  met_vec2_up = ROOT.TVector2()
  met_vec2_down = ROOT.TVector2()

  #goodevent = (chain.double_is_ee or chain.double_is_mm or chain.double_is_em) and (chain.Double_Res_1b or chain.Double_Res_2b or chain.Double_HbbFat)
  goodevent = (chain.double_is_ee or chain.double_is_mm or chain.double_is_em) and (chain.Double_Res_1b or chain.Double_Res_2b)
  if not goodevent:
      tree.Fill() 
      continue


  #lep1_p4.SetPtEtaPhiM(chain.lep1_pt, chain.lep1_eta, chain.lep1_phi, 0)
  #lep2_p4.SetPtEtaPhiM(chain.lep2_pt, chain.lep2_eta, chain.lep2_phi, 0)
  #jet1_p4.SetPtEtaPhiE(chain.jet1_pt, chain.jet1_eta, chain.jet1_phi, chain.jet1_E)
  #jet2_p4.SetPtEtaPhiE(chain.jet2_pt, chain.jet2_eta, chain.jet2_phi, chain.jet2_E)
  #met_vec2.SetMagPhi(chain.MET_pt_nom, chain.MET_phi_nom)
  lep1_p4.SetXYZT(chain.lep0_px, chain.lep0_py, chain.lep0_pz, chain.lep0_E)
  lep2_p4.SetXYZT(chain.lep1_px, chain.lep1_py, chain.lep1_pz, chain.lep1_E)
  #jet1_p4.SetXYZT(chain.ak4_jet0_px, chain.ak4_jet0_py, chain.ak4_jet0_pz, chain.ak4_jet0_E)
  #jet2_p4.SetXYZT(chain.ak4_jet1_px, chain.ak4_jet1_py, chain.ak4_jet1_pz, chain.ak4_jet1_E)
  jet1_p4.SetPtEtaPhiE(chain.ak4_jet0_pt_original, chain.ak4_jet0_eta, chain.ak4_jet0_phi, chain.ak4_jet0_E)
  jet2_p4.SetPtEtaPhiE(chain.ak4_jet1_pt_original, chain.ak4_jet1_eta, chain.ak4_jet1_phi, chain.ak4_jet1_E)
  jet1_p4_JER_up.SetPtEtaPhiE(  chain.ak4_jet0_pt_JER_up, chain.ak4_jet0_eta, chain.ak4_jet0_phi, chain.ak4_jet0_E)
  jet2_p4_JER_up.SetPtEtaPhiE(  chain.ak4_jet1_pt_JER_up, chain.ak4_jet1_eta, chain.ak4_jet1_phi, chain.ak4_jet1_E)
  jet1_p4_JER_down.SetPtEtaPhiE(chain.ak4_jet0_pt_JER_down, chain.ak4_jet0_eta, chain.ak4_jet0_phi, chain.ak4_jet0_E)
  jet2_p4_JER_down.SetPtEtaPhiE(chain.ak4_jet1_pt_JER_down, chain.ak4_jet1_eta, chain.ak4_jet1_phi, chain.ak4_jet1_E)
  jet1_p4_JES_up.SetPtEtaPhiE(  chain.ak4_jet0_pt_JES_up, chain.ak4_jet0_eta, chain.ak4_jet0_phi, chain.ak4_jet0_E)
  jet2_p4_JES_up.SetPtEtaPhiE(  chain.ak4_jet1_pt_JES_up, chain.ak4_jet1_eta, chain.ak4_jet1_phi, chain.ak4_jet1_E)
  jet1_p4_JES_down.SetPtEtaPhiE(chain.ak4_jet0_pt_JES_down, chain.ak4_jet0_eta, chain.ak4_jet0_phi, chain.ak4_jet0_E)
  jet2_p4_JES_down.SetPtEtaPhiE(chain.ak4_jet1_pt_JES_down, chain.ak4_jet1_eta, chain.ak4_jet1_phi, chain.ak4_jet1_E)
  if chain.Double_HbbFat:
    jet1_p4.SetXYZT(chain.ak8_jet0_subjet0_px, chain.ak8_jet0_subjet0_py, chain.ak8_jet0_subjet0_pz, chain.ak8_jet0_subjet0_E)
    jet2_p4.SetXYZT(chain.ak8_jet0_subjet1_px, chain.ak8_jet0_subjet1_py, chain.ak8_jet0_subjet1_pz, chain.ak8_jet0_subjet1_E)

  met_vec2.SetX(chain.met_px)
  met_vec2.SetY(chain.met_py)
  met_vec2_up.SetX(chain.met_px+chain.met_unclust_energy_up_x)
  met_vec2_up.SetY(chain.met_py+chain.met_unclust_energy_up_y)
  met_vec2_down.SetX(chain.met_px-chain.met_unclust_energy_up_x)
  met_vec2_down.SetY(chain.met_py-chain.met_unclust_energy_up_y)

  hme = HeavyMassEstimator()
  hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
  #hme.showKinematic()
  hme.setIterations(iterations)
  hme.runHME()
  if hme.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      hme_mass_peak[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
      hme_mass_peak_divSol[0] = hme.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme.hme_h2Mass_divSolutions.GetMaximumBin())


  hme_metunclusterup = HeavyMassEstimator()
  hme_metunclusterup.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2_up, 0)
  hme_metunclusterup.showKinematic()
  hme_metunclusterup.setIterations(iterations)
  hme_metunclusterup.runHME()
  if hme_metunclusterup.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      maxbin1 = hme_metunclusterup.hme_h2Mass.GetMaximumBin()
      hme_mass_peak_metunclusterup[0]        = hme_metunclusterup.hme_h2Mass.GetXaxis().GetBinCenter(maxbin1)
      hme_mass_peak_metunclusterup_divSol[0] = hme_metunclusterup.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(maxbin1)

  hme_metunclusterdown = HeavyMassEstimator()
  hme_metunclusterdown.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2_down, 0)
  hme_metunclusterdown.showKinematic()
  hme_metunclusterdown.setIterations(iterations)
  hme_metunclusterdown.runHME()
  if hme_metunclusterdown.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      maxbin2 = hme_metunclusterdown.hme_h2Mass.GetMaximumBin()
      hme_mass_peak_metunclusterdown[0]        = hme_metunclusterdown.hme_h2Mass.GetXaxis().GetBinCenter(maxbin2)
      hme_mass_peak_metunclusterdown_divSol[0] = hme_metunclusterdown.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(maxbin2)

  print("ievent ", nEv," hme gave max mass = ", hme_mass_peak[0], " metuncluset up ", hme_mass_peak_metunclusterup[0], " down ", hme_mass_peak_metunclusterdown[0])

  hme_JER_up = HeavyMassEstimator()
  hme_JER_up.setKinematic(lep1_p4, lep2_p4, jet1_p4_JER_up, jet2_p4_JER_up, met_vec2, 0)
  hme_JER_up.showKinematic()
  hme_JER_up.setIterations(iterations)
  hme_JER_up.runHME()
  if hme_JER_up.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      maxbin3                                = hme_JER_up.hme_h2Mass.GetMaximumBin()
      hme_mass_peak_JER_up[0]        = hme_JER_up.hme_h2Mass.GetXaxis().GetBinCenter(maxbin3)
      hme_mass_peak_JER_up_divSol[0] = hme_JER_up.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(maxbin3)

  hme_JER_down = HeavyMassEstimator()
  hme_JER_down.setKinematic(lep1_p4, lep2_p4, jet1_p4_JER_down, jet2_p4_JER_down, met_vec2, 0)
  hme_JER_down.showKinematic()
  hme_JER_down.setIterations(iterations)
  hme_JER_down.runHME()
  if hme_JER_down.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      maxbin4                                = hme_JER_down.hme_h2Mass.GetMaximumBin()
      hme_mass_peak_JER_down[0]              = hme_JER_down.hme_h2Mass.GetXaxis().GetBinCenter(maxbin4)
      hme_mass_peak_JER_down_divSol[0]       = hme_JER_down.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(maxbin4)


  hme_JES_up = HeavyMassEstimator()
  hme_JES_up.setKinematic(lep1_p4, lep2_p4, jet1_p4_JES_up, jet2_p4_JES_up, met_vec2, 0)
  hme_JES_up.showKinematic()
  hme_JES_up.setIterations(iterations)
  hme_JES_up.runHME()
  if hme_JES_up.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      maxbin5                        = hme_JES_up.hme_h2Mass.GetMaximumBin()
      hme_mass_peak_JES_up[0]        = hme_JES_up.hme_h2Mass.GetXaxis().GetBinCenter(maxbin5)
      hme_mass_peak_JES_up_divSol[0] = hme_JES_up.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(maxbin5)

  hme_JES_down = HeavyMassEstimator()
  hme_JES_down.setKinematic(lep1_p4, lep2_p4, jet1_p4_JES_down, jet2_p4_JES_down, met_vec2, 0)
  hme_JES_down.showKinematic()
  hme_JES_down.setIterations(iterations)
  hme_JES_down.runHME()
  if hme_JES_down.hme_h2Mass.GetEntries() >0:
      ##fill TTree if it return valid HME value, otherwise fill 0.0
      maxbin6                                = hme_JES_down.hme_h2Mass.GetMaximumBin()
      hme_mass_peak_JES_down[0]              = hme_JES_down.hme_h2Mass.GetXaxis().GetBinCenter(maxbin6)
      hme_mass_peak_JES_down_divSol[0]       = hme_JES_down.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(maxbin6)


  tree.Fill() 

tree.Write()
out.Close()
print("========== HME is Done! ==========")


