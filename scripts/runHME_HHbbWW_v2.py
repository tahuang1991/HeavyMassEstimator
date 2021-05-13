import ROOT
from HeavyMassEstimator import *
from array import array
import os
import sys

print sys.argv
for i in range(len(sys.argv)):
  if sys.argv[i] == "-i":
    infile_name = sys.argv[i+1]
  if sys.argv[i] == "-o":
    outfile_name = sys.argv[i+1]

if ("-i" or "-o") not in sys.argv:
  print "No input/output file given"
  sys.exit("Use 'python runHME_HHbbWW_v2.py -i <inputFile> -o <outputFile>'")

#f = ROOT.TFile("testfiles/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2_HME_Friends.root")
f = ROOT.TFile(infile_name)
chain = f.Get("Friends")

os.system("mkdir -p testfiles/output")

#out = ROOT.TFile("testfiles/output/m260.root", 'recreate')
out = ROOT.TFile(outfile_name, 'recreate')
out.cd()
tree = chain.CloneTree(10)
out.Write()

f.Close()

iterations = 10000
TotalEv = tree.GetEntries()
nStart = 1
nEnd = TotalEv + 1
print "Total events = ", TotalEv

hme_mass_peak = array('f', [0.])
hme_mass_peak_divSol = array('f', [0.])
hme_mass_peak_boosted = array('f', [0.])
hme_mass_peak_boosted_divSol = array('f', [0.])

massBranch = tree.Branch("hme_mass_peak", hme_mass_peak, "hme_mass_peak/F")
massBranch_divSol = tree.Branch("hme_mass_peak_divSol", hme_mass_peak_divSol, "hme_mass_peak_divSol/F")
massBranch_boosted = tree.Branch("hme_mass_peak_boosted", hme_mass_peak_boosted, "hme_mass_peak_boosted/F")
massBranch_boosted_divSol = tree.Branch("hme_mass_peak_boosted_divSol", hme_mass_peak_boosted_divSol, "hme_mass_peak_boosted_divSol/F")

for nEv in range(nStart, nEnd):
  tree.GetEntry(nEv)
  lep1_p4 = ROOT.TLorentzVector()
  lep2_p4 = ROOT.TLorentzVector()
  jet1_p4 = ROOT.TLorentzVector()
  jet2_p4 = ROOT.TLorentzVector()
  dijet_p4 = ROOT.TLorentzVector()
  met_vec2 = ROOT.TVector2()

  lep1_p4.SetPtEtaPhiM(tree.lep1_pt, tree.lep1_eta, tree.lep1_phi, 0)
  lep2_p4.SetPtEtaPhiM(tree.lep2_pt, tree.lep2_eta, tree.lep2_phi, 0)
  jet1_p4.SetPtEtaPhiE(tree.jet1_pt, tree.jet1_eta, tree.jet1_phi, tree.jet1_E)
  jet2_p4.SetPtEtaPhiE(tree.jet2_pt, tree.jet2_eta, tree.jet2_phi, tree.jet2_E)
  dijet_p4 = jet1_p4 + jet2_p4
  met_vec2.SetMagPhi(tree.MET_pt_nom, tree.MET_phi_nom)


  hme = HeavyMassEstimator()
  hme_boosted = HeavyMassEstimator()
  hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
  hme.showKinematic()
  hme.setIterations(iterations)
  hme.runHME()
  hme_mass_peak[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
  hme_mass_peak_divSol[0] = hme.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme.hme_h2Mass_divSolutions.GetMaximumBin())
  #print "hme gave max mass = ", hme_mass_peak[0]

  hme_boosted = HeavyMassEstimator()
  hme_boosted.setKinematic_boosted(lep1_p4, lep2_p4, dijet_p4, met_vec2, 0.0)
  hme_boosted.showKinematic()
  hme_boosted.setIterations(iterations)
  hme_boosted.runHME()
  hme_mass_peak_boosted[0] = hme_boosted.hme_h2Mass.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass.GetMaximumBin())
  hme_mass_peak_boosted_divSol[0] = hme_boosted.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass_divSolutions.GetMaximumBin())
  #print "hme boosted gave max mass = ", hme_boosted_maxMass

  massBranch.Fill()
  massBranch_divSol.Fill()
  massBranch_boosted.Fill()
  massBranch_boosted_divSol.Fill()

out.Write()


