import ROOT
from HeavyMassEstimator import *
from array import array
import os
import sys

print sys.argv
nbins_rebin = 1
iterations = 10000
for i in range(len(sys.argv)):
  if sys.argv[i] == "-i":
    infile_name = sys.argv[i+1]
  if sys.argv[i] == "-o":
    outfile_name = sys.argv[i+1]
  if sys.argv[i] == "--rebin":
    nbins_rebin = int(sys.argv[i+1])
  if sys.argv[i] == "--iteration":
    iterations  = int(sys.argv[i+1])

if ("-i" or "-o") not in sys.argv:
  print "No input/output file given"
  sys.exit("Use 'python runHME_HHbbWW_v2.py -i <inputFile> -o <outputFile>'")

#f = ROOT.TFile("testfiles/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2_HME_Friends.root")
f = ROOT.TFile(infile_name)
#chain = f.Get("Friends")
chain = f.Get("Events")

#os.system("mkdir -p output/")

#out = ROOT.TFile("testfiles/output/m260.root", 'recreate')
#out = ROOT.TFile("output/"+outfile_name, 'recreate')
out = ROOT.TFile(outfile_name, 'recreate')
out.cd()
#tree = chain.CloneTree(10)
tree = chain.CloneTree(-1)
out.Write()

f.Close()

TotalEv = tree.GetEntries()
nStart = 0
nEnd = TotalEv 
print "Total events = ", TotalEv

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



for nEv in range(nStart, nEnd):
  tree.GetEntry(nEv)
  print ("nEv ", nEv, " l1px ",tree.l1_Px, " py ",tree.l1_Py," pz ",tree.l1_Pz, " E ",tree.l1_E )
  #continue
  lep1_p4 = ROOT.TLorentzVector()
  lep2_p4 = ROOT.TLorentzVector()
  jet1_p4 = ROOT.TLorentzVector()
  jet2_p4 = ROOT.TLorentzVector()
  dijet_p4 = ROOT.TLorentzVector()
  met_vec2 = ROOT.TVector2()

  lep1_p4.SetPxPyPzE(tree.l1_Px,   tree.l1_Py,    tree.l1_Pz,    tree.l1_E)
  lep2_p4.SetPxPyPzE(tree.l2_Px,   tree.l2_Py,    tree.l2_Pz,    tree.l2_E)
  jet1_p4.SetPxPyPzE(tree.fatbjet_subjet1_Px,   tree.fatbjet_subjet1_Py,    tree.fatbjet_subjet1_Pz,    tree.fatbjet_subjet1_E)
  jet2_p4.SetPxPyPzE(tree.fatbjet_subjet2_Px,   tree.fatbjet_subjet2_Py,    tree.fatbjet_subjet2_Pz,    tree.fatbjet_subjet2_E)
  dijet_p4.SetPxPyPzE(tree.fatbjet_Px,   tree.fatbjet_Py,    tree.fatbjet_Pz,    tree.fatbjet_E)
  #dijet_p4 = jet1_p4 + jet2_p4
  #met_vec2.SetMagPhi(tree.MET_pt_nom, tree.MET_phi_nom)
  met_vec2.Set(tree.met_Px, tree.met_Py)


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
  print "hme gave max mass = ", hme_mass_peak[0]," bincontent 3 ",hme.hme_nsolutions.GetBinContent(3)," nsolution2[0] ",nsolution2[0]
  #hme.hme_nsolutions.Print("ALL")

  hme_boosted = HeavyMassEstimator()
  hme_boosted.setKinematic_boosted(lep1_p4, lep2_p4, dijet_p4, met_vec2, 0.0)
  hme_boosted.showKinematic()
  hme_boosted.setIterations(iterations)
  hme_boosted.runHME()
  hme_boosted.hme_h2Mass.Rebin(nbins_rebin)
  hme_mass_peak_boosted[0] = hme_boosted.hme_h2Mass.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass.GetMaximumBin())
  hme_mass_peak_boosted_divSol[0] = hme_boosted.hme_h2Mass_divSolutions.GetXaxis().GetBinCenter(hme_boosted.hme_h2Mass_divSolutions.GetMaximumBin())
  print "hme boosted gave max mass = ", hme_mass_peak_boosted[0]

  massBranch.Fill()
  massBranch_divSol.Fill()
  massBranch_boosted.Fill()
  massBranch_boosted_divSol.Fill()
  br_nsolution0.Fill()
  br_nsolution1.Fill()
  br_nsolution2.Fill()
  br_nsolution3.Fill()
  br_nsolution4.Fill()


out.Write()


