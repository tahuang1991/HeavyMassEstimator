import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
from math import *
import argparse
import numpy as np
from array import array

import warnings
sys.argv.append( '-b' )
sys.argv.append( '-q' )

warnings.filterwarnings("ignore", category=RuntimeWarning)

from HeavyMassEstimator import *

doTest = False
doHME = True
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

tree_name="treeMaker/Events"
TCha = ROOT.TChain(tree_name)
TCha.Add(args.infile)
TotalEv = TCha.GetEntries()

f = ROOT.TFile(args.outfile, 'recreate'); f.cd()
TCha2 = TCha.CloneTree(0)

nStart = args.nStart
nEnd = TotalEv
if args.nEnd > 0 and args.nEnd <= TotalEv:
  nEnd = args.nEnd
if nStart >= nEnd:
  sys.exit("Error! first event index to process is larger than last one: nStart ", nStart, " nEnd ", nEnd)
print("Total events = ", TotalEv, " nStart ", nStart, " nEnd ", nEnd, " iteration for HME ", iterations)

### add HME information to TCha2
maxn = 1 
ak8jetindex                 = array( 'i', maxn*[ 0 ] ) #np.zeros(1, dtype=float)
lep1index                   = array( 'i', maxn*[ 0 ] ) #np.zeros(1, dtype=float)
lep2index                   = array( 'i', maxn*[ 0 ] ) #np.zeros(1, dtype=float)
leptonpairtype              = array( 'i', maxn*[ 0 ] ) #np.zeros(1, dtype=float)
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

TCha2.Branch("ak8jetindex",                   ak8jetindex,                  "ak8jetindex/I")
TCha2.Branch("lep1index",                     lep1index,                    "lep1index/I")
TCha2.Branch("lep2index",                     lep2index,                    "lep2index/I")
TCha2.Branch("leptonpairtype",                leptonpairtype,               "leptonpairtype/I")
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

hme_min = 200.0; hme_max = 4000; hme_nbin = 3800
h_h2mass_weight_gen      = ROOT.TH1F("h_h2mass_weight_gen","", hme_nbin, hme_min, hme_max);  h_h2mass_weight_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_gen_sum  = ROOT.TH1F("h_h2mass_weight_gen_sum","", hme_nbin, hme_min, hme_max);  h_h2mass_weight_gen_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_gen     = ROOT.TH1F("h_h2mass_weight1_gen","", hme_nbin, hme_min, hme_max);  h_h2mass_weight1_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight2_gen     = ROOT.TH1F("h_h2mass_weight2_gen","", hme_nbin, hme_min, hme_max);  h_h2mass_weight2_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco     = ROOT.TH1F("h_h2mass_weight_reco","", hme_nbin, hme_min, hme_max);  h_h2mass_weight_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco_sum = ROOT.TH1F("h_h2mass_weight_reco_sum","",hme_nbin, hme_min, hme_max);  h_h2mass_weight_reco_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_reco    = ROOT.TH1F("h_h2mass_weight1_reco","",hme_nbin, hme_min, hme_max);  h_h2mass_weight1_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight2_reco     = ROOT.TH1F("h_h2mass_weight2_reco","",hme_nbin, hme_min, hme_max);  h_h2mass_weight2_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
offshellWmass_gen	 = ROOT.TH1F("hme_offshellWmass_gen","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_gen.GetXaxis().SetTitle("HME offshell W mass [GeV]");
offshellWmass_reco	 = ROOT.TH1F("hme_offshellWmass_reco","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_reco.GetXaxis().SetTitle("HME offshell W mass [GeV]");

# adding more branches to TTree
def initbr():
    hme_h2mass_gen[0] = -1
    hme_h2mass_reco[0] = -1
    hme_mean_reco[0] = -1
    hme_stddev_reco[0] = -1
    hme_entry_peak_reco[0] = -1
    hme_h2mass_weight2_reco[0] = -1
    hme_mean_weight2_reco[0] = -1
    hme_stddev_weight2_reco[0] = -1
    hme_entry_peak_weight2_reco[0] = -1
    hme_entries_reco[0] = -1

useGenParticles  = False

stop_watch = ROOT.TStopwatch()
stop_watch2 = ROOT.TStopwatch()

stop_watch2.Start()
for nEv in range(nStart, nEnd):
  if (doTest and nEv%10 == 0 ):
      print("nEv ",nEv)
  elif (nEv%10000 == 0):
      print("nEv ",nEv)
  if (doTest and nEv-nStart>=100):
      break

  if nEv == nStart:
      print("First event to run: ",nEv, " last event ", nEnd," HME iterations ",args.iterations)
  TCha.GetEntry(nEv)
  initbr()


  leadinglepPt = {"Muon": 26, "Electron": 30}
  subleadinglepPt = {"Muon": 10, "Electron": 10}
  allmu_index = []; leadingmus_index = []
  allel_index = []; leadingels_index = []
  leadingleppt = 26.0; subleadingleppt = 10.0;
  leadingmu_index = -1; subleadingmu_index = -1; 
  leadingel_index = -1; subleadingel_index = -1;
  leptonpair_type = -1
  lep1_pt = 0.0; lep1_eta = 0.0; lep1_phi = 0.0;
  lep2_pt = 0.0; lep2_eta = 0.0; lep2_phi = 0.0;
  for i, mupt in enumerate(TCha.muon_pt):
      if mupt > 10.0 and abs(TCha.muon_eta[i])<2.4 and abs(TCha.muon_d0[i]) < 0.1 and abs(TCha.muon_dz[i]) < 0.1 and abs(TCha.muon_sip3D[i]) < 4 and TCha.muon_miniIso[i] < 0.2:#add id cut
          #print "Muon index ",i, " mupt ",mupt," TCha.muon_pt ",TCha.muon_pt[i]
          allmu_index.append(i)
          if mupt > 26.0: leadingmus_index.append(i)
  for i, elpt in enumerate(TCha.electron_pt):
      if elpt > 10.0 and abs(TCha.electron_scEta[i])<2.5 and abs(TCha.electron_d0[i]) < 0.1 and abs(TCha.electron_dz[i]) < 0.1 and abs(TCha.electron_sip3D[i]) < 4 and TCha.electron_miniIso[i] < 0.2:#add id cut
          #print "El index ",i," elpt ",elpt," TCha.muon_pt ",TCha.electron_pt[i]
          allel_index.append(i)
          if elpt > 30.0: leadingels_index.append(i) 

  #print "allmu index ",allmu_index, " allel_index ",allel_index
  if len(leadingmus_index) >= 1 and len(allmu_index) >= 2:
      leptonpair_type = 1
      for i in leadingmus_index:
          if TCha.muon_pt[i] > leadingleppt:
              leadingleppt = TCha.muon_pt[i]
              leadingmu_index = i
      for j in allmu_index:
          if j != leadingmu_index and TCha.muon_pt[j] > subleadingleppt:
              subleadingleppt = TCha.muon_pt[j]
              subleadingmu_index = j
      lep1_pt = TCha.muon_pt[leadingmu_index]; lep1_eta = TCha.muon_eta[leadingmu_index]; lep1_phi = TCha.muon_phi[leadingmu_index];
      lep2_pt = TCha.muon_pt[subleadingmu_index]; lep2_eta = TCha.muon_eta[subleadingmu_index]; lep2_phi = TCha.muon_phi[subleadingmu_index];
      lep1index = leadingmu_index; lep2index = subleadingmu_index
  elif len(leadingmus_index) >= 1 and len(allel_index) >= 1:
      leptonpair_type = 2
      for i in leadingmus_index:
          if TCha.muon_pt[i] > leadingleppt:
              leadingleppt = TCha.muon_pt[i]
              leadingmu_index = i
      for j in allel_index:
          if TCha.electron_pt[j] > subleadingleppt:
              subleadingleppt = TCha.electron_pt[j]
              subleadingel_index = j
      lep1_pt = TCha.muon_pt[leadingmu_index]; lep1_eta = TCha.muon_eta[leadingmu_index]; lep1_phi = TCha.muon_phi[leadingmu_index];
      lep2_pt = TCha.electron_pt[subleadingel_index]; lep2_eta = TCha.electron_eta[subleadingel_index]; lep2_phi = TCha.electron_phi[subleadingel_index];
      lep1index = leadingmu_index; lep2index = subleadingel_index
  elif len(leadingels_index) >= 1 and len(allmu_index) >= 1:
      leptonpair_type = 3
      for i in leadingels_index:
          if TCha.electron_pt[i] > leadingleppt:
              leadingleppt = TCha.electron_pt[i]
              leadingel_index = i
      for j in allmu_index:
          if TCha.muon_pt[j] > subleadingleppt:
              subleadingleppt = TCha.muon_pt[j]
              subleadingmu_index = j
      lep1_pt = TCha.electron_pt[leadingel_index]; lep1_eta = TCha.electron_eta[leadingel_index]; lep1_phi = TCha.electron_phi[leadingel_index];
      lep2_pt = TCha.muon_pt[subleadingmu_index]; lep2_eta = TCha.muon_eta[subleadingmu_index]; lep2_phi = TCha.muon_phi[subleadingmu_index];
      lep1index = leadingel_index; lep2index = subleadingmu_index

  elif len(leadingels_index) >= 1 and len(allel_index) >= 2:
      leptonpair_type = 4
      for i in leadingels_index:
          if TCha.electron_pt[i] > leadingleppt:
              leadingleppt = TCha.electron_pt[i]
              leadingel_index = i
      for j in allel_index:
          if j != leadingel_index and TCha.electron_pt[j] > subleadingleppt:
              subleadingleppt = TCha.electron_pt[j]
              subleadingel_index = j
      lep1_pt = TCha.electron_pt[leadingel_index]; lep1_eta = TCha.electron_eta[leadingel_index]; lep1_phi = TCha.electron_phi[leadingel_index];
      lep2_pt = TCha.electron_pt[subleadingel_index]; lep2_eta = TCha.electron_eta[subleadingel_index]; lep2_phi = TCha.electron_phi[subleadingel_index];
      lep1index = leadingel_index; lep2index = subleadingel_index
  else:
      print("no dilepton pair is found")
      continue
  
  print "leptonpair type ",leptonpair_type," lep1 index ",lep1index," pt ",lep1_pt," eta ",lep1_eta," lep2 index ",lep2index," pt ",lep2_pt," eta ",lep2_eta
  leptonpairtype[0] = leptonpair_type
  ll_dR = sqrt((lep1_eta-lep2_eta)*(lep1_eta-lep2_eta) + (lep1_phi-lep2_phi)*(lep1_phi-lep2_phi))
  lep1_p4 	  = ROOT.TLorentzVector(); lep1_p4.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, 0)
  lep2_p4 	  = ROOT.TLorentzVector(); lep2_p4.SetPtEtaPhiM(lep2_pt, lep2_eta, lep2_phi, 0) 
  ll_p4 = lep1_p4 + lep2_p4
  ll_M = ll_p4.M()
  cleancut        = (ll_M<(91-15) and ll_M > 12 and ll_dR < 1.6)
  if not cleancut:
      print("failed in dilepton clean cut, ll_M ",ll_M," ll_dR ",ll_dR)
      TCha2.Fill()
      continue

  fatjet_index = -1; fatjet_pt = 200.0
  def deltaPhi(phi1, phi2):
      dphi = phi1-phi2
      if dphi > pi:
          dphi = dphi - 2*pi
      elif dphi < pi* (-1.0):
          dphi = dphi + pi
      return dphi
  for i in range(len(TCha.ak8PuppiJet_pt)):
      if fatjet_index >= 0 and TCha.ak8PuppiJet_pt[i]  <= fatjet_pt:
          continue
      bblep1_dR = sqrt((TCha.ak8PuppiJet_phi[i]-lep1_phi)*(TCha.ak8PuppiJet_phi[i]-lep1_phi)+(TCha.ak8PuppiJet_eta[i]-lep1_eta)*(TCha.ak8PuppiJet_eta[i]-lep1_eta))
      bblep2_dR = sqrt((TCha.ak8PuppiJet_phi[i]-lep2_phi)*(TCha.ak8PuppiJet_phi[i]-lep2_phi)+(TCha.ak8PuppiJet_eta[i]-lep2_eta)*(TCha.ak8PuppiJet_eta[i]-lep2_eta))
      subjetscut = (TCha.ak8PuppiJet_sj1_pt[i] > 20 and abs(TCha.ak8PuppiJet_sj1_eta[i])<2.4 and TCha.ak8PuppiJet_sj2_pt[i] > 20 and abs(TCha.ak8PuppiJet_sj2_eta[i])<2.4 and (TCha.ak8PuppiJet_sj1_csv[i] > 0.4941 or TCha.ak8PuppiJet_sj2_csv[i] > 0.4941))
      #print "ak8jet pt ",TCha.ak8PuppiJet_pt[i]," mass ",TCha.ak8PuppiJet_mass[i]," bblep1_dR ",bblep1_dR," bblep2_dR ",bblep2_dR, " subjetcut ",subjetscut," subjet1 b-tagging ",TCha.ak8PuppiJet_sj1_csv[i]," subjet2 b-tagging ",TCha.ak8PuppiJet_sj2_csv[i]
      if TCha.ak8PuppiJet_pt[i] > 200.0 and abs(deltaPhi(TCha.ak8PuppiJet_phi[i], ll_p4.Phi())) > 2.0 and bblep1_dR > 0.8 and bblep2_dR > 0.8 and subjetscut:   
          fatjet_index = i; fatjet_pt = TCha.ak8PuppiJet_pt[i] 
          
  if fatjet_index < 0:
      print("fatjet is not found!!! ",list(TCha.ak8PuppiJet_pt)," subjet1 ",list(TCha.ak8PuppiJet_sj1_pt)," subjet2 ",list(TCha.ak8PuppiJet_sj2_pt)()
      TCha2.Fill()
      continue
  if TCha.event_met_pt <= 40:
      print("MET is ",TCha.event_met_pt," less than 40 ")
      continue

  ak8jetindex[0] = fatjet_index 
  print("nEv ",nEv," ak8jets mass ",TCha.ak8PuppiJet_mass[fatjet_index]," 125.0/mass ", 125.0/TCha.ak8PuppiJet_mass[fatjet_index], " metpt ",TCha.event_met_pt)
  jet1_pt = TCha.ak8PuppiJet_pt[fatjet_index]; jet1_eta = TCha.ak8PuppiJet_eta[fatjet_index]; jet1_phi = TCha.ak8PuppiJet_phi[fatjet_index]; jet1_mass = TCha.ak8PuppiJet_mass[fatjet_index]
  jet2_pt = 0.0; jet2_eta = 0.0; jet2_phi = 0.0; jet2_mass = 0.0;
  vetoEvent  = False
  for i in range(len(TCha.ak4PuppiJet_pt)):
      ak4jet_ak8jet_dR = sqrt((TCha.ak4PuppiJet_eta[i]-jet1_eta)*(TCha.ak4PuppiJet_eta[i]-jet1_eta)+(TCha.ak4PuppiJet_phi[i]-jet1_phi)*(TCha.ak4PuppiJet_phi[i]-jet1_phi))
      if TCha.ak4PuppiJet_csv[i] > 0.4941 and ak4jet_ak8jet_dR > 1.2:
          vetoEvent = True
  if vetoEvent:
      print("veto this event because extra ak4jet with medium btag ")
      TCha2.Fill()
      continue

  jet1_p4 	  = ROOT.TLorentzVector(); jet1_p4.SetPtEtaPhiM(jet1_pt, jet1_eta, jet1_phi, jet1_mass)
  jet2_p4 	  = ROOT.TLorentzVector(); jet2_p4.SetPtEtaPhiM(jet2_pt, jet2_eta, jet2_phi, jet2_mass)
  met_vec2    = ROOT.TVector2();       met_vec2.SetMagPhi(TCha.event_met_pt, TCha.event_met_phi)


  
  if (cleancut and doHME):

      #print "MT2_1 ", MT2_1," MT2_2 ",MT2_2, " MT2 final ", mt2[0]," mt2_ll ",mt2_ll[0]," mt2_jj ",mt2_bb[0]
      hme = HeavyMassEstimator()
      hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
      hme.setIterations(args.iterations)
      hme.setBjetCorrectionType(0)
      hme.setMETCorrectionType(6)
      hme.showKinematic()
      hme.setDebug(False)
      hme.runHME()
      #hme.hme_offshellWmass.SetName("hme_offshellWmass_TCha.d_genlTCha.e"%nEv)
      if hme.hme_h2Mass.GetEntries() <= 0:
          print("NO solution found!!!!! ")
      elif hme.hme_h2Mass.GetEntries() >0 and  hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin()) < 249.0 :
          print("Num solutions ",hme.hme_h2Mass.GetEntries()," BUT the maximum is ",hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin()))
      #hme.hme_h2Mass.Print("ALL")

      if hme.hme_h2Mass.GetEntries()>0 and hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())>=250.0 :
          hme.hme_h2Mass.SetName("hme_h2Mass_ev%d_recolevel"%nEv)
          hme.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_recolevel"%nEv)
          hme.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_recolevel"%nEv)
          h_h2mass_weight_reco_sum.Add(hme.hme_h2Mass)
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
          print("hme_h2mass_reco[0] ",hme_h2mass_reco[0]," hme_h2mass_weight2_reco[0] ",hme_h2mass_weight2_reco[0]," entries ", hme_entries_reco[0])
          #offshell Wmass
          h_offshellWmass_recoh2mass = hme.hme_h2MassAndoffshellWmass.ProjectionY("h_offshellWmass_recoh2mass",hme.hme_h2Mass.GetMaximumBin()-5, hme.hme_h2Mass.GetMaximumBin()+5)
          h_offshellWmass_recoh2mass_weight2 = hme.hme_h2MassAndoffshellWmass_weight2.ProjectionY("h_offshellWmass_recoh2mass_weight2",hme.hme_h2MassWeight2.GetMaximumBin()-5, hme.hme_h2MassWeight2.GetMaximumBin()+5)
          hme_mostprob_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetXaxis().GetBinCenter(h_offshellWmass_recoh2mass.GetMaximumBin())
          hme_mean_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetMean()
          #hme_stddev_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetStdDev(1)
          hme_mostprob_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetXaxis().GetBinCenter(h_offshellWmass_recoh2mass_weight2.GetMaximumBin())

        
  # WEIGHTS
  hmecputime[0] =  stop_watch.CpuTime()
  stop_watch.Stop()
  nEv = nEv +1
  TCha2.Fill()


TCha2.Write()
f.Close()
stop_watch2.Stop()
print("stop_watch2 Cputime ",stop_watch2.CpuTime()," realtime ",stop_watch2.RealTime())
