import ROOT
import numpy as np
from math import *
import HardcodeREFPDF  as REFPDF


class HeavyMassEstimator(object):
    """ Heavy  Mass Estimator 
    Attributes: 
    	two muons lorentz vectors
        two bjet lorentz vectors 
        missing ET
        root fils contains all PDFs
    hmetree should contain full information from  HME
    """
    iterations = 100
    onshellnuptpdf = ROOT.TH1F()
    onshellWmasspdf = ROOT.TH1F()
    offshellWmasspdf = ROOT.TH1F()
    recobjetrescalec1pdf = ROOT.TH1F()
    onshellnuptpdf_flag = False
    onshellWmasspdf_flag = False
    offshellWmasspdf_flag = False
    recobjetrescalec1pdf_flag = False
    eta_gen  = np.zeros(1, dtype=float);   phi_gen  = np.zeros(1, dtype=float)
    wmass_gen =  np.zeros(1, dtype=float); hmass_gen = np.zeros(1, dtype=float)
    metpx_corr = np.zeros(1, dtype=float);  metpy_corr = np.zeros(1, dtype=float);
    nsolutions = np.zeros(1, dtype=int)
    leptonNuPair = np.zeros(1, dtype=bool)
    weight = np.zeros(1, dtype=float)
    weight1 = np.zeros(1, dtype=float)
    weight2 = np.zeros(1, dtype=float)
    weight3 = np.zeros(1, dtype=float)
    weight4 = np.zeros(1, dtype=float)
    l_onshellW_eta = np.zeros(1, dtype=float); l_onshellW_phi = np.zeros(1, dtype=float); l_onshellW_pt = np.zeros(1, dtype=float); l_onshellW_energy = np.zeros(1, dtype=float)
    l_offshellW_eta = np.zeros(1, dtype=float); l_offshellW_phi = np.zeros(1, dtype=float); l_offshellW_pt = np.zeros(1, dtype=float); l_offshellW_energy = np.zeros(1, dtype=float)
    nu_onshellW_eta = np.zeros(1, dtype=float); nu_onshellW_phi = np.zeros(1, dtype=float); nu_onshellW_pt = np.zeros(1, dtype=float); nu_onshellW_energy = np.zeros(1, dtype=float)
    nu_offshellW_eta = np.zeros(1, dtype=float); nu_offshellW_phi = np.zeros(1, dtype=float); nu_offshellW_pt = np.zeros(1, dtype=float); nu_offshellW_energy = np.zeros(1, dtype=float)
    onshellW_eta = np.zeros(1, dtype=float); onshellW_phi = np.zeros(1, dtype=float); onshellW_pt = np.zeros(1, dtype=float); onshellW_energy = np.zeros(1, dtype=float); onshellW_mass = np.zeros(1, dtype=float)
    offshellW_eta = np.zeros(1, dtype=float); offshellW_phi = np.zeros(1, dtype=float); offshellW_pt = np.zeros(1, dtype=float); offshellW_energy = np.zeros(1, dtype=float); offshellW_mass = np.zeros(1, dtype=float)
    met_pt = np.zeros(1, dtype=float); met_phi = np.zeros(1, dtype=float); met_px = np.zeros(1, dtype=float); met_py = np.zeros(1, dtype=float) 
    b1jet_eta = np.zeros(1, dtype=float); b1jet_phi = np.zeros(1, dtype=float); b1jet_pt = np.zeros(1, dtype=float); b1jet_energy = np.zeros(1, dtype=float)
    b2jet_eta = np.zeros(1, dtype=float); b2jet_phi = np.zeros(1, dtype=float); b2jet_pt = np.zeros(1, dtype=float); b2jet_energy = np.zeros(1, dtype=float)
    b1rescalefactor = np.zeros(1, dtype=float); b2rescalefactor = np.zeros(1, dtype=float)
    htoBB_eta = np.zeros(1, dtype=float); htoBB_phi = np.zeros(1, dtype=float); htoBB_pt = np.zeros(1, dtype=float); htoBB_energy = np.zeros(1, dtype=float)
    htoWW_eta = np.zeros(1, dtype=float); htoWW_phi = np.zeros(1, dtype=float); htoWW_pt = np.zeros(1, dtype=float); htoWW_energy = np.zeros(1, dtype=float)
    htoBB_mass = np.zeros(1, dtype=float); htoWW_mass = np.zeros(1, dtype=float);
    h2tohh_eta = np.zeros(1, dtype=float); h2tohh_phi = np.zeros(1, dtype=float); h2tohh_pt = np.zeros(1, dtype=float); h2tohh_energy = np.zeros(1, dtype=float); h2tohh_mass = np.zeros(1, dtype=float); h2tohh_mass_total = np.zeros(1, dtype=float)

    """
    hmetree.Branch('nsolutions', nsolutions, 'nsolutions/I')
    hmetree.Branch('l_onshellW_eta', l_onshellW_eta, 'l_onshellW_eta/D')  
    hmetree.Branch('l_onshellW_phi', l_onshellW_phi, 'l_onshellW_phi/D')  
    hmetree.Branch('l_onshellW_pt', l_onshellW_pt, 'l_onshellW_pt/D')  
    hmetree.Branch('l_onshellW_energy', l_onshellW_energy, 'l_onshellW_energy/D')  
    hmetree.Branch('l_offshellW_eta', l_offshellW_eta, 'l_offshellW_eta/D')  
    hmetree.Branch('l_offshellW_phi', l_offshellW_phi, 'l_offshellW_phi/D')  
    hmetree.Branch('l_offshellW_pt', l_offshellW_pt, 'l_offshellW_pt/D')  
    hmetree.Branch('l_offshellW_energy', l_offshellW_energy, 'l_offshellW_energy/D')  
    hmetree.Branch('nu_onshellW_eta', nu_onshellW_eta, 'nu_onshellW_eta/D')  
    hmetree.Branch('nu_onshellW_phi', nu_onshellW_phi, 'nu_onshellW_phi/D')  
    hmetree.Branch('nu_onshellW_pt', nu_onshellW_pt, 'nu_onshellW_pt/D')  
    hmetree.Branch('nu_onshellW_energy', nu_onshellW_energy, 'nu_onshellW_energy/D')  
    hmetree.Branch('nu_offshellW_eta', nu_offshellW_eta, 'nu_offshellW_eta/D')  
    hmetree.Branch('nu_offshellW_phi', nu_offshellW_phi, 'nu_offshellW_phi/D')  
    hmetree.Branch('nu_offshellW_pt', nu_offshellW_pt, 'nu_offshellW_pt/D')  
    hmetree.Branch('nu_offshellW_energy', nu_offshellW_energy, 'nu_offshellW_energy/D')  
    hmetree.Branch('onshellW_eta', onshellW_eta, 'onshellW_eta/D')  
    hmetree.Branch('onshellW_phi', onshellW_phi, 'onshellW_phi/D')  
    hmetree.Branch('onshellW_pt', onshellW_pt, 'onshellW_pt/D')  
    #hmetree.Branch('onshellW_energy', onshellW_energy, 'onshellW_energy/D')  
    hmetree.Branch('onshellW_mass', onshellW_mass, 'onshellW_mass/D')  
    hmetree.Branch('offshellW_eta', offshellW_eta, 'offshellW_eta/D')  
    hmetree.Branch('offshellW_phi', offshellW_phi, 'offshellW_phi/D')  
    hmetree.Branch('offshellW_pt', offshellW_pt, 'offshellW_pt/D')  
    #hmetree.Branch('offshellW_energy', offshellW_energy, 'offshellW_energy/D')  
    hmetree.Branch('offshellW_mass', offshellW_mass, 'offshellW_mass/D')  
    hmetree.Branch('met_pt', met_pt, 'met_pt/D')
    hmetree.Branch('met_phi', met_phi, 'met_phi/D')
    hmetree.Branch('met_px', met_px, 'met_px/D')
    hmetree.Branch('met_py', met_py, 'met_py/D')
    hmetree.Branch('metpx_corr', metpx_corr, 'metpx_corr/D')
    hmetree.Branch('metpy_corr', metpy_corr, 'metpy_corr/D')
    #hmetree.Branch('b1jet_eta', b1jet_eta, 'b1jet_eta/D')
    #hmetree.Branch('b1jet_phi', b1jet_phi, 'b1jet_phi/D')
    #hmetree.Branch('b1jet_pt', b1jet_pt, 'b1jet_pt/D')
    #hmetree.Branch('b1jet_energy', b1jet_energy, 'b1jet_energy/D')
    #hmetree.Branch('b2jet_eta', b2jet_eta, 'b2jet_eta/D')
    #hmetree.Branch('b2jet_phi', b2jet_phi, 'b2jet_phi/D')
    #hmetree.Branch('b2jet_pt', b2jet_pt, 'b2jet_pt/D')
    #hmetree.Branch('b2jet_energy', b2jet_energy, 'b2jet_energy/D')
    #hmetree.Branch('htoBB_eta', htoBB_eta, 'htoBB_eta/D')
    #hmetree.Branch('htoBB_phi', htoBB_phi, 'htoBB_phi/D')
    #hmetree.Branch('htoBB_pt', htoBB_pt, 'htoBB_pt/D')
    #hmetree.Branch('htoBB_energy', htoBB_energy, 'htoBB_energy/D')
    hmetree.Branch('htoBB_mass', htoBB_mass, 'htoBB_mass/D')
    #hmetree.Branch('htoWW_eta', htoWW_eta, 'htoWW_eta/D')
    #hmetree.Branch('htoWW_phi', htoWW_phi, 'htoWW_phi/D')
    #hmetree.Branch('htoWW_pt', htoWW_pt, 'htoWW_pt/D')
    #hmetree.Branch('htoWW_energy', htoWW_energy, 'htoWW_energy/D')
    hmetree.Branch('htoWW_mass', htoWW_mass, 'htoWW_mass/D')
    hmetree.Branch('h2tohh_eta', h2tohh_eta, 'h2tohh_eta/D')
    hmetree.Branch('h2tohh_phi', h2tohh_phi, 'h2tohh_phi/D')
    hmetree.Branch('h2tohh_pt', h2tohh_pt, 'h2tohh_pt/D')
    hmetree.Branch('h2tohh_energy', h2tohh_energy, 'h2tohh_energy/D')
    hmetree.Branch('h2tohh_mass', h2tohh_mass, 'h2tohh_mass/D')
    hmetree.Branch('leptonNuPair', leptonNuPair, 'leptonNuPair/B')

    #hmetree.Branch('htoBB_phi', htoBB_phi, 'htoBB_phi/D')
    #hmetree.Branch('htoBB_pt', htoBB_pt, 'htoBB_pt/D')
    #hmetree.Branch('htoBB_energy', htoBB_energy, 'htoBB_energy/D')
    hmetree.Branch('htoBB_mass', htoBB_mass, 'htoBB_mass/D')
    #hmetree.Branch('htoWW_eta', htoWW_eta, 'htoWW_eta/D')
    #hmetree.Branch('htoWW_phi', htoWW_phi, 'htoWW_phi/D')
    #hmetree.Branch('htoWW_pt', htoWW_pt, 'htoWW_pt/D')
    #hmetree.Branch('htoWW_energy', htoWW_energy, 'htoWW_energy/D')
    hmetree.Branch('htoWW_mass', htoWW_mass, 'htoWW_mass/D')
    hmetree.Branch('h2tohh_eta', h2tohh_eta, 'h2tohh_eta/D')
    hmetree.Branch('h2tohh_phi', h2tohh_phi, 'h2tohh_phi/D')
    hmetree.Branch('h2tohh_pt', h2tohh_pt, 'h2tohh_pt/D')
    hmetree.Branch('h2tohh_energy', h2tohh_energy, 'h2tohh_energy/D')
    hmetree.Branch('h2tohh_mass', h2tohh_mass, 'h2tohh_mass/D')
    hmetree.Branch('leptonNuPair', leptonNuPair, 'leptonNuPair/B')
    """

    def __init__(self):
	
	#print "  create a HeavyMassEstimator object "
	try:
	    self.onshellWmasspdf = REFPDF.onshellWmasspdf
	    self.onshellWmasspdf_flag = True
	except NameError:
	    print "Failed to find onshellWmasspdf from hardcorded REFPDF"
	    self.onshellWmasspdf_flag = False

	try:
	    self.offshellWmasspdf = REFPDF.offshellWmasspdf
	    self.offshellWmasspdf_flag = True
	except NameError:
	    print "Failed to find offshellWmasspdf from hardcorded REFPDF"
	    self.offshellWmasspdf_flag = False

	try:
	    self.onshellnuptpdf = REFPDF.onshellnuptpdf
	    self.onshellnuptpdf_flag = True
	except NameError:
	    print "Failed to find onshellnuptpdf from hardcorded REFPDF"
	    self.onshellnuptpdf_flag = False

	try:
	    self.recobjetrescalec1pdf = REFPDF.recobjetrescalec1pdfPU40
	    self.recobjetrescalec1pdf_flag = True
	except NameError:
	    print "Failed to find recobjetrescalec1pdf from hardcorded REFPDF"
	    self.recobjetrescalec1pdf_flag = False


        #PUSample: 25.2, PU0: 14.8
        self.met_sigma = 25.2
	self.met_covxx = 1.0; self.met_covyy = 1.0;  self.met_covxy = 0.0; ## by default, it is 
	self.met_covcorrection = False;


        minMass = 200.0;  maxMass = 4000.0; nbinsMass = 3800		
	self.debug = False
        self.hme_h2Mass = ROOT.TH1F("hme_h2Mass","h2 mass from HME", nbinsMass, minMass, maxMass)
        self.hme_h2Mass_divSolutions = ROOT.TH1F("hme_h2Mass_divSolutions", "h2 mass from HME / nSolutions", nbinsMass, minMass, maxMass)
        self.hme_h2Mass_correctmunupair = ROOT.TH1F("hme_h2Mass_correctmunupair","h2 mass from HME", nbinsMass, minMass, maxMass)
        self.hme_h2Mass_incorrectmunupair = ROOT.TH1F("hme_h2Mass_incorrectmunupair","h2 mass from HME", nbinsMass, minMass, maxMass)
	self.hme_h2MassWeight1 = ROOT.TH1F("hme_h2MassWeight1","h2 mass from HME", nbinsMass, minMass, maxMass)
	self.hme_h2MassWeight2 = ROOT.TH1F("hme_h2MassWeight2","h2 mass from HME", nbinsMass, minMass, maxMass)
	self.hme_h2MassWeight3 = ROOT.TH1F("hme_h2MassWeight3","h2 mass from HME", nbinsMass, minMass, maxMass)
	#self.hme_h2MassWeight4 = ROOT.TH1F("hme_h2MassWeight4","h2 mass from HME",nbinsMass, minMass, maxMass)
	self.hme_offshellWmass = ROOT.TH1F("hme_offshellWmass","offshell W mass from HME", 100, 0.0, 100.0)
	### offshell Wmass(y-axis) Vs HME 
	self.hme_h2MassAndoffshellWmass = ROOT.TH2F("hme_h2MassAndoffshellWmass","h2 Mass and offshell W mass from HME", nbinsMass, minMass, maxMass, 100, 0.0, 100.0)
	self.hme_h2MassAndoffshellWmass_correctmunupair = ROOT.TH2F("hme_h2MassAndoffshellWmass_correctmunupair","h2 Mass and offshell W mass from HME", nbinsMass, minMass, maxMass, 100, 0.0, 100.0)
	self.hme_h2MassAndoffshellWmass_weight1 = ROOT.TH2F("hme_h2MassAndoffshellWmass_weight1","h2 Mass and offshell W mass from HME", nbinsMass, minMass, maxMass, 100, 0.0, 100.0)
	self.hme_h2MassAndoffshellWmass_weight2 = ROOT.TH2F("hme_h2MassAndoffshellWmass_weight2","h2 Mass and offshell W mass from HME", nbinsMass, minMass, maxMass, 100, 0.0, 100.0)
	
	self.hme_nsolutions = ROOT.TH1F("hme_nsolutions","num of solutions", 5, 0, 5.0)

	self.lepton1_p4  = ROOT.TLorentzVector()
	self.lepton2_p4  = ROOT.TLorentzVector()
	self.b1jet_p4  = ROOT.TLorentzVector()
	self.b2jet_p4  = ROOT.TLorentzVector()
	self.dijet_p4  = ROOT.TLorentzVector()
	self.met = ROOT.TVector2()
	self.wmasshist = ROOT.TH1F()
	self.onshellnupthist = ROOT.TH1F()
	self.lepton1_onshellW_p4 = ROOT.TLorentzVector()
	self.lepton1_offshellW_p4 = ROOT.TLorentzVector()
	self.nu_onshellW_p4 = ROOT.TLorentzVector()
	self.nu_offshellW_p4 = ROOT.TLorentzVector()
	self.onshellW_p4 = ROOT.TLorentzVector()
	self.onshellW_p4 = ROOT.TLorentzVector()
	self.htoWW_p4 =  ROOT.TLorentzVector()
	self.htoBB_p4 =  ROOT.TLorentzVector()
	self.h2tohh_p4 = ROOT.TLorentzVector()

        self.hmetree = ROOT.TTree("hmetree","HME Tree") 
        self.hmetree.Branch('nsolutions',            self.nsolutions, 'nsolutions/I')
        self.hmetree.Branch('l_onshellW_eta',        self.l_onshellW_eta, 'l_onshellW_eta/D')  
        self.hmetree.Branch('l_onshellW_phi',        self.l_onshellW_phi, 'l_onshellW_phi/D')  
        self.hmetree.Branch('l_onshellW_pt',         self.l_onshellW_pt, 'l_onshellW_pt/D')  
        self.hmetree.Branch('l_onshellW_energy',     self.l_onshellW_energy, 'l_onshellW_energy/D')  
        self.hmetree.Branch('l_offshellW_eta',       self.l_offshellW_eta, 'l_offshellW_eta/D')  
        self.hmetree.Branch('l_offshellW_phi',       self.l_offshellW_phi, 'l_offshellW_phi/D')  
        self.hmetree.Branch('l_offshellW_pt',        self.l_offshellW_pt, 'l_offshellW_pt/D')  
        self.hmetree.Branch('l_offshellW_energy',    self.l_offshellW_energy, 'l_offshellW_energy/D')  
        self.hmetree.Branch('nu_onshellW_eta',       self.nu_onshellW_eta, 'nu_onshellW_eta/D')  
        self.hmetree.Branch('nu_onshellW_phi',       self.nu_onshellW_phi, 'nu_onshellW_phi/D')  
        self.hmetree.Branch('nu_onshellW_pt',        self.nu_onshellW_pt, 'nu_onshellW_pt/D')  
        self.hmetree.Branch('nu_onshellW_energy',    self.nu_onshellW_energy, 'nu_onshellW_energy/D')  
        self.hmetree.Branch('nu_offshellW_eta',      self.nu_offshellW_eta, 'nu_offshellW_eta/D')  
        self.hmetree.Branch('nu_offshellW_phi',      self.nu_offshellW_phi, 'nu_offshellW_phi/D')  
        self.hmetree.Branch('nu_offshellW_pt',       self.nu_offshellW_pt, 'nu_offshellW_pt/D')  
        self.hmetree.Branch('nu_offshellW_energy',   self.nu_offshellW_energy, 'nu_offshellW_energy/D')  
        self.hmetree.Branch('onshellW_eta',          self.onshellW_eta, 'onshellW_eta/D')  
        self.hmetree.Branch('onshellW_phi',          self.onshellW_phi, 'onshellW_phi/D')  
        self.hmetree.Branch('onshellW_pt',           self.onshellW_pt, 'onshellW_pt/D')  
        #self.hmetree.Branch('onshellW_energ         self.y', onshellW_energy, 'onshellW_energy/D')  
        self.hmetree.Branch('onshellW_mass',         self.onshellW_mass, 'onshellW_mass/D')  
        self.hmetree.Branch('offshellW_eta',         self.offshellW_eta, 'offshellW_eta/D')  
        self.hmetree.Branch('offshellW_phi',         self.offshellW_phi, 'offshellW_phi/D')  
        self.hmetree.Branch('offshellW_pt',          self.offshellW_pt, 'offshellW_pt/D')  
        #self.hmetree.Branch('offshellW_energ        self.y', offshellW_energy, 'offshellW_energy/D')  
        self.hmetree.Branch('offshellW_mass',        self.offshellW_mass, 'offshellW_mass/D')  
        self.hmetree.Branch('met_pt',                self.met_pt, 'met_pt/D')
        self.hmetree.Branch('met_phi',               self.met_phi, 'met_phi/D')
        self.hmetree.Branch('met_px',                self.met_px, 'met_px/D')
        self.hmetree.Branch('met_py',                self.met_py, 'met_py/D')
        self.hmetree.Branch('metpx_corr',            self.metpx_corr, 'metpx_corr/D')
        self.hmetree.Branch('metpy_corr',            self.metpy_corr, 'metpy_corr/D')
        #self.hmetree.Branch('b1jet_eta', b1jet_eta, 'b1jet_eta/D')
        #self.hmetree.Branch('b1jet_phi', b1jet_phi, 'b1jet_phi/D')
        #self.hmetree.Branch('b1jet_pt', b1jet_pt, 'b1jet_pt/D')
        #self.hmetree.Branch('b1jet_energy', b1jet_energy, 'b1jet_energy/D')
        #self.hmetree.Branch('b2jet_eta', b2jet_eta, 'b2jet_eta/D')
        #self.hmetree.Branch('b2jet_phi', b2jet_phi, 'b2jet_phi/D')
        #self.hmetree.Branch('b2jet_pt', b2jet_pt, 'b2jet_pt/D')
        #self.hmetree.Branch('b2jet_energy', b2jet_energy, 'b2jet_energy/D')
        #self.hmetree.Branch('htoBB_eta', htoBB_eta, 'htoBB_eta/D')
        #self.hmetree.Branch('htoBB_phi', htoBB_phi, 'htoBB_phi/D')
        #self.hmetree.Branch('htoBB_pt', htoBB_pt, 'htoBB_pt/D')
        #self.hmetree.Branch('htoBB_energy', htoBB_energy, 'htoBB_energy/D')
        #self.hmetree.Branch('htoWW_eta', htoWW_eta, 'htoWW_eta/D')
        #self.hmetree.Branch('htoWW_phi', htoWW_phi, 'htoWW_phi/D')
        #self.hmetree.Branch('htoWW_pt', htoWW_pt, 'htoWW_pt/D')
        #self.hmetree.Branch('htoWW_energy', htoWW_energy, 'htoWW_energy/D')
        self.hmetree.Branch('htoBB_mass',            self.htoBB_mass, 'htoBB_mass/D')
        self.hmetree.Branch('htoWW_mass',            self.htoWW_mass, 'htoWW_mass/D')
        self.hmetree.Branch('h2tohh_eta',            self.h2tohh_eta, 'h2tohh_eta/D')
        self.hmetree.Branch('h2tohh_phi',            self.h2tohh_phi, 'h2tohh_phi/D')
        self.hmetree.Branch('h2tohh_pt',             self.h2tohh_pt, 'h2tohh_pt/D')
        self.hmetree.Branch('h2tohh_energy',         self.h2tohh_energy, 'h2tohh_energy/D')
        self.hmetree.Branch('h2tohh_mass',           self.h2tohh_mass, 'h2tohh_mass/D')
        self.hmetree.Branch('h2tohh_mass_total',     self.h2tohh_mass_total, 'h2tohh_mass_total/D')
        self.hmetree.Branch('leptonNuPair',          self.leptonNuPair, 'leptonNuPair/B')

	#print "self.hme_h2Mass entries ",self.hme_h2Mass.GetEntries()," hmetree entries ",self.hmetree.GetEntries()
	#print "hme_h2MassWeight1 entries ",self.hme_h2MassWeight1.GetEntries()
	

    def setDebug(self, x):
	self.debug = x

    def setIterations(self, n):
	self.iterations = int(n)
    def setMETResolution(self, x):
	##if not set, then use the default value!
        #PUSample: 25.2, PU0: 14.8
	print("Setting the MET resolution to be: ", x, " the default is ",self.met_sigma)
        self.met_sigma = x 

    def setMETCovMatrix(self, covxx, covyy, covxy, covcorrection):
	self.met_covxx = covxx
	self.met_covyy = covyy
	self.met_covxy = covxy
	self.met_covcorrection = covcorrection ## True or False
	print("Setting MET CovMatrix to correct MET: covxx ", covxx," covxy ", covxy, " covyy ", covyy)

    def setKinematic(self, lepton1_p4, lepton2_p4, jet1_p4, jet2_p4, met, onshellW_mu):
	self.lepton1_p4 = lepton1_p4
	self.lepton2_p4 = lepton2_p4
	self.b1jet_p4 = jet1_p4
	self.b2jet_p4 = jet2_p4
	self.dijet_p4 = self.b1jet_p4 + self.b2jet_p4
	self.met = met
	self.onshellW_mu = onshellW_mu
	#self.RefPDFFileName = RefPDFFileName
	#self.RefPDFFile = ROOT.TFile(RefPDFFileName,"READ")

    def setKinematic_boosted(self, lepton1_p4, lepton2_p4, dijet_p4, met, onshellW_mu):
	self.lepton1_p4 = lepton1_p4
	self.lepton2_p4 = lepton2_p4
	self.dijet_p4 = dijet_p4
	self.b1jet_p4 = None
	self.b2jet_p4 = None
	self.met = met
	self.onshellW_mu = onshellW_mu
	#self.RefPDFFileName = RefPDFFileName
	#self.RefPDFFile = ROOT.TFile(RefPDFFileName,"READ")


    def setonshellWmasspdf(self, hist):
	self.onshellWmasspdf = hist
	self.onshellWmasspdf_flag = True

    def setoffshellWmasspdf(self, hist):
	self.offshellWmasspdf = hist
	self.offshellWmasspdf_flag = True
    
    def setonshellnuptpdf(self, hist):
	self.onshellnuptpdf = hist
	self.onshellnuptpdf_flag = True

    def setrecobjetrescalec1pdf(self, hist):
	self.recobjetrescalec1pdf = hist
	self.recobjetrescalec1pdf_flag = True

    def showKinematic(self):
	print "lepton1 ",self.lepton1_p4.Print()
	print "lepton2 ",self.lepton2_p4.Print()
	if (self.b1jet_p4 and self.b2jet_p4):
	    print "b1jet ",self.b1jet_p4.Print()
	    print "b2jet ",self.b2jet_p4.Print()
	print "dijet ",self.dijet_p4.Print()
	print "Met ",self.met.Print()

    def initHMETree(self):
	""" intialize the HME tree """
	self.l_onshellW_eta[0] = -9.0; self.l_onshellW_phi[0] = -9.0; self.l_onshellW_pt[0] = -1.0; self.l_onshellW_energy[0] = -1.0
	self.l_offshellW_eta[0] = -9.0; self.l_offshellW_phi[0] = -9.0; self.l_offshellW_pt[0] = -1.0; self.l_offshellW_energy[0] = -1.0
	self.nu_onshellW_eta[0] = -9.0; self.nu_onshellW_phi[0] = -9.0; self.nu_onshellW_pt[0] = -1.0; self.nu_onshellW_energy[0] = -1.0
	self.nu_offshellW_eta[0] = -9.0; self.nu_offshellW_phi[0] = -9.0; self.nu_offshellW_pt[0] = -1.0; self.nu_offshellW_energy[0] = -1.0
	self.onshellW_eta[0] = -9.0; self.onshellW_phi[0] = -9.0; self.onshellW_pt[0] = -1.0; self.onshellW_energy[0] = -1.0; self.onshellW_mass[0] = -1.0
	self.offshellW_eta[0] = -9.0; self.offshellW_phi[0] = -9.0; self.offshellW_pt[0] = -1.0; self.offshellW_energy[0] = -1.0;self.offshellW_mass[0] = -1.0
	self.htoWW_eta[0] = -9.0; self.htoWW_phi[0] = -9.0; self.htoWW_pt[0] = -1.0; self.htoWW_energy[0] = -1.0; self.htoWW_mass[0] = -1.0
	self.b1jet_eta[0] = -9.0; self.b1jet_phi[0] = -9.0; self.b1jet_pt[0] = -1.0; self.b1jet_energy[0] = -1.0
	self.b2jet_eta[0] = -9.0; self.b2jet_phi[0] = -9.0; self.b2jet_pt[0] = -1.0; self.b2jet_energy[0] = -1.0
	self.htoBB_eta[0] = -9.0; self.htoBB_phi[0] = -9.0; self.htoBB_pt[0] = -1.0; self.htoBB_energy[0] = -1.0; self.htoBB_mass[0] = -1.0
	self.h2tohh_eta[0] = -9.0; self.h2tohh_phi[0] = -9.0; self.h2tohh_pt[0] = -1.0; self.h2tohh_energy[0] = -1.0; self.h2tohh_mass[0] = -1.0; self.h2tohh_mass_total[0] = 0.0
	self.met_pt[0] = -1.0; self.met_px[0] = -99999.0; self.met_py[0] = -99999.0; self.met_phi[0] = -99999.0
	self.weight[0] = 1.0; self.weight1[0] = 1.0;  self.weight2[0] = 1.0; self.weight3[0] = 1.0; self.weight4[0] = 1.0
	self.b1rescalefactor[0] = 1.0; self.b2rescalefactor[0] = 1.0
	self.dijetrescalefactor = 1.0
	self.nsolutions[0] = 0

    def getWeightFromHist(self, hist, x):
        binx = hist.FindBin(x)
        if binx == 0 or binx == hist.GetNbinsX() + 1:
       	   return 0.0
        return hist.Interpolate(x)
   
    def getOnshellWMass(self, x0, step, random):
	### in higgs->WW, one W is close to onshell and with mass peaking at 80 GeV 
	### but the Wmass distribution has one long tail from 50GeV -80GeV
	### generate onshell W mass from profiled Wmass distribution for on shell W in Higgs->WW
	### profiled W mass distribution in self.onshellWmasspdf
	
	xmin = 50.0; xmax = 90.0
        while (x0 > xmax or x0 < xmin):
	    if x0 > xmax:
	    	x0 = x0 - xmax + xmin
	    if x0 < xmin:
	        x0 = xmax - (xmin - x0)
	x1 = x0 + step
	while (x1 > xmax or x1 < xmin):
	    if x1 > xmax:
	        x1 = x1 - xmax + xmin
	    if x1 < xmin:
	        x1 = xmax - (xmin - x1)
	w0  = self.onshellWmasspdf.Interpolate(x0)
	w1  = self.onshellWmasspdf.Interpolate(x1)
	#print "w0 ",w0," interpolate x0 ",self.onshellWmasspdf.Interpolate(x0)," w1 ",w1, " interpolate x1 ",self.onshellWmasspdf.Interpolate(x1)
	#w1/w0: transition probability 
	if (w1/w0 >= random):
	    return x1
	elif (w1/w0 < random):
	    return x0
        else:	
            print "error in getOnshellWMass "
            return 80.3
    
   
    def bjetsCorrection_dijet(self):
	m = self.dijet_p4.M()
	ratio = 125.0/m
	self.dijetrescalefactor  = ratio

    def metCorrection_dijet(self):
        if self.dijetrescalefactor < 0.0:
	    print "dibjet correction is working properly, self.dijetrescalefactor ",self.dijetrescalefactor
	    return ROOT.TVector2(0.0, 0.0)	
	metpx_tmp = - (self.dijetrescalefactor - 1.0)*self.dijet_p4.Px()
	metpy_tmp = - (self.dijetrescalefactor - 1.0)*self.dijet_p4.Py()
	return ROOT.TVector2(metpx_tmp, metpy_tmp)
    

    def bjetsCorrection(self):
	if not (self.recobjetrescalec1pdf_flag):
	    #print "failed to have recobjetrescalec1pdf_flag , no Correction"
	    return True
	rescalec1 = self.recobjetrescalec1pdf.GetRandom()
        leadingbjet_p4 = self.b1jet_p4
        trailingbjet_p4 = self.b2jet_p4
        b1jetleadingjet = True
        if self.b1jet_p4.Pt() < self.b2jet_p4.Pt():
	    b1jetleadingjet = False
	    leadingbjet_p4 = self.b2jet_p4
	    trailingbjet_p4 = self.b1jet_p4
	x1 = trailingbjet_p4.M2()
	x2 = 2*rescalec1*(leadingbjet_p4*trailingbjet_p4)
	x3 = rescalec1*rescalec1*leadingbjet_p4.M2() - 125.0*125.0
	if x2<0:
	    print "error bjets lorentzvector dot product less than 0"
	    return False
        if ((x2*x2 - 4*x1*x3) < 0 or x1 == 0):
	    print "error ! there is no soluations for bjetsCorrection "
	    return False
	rescalec2 = (-x2 + sqrt(x2*x2 - 4*x1*x3))/(2*x1)
	if rescalec2 < .0:
	    #print "negative rescalec2: ",rescalec2
	    return False
	#print "rescalec1 ",rescalec1," rescalec2 ",rescalec2
	if b1jetleadingjet:
	    self.b1rescalefactor[0] = rescalec1
	    self.b2rescalefactor[0] = rescalec2
	else:
	    self.b1rescalefactor[0] = rescalec2
	    self.b2rescalefactor[0] = rescalec1
        return True	

    def metCorrection(self):
        if self.b1rescalefactor[0] < 0.0 or self.b2rescalefactor[0] < 0.0:
	    print "bjet elaborate correction is working properly, b1rescalefactor ",self.b1rescalefactor[0]," b2rescalefactor ",self.b2rescalefactor[0]
	    return ROOT.TVector2(0.0, 0.0)	
	metpx_tmp = - (self.b1rescalefactor[0] - 1.0)*self.b1jet_p4.Px() - (self.b2rescalefactor[0] - 1.0)*self.b2jet_p4.Px()
	metpy_tmp = - (self.b1rescalefactor[0] - 1.0)*self.b1jet_p4.Py() - (self.b2rescalefactor[0] - 1.0)*self.b2jet_p4.Py()
	return ROOT.TVector2(metpx_tmp, metpy_tmp)

    def metSmearing_Cov(self, metx, mety, covxx, covyy, covxy, Niter):
	##reference: https://juanitorduz.github.io/multivariate_normal/
	d = 2
	mean = np.array([metx, mety]).reshape(2,1)
	covmatrix = np.array([[covxx, covxy], [covxy, covyy]])

	np.linalg.eigvals(covmatrix)

	epsilon = 0.01


	K = covmatrix + epsilon*np.identity(d)
	L = np.linalg.cholesky(K)
	testK = np.dot(L, np.transpose(L))
	print("K ", K, " testK ", testK)

	#Niter = 100
	u = np.random.normal(loc=0, scale=1, size=d*Niter).reshape(d, Niter)
	dx = np.dot(L, u)## met after smearing with covariant method


	for i in range(Niter):
		print("x : ", dx[0][i]," y: ", dx[1][i])
    
	return dx

    def assignMuP4(self, case):
	"""lepton+nu permutation
	    in simulation:
	    case 0: assgin lepton from onshell W to l_onshellW_p4
	    case 1: assgin lepton from offshell W to l_onshellW_p4
	   not in simuation: what ever comes first is assign to l_onshellW_p4
	"""
	if case == 0:
	    self.lepton_onshellW_p4 = self.lepton1_p4
	    self.lepton_offshellW_p4 = self.lepton2_p4
	    if self.onshellW_mu == 1:
	    	self.correctmunupair = True
	    else:
	    	self.correctmunupair = False

	elif case == 1:
	    self.lepton_offshellW_p4 = self.lepton1_p4
	    self.lepton_onshellW_p4 = self.lepton2_p4
	    if self.onshellW_mu == 2:
	    	self.correctmunupair = True
	    else:
	    	self.correctmunupair = False
	
    def nuPtFromOnshellW(self, nu_eta, nu_phi, lepton_p4, wMass):
	deta = nu_eta-lepton_p4.Eta()
	dphi = nu_phi-lepton_p4.Phi()
	nuPt = wMass*wMass/(2*lepton_p4.Pt()*(cosh(deta)-cos(dphi)))
	#nuPt = wMass*wMass/(2*lepton_p4.Pt()*(cosh(deta)-cosh(dphi)))

        return nuPt

    def nuP4FromOffshellW(self, met, lepton1_p4, lepton2_p4, nu1_p4, nu2_p4, case, hMass):
	tmp_p4 = lepton1_p4 + lepton2_p4 + nu1_p4
        tmp_nu_px = met.Px() - nu1_p4.Px()
        tmp_nu_py = met.Py() - nu1_p4.Py()
        nu_pxpy =  ROOT.TVector2(tmp_nu_px, tmp_nu_py)
        tmp_nu_pt = nu_pxpy.Mod()
        tmp_p4_v2 = ROOT.TLorentzVector(sqrt(pow(tmp_p4.Pt(), 2) + pow(tmp_p4.M(), 2)), 0, tmp_p4.Pz(), tmp_p4.Energy())
        chdeta = (pow(hMass, 2) + 2*(nu_pxpy.Px()*tmp_p4.Px() + nu_pxpy.Py()*tmp_p4.Py()) - pow(tmp_p4.M(), 2))/(2.0*tmp_p4_v2.Pt()*tmp_nu_pt)
        if chdeta < 1.0:
        #no solution if chdeta<1.0
    	    #print "no solution since chdeta<1.0, chdeta ",chdeta
	    nu2_p4.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0)
            return False 
	tmp_nu_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi())
        deta = acosh( chdeta )
        tmp_nu_eta = 0.0
        if case == 1:
	    tmp_nu_eta = tmp_p4_v2.Eta() - deta
	else :
	    tmp_nu_eta = tmp_p4_v2.Eta() + deta
        if (abs(tmp_nu_eta) > 6.0):
        #very unlikely solution
	    print "tmp_nu_eta ",tmp_nu_eta, " very unlikely solution, pass"
	    nu2_p4.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0)
            return False 
	nu2_p4.SetPtEtaPhiM(tmp_nu_pt, tmp_nu_eta, tmp_nu_phi, 0.0)	
        htoWW_tmp = tmp_p4 + nu2_p4
	if abs(htoWW_tmp.M() - hMass) > 1.0:
	    print "Warning!!! gen hmass ", hMass, " HME htoWW mass ", htoWW_tmp.M()
	return True

    def runHME(self):
	    
	if self.debug:
	    self.showKinematic()

	if not(self.onshellWmasspdf_flag):
	    print "no onshellWmasspdf, error!!! "
	    return  False
	self.eta_gen[0] = 0.0 ; self.phi_gen[0] = 0.0
	#initial wmass_gen
	self.wmass_gen[0] = 70.0
	self.met_px[0] = self.met.Px(); self.met_py[0] = self.met.Py()
	it = 0
	genRandom = ROOT.TRandom3(0)
	boosted_flg = (self.b1jet_p4 == None and self.b2jet_p4 == None and self.dijet_p4 != None)
        if not self.recobjetrescalec1pdf_flag and  not boosted_flg:
	    if self.debug:	print "self.recobjetrescalec1pdf_flag is False! and not boosted topology"
	    self.met_sigma = 0.0	
	
	hme_binwidth = self.hme_h2Mass.GetBinWidth(1)
        dxy_covcorrection = None
        if self.met_covcorrection:
	     dxy_covcorrection = self.metSmearing_Cov(self.met_px[0], self.met_py[0], self.met_covxx, self.met_covyy, self.met_covxy, self.iterations)
        #genRandom.SetSeed()
	while (it < self.iterations ):
	    it += 1
	    self.initHMETree()
	    self.eta_gen[0] = genRandom.Uniform(-6, 6)
	    self.phi_gen[0] = genRandom.Uniform(-3.1415926, 3.1415926)
	    self.hmass_gen[0] = genRandom.Gaus(125.03, 0.004)
	    rand01 = genRandom.Uniform(0., 1.0)
	    step = genRandom.Uniform(-10.0, 10.0)
	    self.wmass_gen[0] = self.getOnshellWMass(self.wmass_gen[0], step, rand01)
	    if self.debug:
		print "it ",it," self.eta_gen[0] ",self.eta_gen[0]," wmass_gen ",self.wmass_gen[0]
	    #smearing met 
            met_dpx = None; met_dpy = None
	    if self.met_covcorrection:
	        met_dpx = dxy_covcorrection[0][it]
	        met_dpy = dxy_covcorrection[1][it]
	    else:
		met_dpx = genRandom.Gaus(0.0, self.met_sigma)
		met_dpy = genRandom.Gaus(0.0, self.met_sigma)

            if boosted_flg:
	        self.bjetsCorrection_dijet()
		met_corr = self.met + ROOT.TVector2(met_dpx, met_dpy)+ self.metCorrection_dijet()
		self.htoBB_p4 = self.dijet_p4 * self.dijetrescalefactor
	    else:
	        ### resolved case
		if self.recobjetrescalec1pdf_flag:
		    while not self.bjetsCorrection():
			#print "fail to get bjetcorrection, try to get next one "
			pass
		    met_corr = self.met + ROOT.TVector2(met_dpx, met_dpy)+ self.metCorrection()
		else:
		    met_corr = self.met
		    self.b1rescalefactor[0] = 1.0
		    self.b2rescalefactor[0] = 1.0
		self.htoBB_p4 = self.b1jet_p4 * self.b1rescalefactor[0] + self.b2jet_p4 * self.b2rescalefactor[0]



	    self.metpx_corr[0]= met_corr.Px()
	    self.metpy_corr[0] = met_corr.Py()
    	    if self.debug:
		print "met_px ",self.met.Px()," met_py ",self.met.Py()," after correction px ",met_corr.Px()," py ",met_corr.Py()
	    self.nsolutions[0] = 0
	    isolution = 0
	    solutions = [False, False, False, False]
	    hme_values = [-1.0, -1.0, -1.0, -1.0]
	    trueweights = [1.0, 1.0, 1.0, 1.0]
	    #1. permutation 
	    #2. check nu_onshell_W pt
	    #3. solve the kinematics
	    #4. mark solution is ture if it is solved
	    #5. dump information into tree
	    while isolution < len(solutions):
		case = isolution/2
		#print "case ",case," isolution ",isolution
		self.assignMuP4(case)
	        nu_onshellW_pt_tmp = self.nuPtFromOnshellW(self.eta_gen[0], self.phi_gen[0], self.lepton_onshellW_p4, self.wmass_gen[0])
	        nu_onshellW_p4_tmp = ROOT.TLorentzVector()
	        nu_onshellW_p4_tmp.SetPtEtaPhiM(nu_onshellW_pt_tmp, self.eta_gen[0], self.phi_gen[0], 0)
	        nu_offshellW_p4_tmp = ROOT.TLorentzVector()
	        solutions[isolution] = self.nuP4FromOffshellW(met_corr, self.lepton_onshellW_p4, self.lepton_offshellW_p4, nu_onshellW_p4_tmp, nu_offshellW_p4_tmp, isolution%2, self.hmass_gen[0])
	 	if solutions[isolution]:
		    self.nsolutions[0] += 1
	 	isolution += 1 
	    isolution = 0
	    if self.nsolutions[0] == 0:#no solution in this iterations
	    	if self.debug:	print " no solution in this iterations "
	 	continue
	    self.weight[0] = 1.0/float(self.nsolutions[0])
	    while isolution < len(solutions):
		if not(solutions[isolution]):	
		    isolution += 1
		    continue
		case = isolution/2
		self.assignMuP4(case)
		self.nu_onshellW_pt[0] = self.nuPtFromOnshellW(self.eta_gen[0], self.phi_gen[0], self.lepton_onshellW_p4, self.wmass_gen[0])
		self.nu_onshellW_p4.SetPtEtaPhiM(self.nu_onshellW_pt[0], self.eta_gen[0], self.phi_gen[0], 0)
		self.nuP4FromOffshellW(met_corr, self.lepton_onshellW_p4, self.lepton_offshellW_p4, self.nu_onshellW_p4, self.nu_offshellW_p4, isolution%2, self.hmass_gen[0])
	        self.onshellW_p4 = self.lepton_onshellW_p4+self.nu_onshellW_p4
	        self.offshellW_p4 = self.lepton_offshellW_p4+self.nu_offshellW_p4
		self.htoWW_p4 = self.onshellW_p4 + self.offshellW_p4
		self.h2tohh_p4 = self.htoWW_p4 + self.htoBB_p4

		if (fabs(self.htoWW_p4.M() - self.hmass_gen[0])>1.0):
		    print "Error!! hmass_gen ", self.hmass_gen[0], " higgs mass from HME htoWW_p4 ", self.htoWW_p4.M()
		if self.debug:
		    print "it ",it, " get this h2tohh_mass ",self.h2tohh_p4.M()," higgs from WW ", self.htoWW_p4.M()," higgs from bb ", self.htoBB_p4.M()
		self.l_onshellW_eta[0] = self.lepton_onshellW_p4.Eta()
		self.l_onshellW_phi[0] = self.lepton_onshellW_p4.Phi()
		self.l_onshellW_pt[0] = self.lepton_onshellW_p4.Pt()
		self.l_onshellW_energy[0] = self.lepton_onshellW_p4.Energy()
		self.l_offshellW_eta[0] = self.lepton_offshellW_p4.Eta()
		self.l_offshellW_phi[0] = self.lepton_offshellW_p4.Phi()
		self.l_offshellW_pt[0] = self.lepton_offshellW_p4.Pt()
		self.l_offshellW_energy[0] = self.lepton_offshellW_p4.Energy()

		self.nu_onshellW_eta[0] = self.nu_onshellW_p4.Eta()
		self.nu_onshellW_phi[0] = self.nu_onshellW_p4.Phi()
		self.nu_onshellW_pt[0] = self.nu_onshellW_p4.Pt()
		self.nu_onshellW_energy[0] = self.nu_onshellW_p4.Energy()
		self.nu_offshellW_eta[0] = self.nu_offshellW_p4.Eta()
		self.nu_offshellW_phi[0] = self.nu_offshellW_p4.Phi()
		self.nu_offshellW_pt[0] = self.nu_offshellW_p4.Pt()
		self.nu_offshellW_energy[0] = self.nu_offshellW_p4.Energy()

		self.onshellW_eta[0] = self.onshellW_p4.Eta()
		self.onshellW_phi[0] = self.onshellW_p4.Phi()
		self.onshellW_pt[0] = self.onshellW_p4.Pt()
		self.onshellW_energy[0] = self.onshellW_p4.Energy()
                self.onshellW_mass[0] = self.wmass_gen[0] 
		self.offshellW_eta[0] = self.offshellW_p4.Eta()
		self.offshellW_phi[0] = self.offshellW_p4.Phi()
		self.offshellW_pt[0] = self.offshellW_p4.Pt()
		self.offshellW_energy[0] = self.offshellW_p4.Energy()
    		self.offshellW_mass[0] = self.offshellW_p4.M()

    		self.htoWW_eta[0] = self.htoWW_p4.Eta()
    		self.htoWW_phi[0] = self.htoWW_p4.Phi()
    		self.htoWW_pt[0] = self.htoWW_p4.Pt()
    		self.htoWW_energy[0] = self.htoWW_p4.Energy()
    		self.htoWW_mass[0] = self.htoWW_p4.M()
    		self.htoBB_eta[0] = self.htoBB_p4.Eta()
    		self.htoBB_phi[0] = self.htoBB_p4.Phi()
    		self.htoBB_pt[0] = self.htoBB_p4.Pt()
    		self.htoBB_energy[0] = self.htoBB_p4.Energy()
    		self.htoBB_mass[0] = self.htoBB_p4.M()
    		self.h2tohh_pt[0] = self.h2tohh_p4.Pt()
    		self.h2tohh_energy[0] = self.h2tohh_p4.Energy()
    		self.h2tohh_mass[0] = self.h2tohh_p4.M()
                self.h2tohh_mass_total[0] += self.h2tohh_mass[0]

                ##check the offshell Wmass and h->WW mass.
    		if (self.offshellW_mass[0]>self.htoWW_mass[0]/2.0):
		    #print "self.htoWW_mass[0] ",self.htoWW_mass[0]," self.offshellW_mass[0] ",self.offshellW_mass[0]," weight from offshellW mass ",self.weight2[0]," onshell wmass ", self.onshellW_mass[0]
		    isolution += 1
		    continue

		if (self.h2tohh_p4.Pt()/self.h2tohh_p4.E() <.00000001):
		    print "Strange case: h2tohh pt ", self.h2tohh_p4.Pt(), " energy ",self.h2tohh_p4.E()
		    self.h2tohh_eta[0] = 1000000.0; self.h2tohh_phi[0] = 0.0
		else:
		    self.h2tohh_eta[0] = self.h2tohh_p4.Eta(); self.h2tohh_phi[0] = self.h2tohh_p4.Phi()


		self.hme_h2Mass.Fill(self.h2tohh_mass[0], self.weight[0])
		hme_values[isolution] = self.h2tohh_mass[0]
		N_solution_withinOnebin = 0
		ii = 0; jj=0
		while ii < isolution:
		    if fabs(self.h2tohh_mass[0] - hme_values[ii]) <= hme_binwidth: N_solution_withinOnebin += 1
		    ii += 1
		trueweight_this = 1.0/(N_solution_withinOnebin+1.0)
		while jj <= isolution:
		    if fabs(self.h2tohh_mass[0] - hme_values[jj]) <= hme_binwidth: trueweights[jj] = trueweight_this
		    jj += 1
		        
    		self.hme_offshellWmass.Fill(self.offshellW_mass[0], self.weight[0])

                ########################################
                #### the following is NOT used in finding most probable mass from HME
                ########################################
            
	        #if self.onshellnuptpdf_flag:
    		self.hme_h2MassAndoffshellWmass.Fill(self.h2tohh_mass[0], self.offshellW_mass[0], self.weight[0])

    		##weight1 results
		self.weight1[0] = self.weight[0] * self.getWeightFromHist(self.onshellnuptpdf, self.nu_onshellW_pt[0]) 
    		self.hme_h2MassWeight1.Fill(self.h2tohh_mass[0], self.weight1[0])
    		self.hme_h2MassAndoffshellWmass_weight1.Fill(self.h2tohh_mass[0], self.offshellW_mass[0], self.weight1[0])

    		##weight2 results
                self.weight2[0] = self.weight[0] * self.getWeightFromHist(self.offshellWmasspdf, self.offshellW_mass[0])
    		self.hme_h2MassWeight2.Fill(self.h2tohh_mass[0], self.weight2[0])
                self.hme_offshellWmass.Fill(self.offshellW_mass[0])
    		self.hme_h2MassAndoffshellWmass_weight2.Fill(self.h2tohh_mass[0], self.offshellW_mass[0], self.weight2[0])


		self.leptonNuPair[0]  = self.correctmunupair
    		if self.correctmunupair:
		    self.hme_h2Mass_correctmunupair.Fill(self.h2tohh_mass[0],  self.weight[0])
		    self.hme_h2MassAndoffshellWmass_correctmunupair.Fill(self.h2tohh_mass[0], self.offshellW_mass[0], self.weight[0])
    		else:
		    self.hme_h2Mass_incorrectmunupair.Fill(self.h2tohh_mass[0],  self.weight[0])



    		#self.hmetree.Fill()
	 	isolution += 1 
	    self.hme_h2Mass_divSolutions.Fill(self.h2tohh_mass_total[0]/self.nsolutions[0])
            self.hme_nsolutions.Fill(self.nsolutions[0]+0.5)
            #################################################################
            ## the following apply the reco mass dependent weight to each solution
            #################################################################
            #isolution = 0
	    #while isolution < len(solutions):
	    #    if not(solutions[isolution]):	
	    #        isolution += 1
	    #        continue
	    #    self.hme_h2Mass.Fill(hme_values[isolution], trueweights[isolution])
	    #    isolution += 1 
            #print("hme_values ", hme_values, " trueweights ", trueweights)
	    #print("self.h2tohh_mass_total ", self.h2tohh_mass_total[0], " nsolutions ", self.nsolutions[0])
	##### end of iteration






