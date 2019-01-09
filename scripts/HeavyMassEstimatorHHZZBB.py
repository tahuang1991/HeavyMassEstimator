import ROOT
import numpy as np
from math import *


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
    onshellZmasspdf = ROOT.TH1F()
    offshellZmasspdf = ROOT.TH1F()
    recobjetrescalec1pdf = ROOT.TH1F()
    onshellnuptpdf_flag = False
    onshellZmasspdf_flag = False
    recobjetrescalec1pdf_flag = False
    hmetree = ROOT.TTree("hmetree","HME Tree") 
    Zmass_gen =  np.zeros(1, dtype=float); hmass_gen = np.zeros(1, dtype=float)
    metpx_corr = np.zeros(1, dtype=float);  metpy_corr = np.zeros(1, dtype=float);
    nsolutions = np.zeros(1, dtype=int)
    weight = np.zeros(1, dtype=float)
    weight1 = np.zeros(1, dtype=float)
    weight2 = np.zeros(1, dtype=float)
    weight3 = np.zeros(1, dtype=float)
    weight4 = np.zeros(1, dtype=float)
    lepton1_eta = np.zeros(1, dtype=float); lepton1_phi = np.zeros(1, dtype=float); lepton1_pt = np.zeros(1, dtype=float); lepton1_energy = np.zeros(1, dtype=float)
    lepton2_eta = np.zeros(1, dtype=float); lepton2_phi = np.zeros(1, dtype=float); lepton2_pt = np.zeros(1, dtype=float); lepton2_energy = np.zeros(1, dtype=float)
    Zll_eta = np.zeros(1, dtype=float); Zll_phi = np.zeros(1, dtype=float); Zll_pt = np.zeros(1, dtype=float); Zll_energy = np.zeros(1, dtype=float); Zll_mass = np.zeros(1, dtype=float)
    Znunu_eta = np.zeros(1, dtype=float); Znunu_phi = np.zeros(1, dtype=float); Znunu_pt = np.zeros(1, dtype=float); Znunu_energy = np.zeros(1, dtype=float); Znunu_mass = np.zeros(1, dtype=float); Znunu_pz = np.zeros(1, dtype=float)
    met_pt = np.zeros(1, dtype=float); met_phi = np.zeros(1, dtype=float); met_px = np.zeros(1, dtype=float); met_py = np.zeros(1, dtype=float) 
    b1jet_eta = np.zeros(1, dtype=float); b1jet_phi = np.zeros(1, dtype=float); b1jet_pt = np.zeros(1, dtype=float); b1jet_energy = np.zeros(1, dtype=float)
    b2jet_eta = np.zeros(1, dtype=float); b2jet_phi = np.zeros(1, dtype=float); b2jet_pt = np.zeros(1, dtype=float); b2jet_energy = np.zeros(1, dtype=float)
    b1rescalefactor = np.zeros(1, dtype=float); b2rescalefactor = np.zeros(1, dtype=float)
    htoBB_eta = np.zeros(1, dtype=float); htoBB_phi = np.zeros(1, dtype=float); htoBB_pt = np.zeros(1, dtype=float); htoBB_energy = np.zeros(1, dtype=float)
    htoZZ_eta = np.zeros(1, dtype=float); htoZZ_phi = np.zeros(1, dtype=float); htoZZ_pt = np.zeros(1, dtype=float); htoZZ_energy = np.zeros(1, dtype=float)
    htoBB_mass = np.zeros(1, dtype=float); htoZZ_mass = np.zeros(1, dtype=float);
    h2tohh_eta = np.zeros(1, dtype=float); h2tohh_phi = np.zeros(1, dtype=float); h2tohh_pt = np.zeros(1, dtype=float); h2tohh_energy = np.zeros(1, dtype=float); h2tohh_mass = np.zeros(1, dtype=float)

    hmetree.Branch('nsolutions', nsolutions, 'nsolutions/I')
    hmetree.Branch('lepton1_eta', lepton1_eta, 'lepton1_eta/D')  
    hmetree.Branch('lepton1_phi', lepton1_phi, 'lepton1_phi/D')  
    hmetree.Branch('lepton1_pt', lepton1_pt, 'lepton1_pt/D')  
    hmetree.Branch('lepton1_energy', lepton1_energy, 'lepton1_energy/D')  
    hmetree.Branch('lepton2_eta', lepton2_eta, 'lepton2_eta/D')  
    hmetree.Branch('lepton2_phi', lepton2_phi, 'lepton2_phi/D')  
    hmetree.Branch('lepton2_pt', lepton2_pt, 'lepton2_pt/D')  
    hmetree.Branch('lepton2_energy', lepton2_energy, 'lepton2_energy/D')  
    hmetree.Branch('Zll_eta', Zll_eta, 'Zll_eta/D')  
    hmetree.Branch('Zll_phi', Zll_phi, 'Zll_phi/D')  
    hmetree.Branch('Zll_pt', Zll_pt, 'Zll_pt/D')  
    hmetree.Branch('Zll_energy', Zll_energy, 'Zll_energy/D')  
    hmetree.Branch('Zll_mass', Zll_mass, 'Zll_mass/D')  
    hmetree.Branch('Znunu_eta', Znunu_eta, 'Znunu_eta/D')  
    hmetree.Branch('Znunu_phi', Znunu_phi, 'Znunu_phi/D')  
    hmetree.Branch('Znunu_pt', Znunu_pt, 'Znunu_pt/D')  
    hmetree.Branch('Znunu_pz', Znunu_pz, 'Znunu_pz/D')  
    hmetree.Branch('Znunu_energy', Znunu_energy, 'Znunu_energy/D')  
    hmetree.Branch('Znunu_mass', Znunu_mass, 'Znunu_mass/D')  
    hmetree.Branch('met_pt', met_pt, 'met_pt/D')
    hmetree.Branch('met_phi', met_phi, 'met_phi/D')
    hmetree.Branch('met_px', met_px, 'met_px/D')
    hmetree.Branch('met_py', met_py, 'met_py/D')
    hmetree.Branch('metpx_corr', metpx_corr, 'metpx_corr/D')
    hmetree.Branch('metpy_corr', metpy_corr, 'metpy_corr/D')
    hmetree.Branch('b1jet_eta', b1jet_eta, 'b1jet_eta/D')
    hmetree.Branch('b1jet_phi', b1jet_phi, 'b1jet_phi/D')
    hmetree.Branch('b1jet_pt', b1jet_pt, 'b1jet_pt/D')
    hmetree.Branch('b1jet_energy', b1jet_energy, 'b1jet_energy/D')
    hmetree.Branch('b2jet_eta', b2jet_eta, 'b2jet_eta/D')
    hmetree.Branch('b2jet_phi', b2jet_phi, 'b2jet_phi/D')
    hmetree.Branch('b2jet_pt', b2jet_pt, 'b2jet_pt/D')
    hmetree.Branch('b2jet_energy', b2jet_energy, 'b2jet_energy/D')
    hmetree.Branch('htoBB_eta', htoBB_eta, 'htoBB_eta/D')
    hmetree.Branch('htoBB_phi', htoBB_phi, 'htoBB_phi/D')
    hmetree.Branch('htoBB_pt', htoBB_pt, 'htoBB_pt/D')
    hmetree.Branch('htoBB_energy', htoBB_energy, 'htoBB_energy/D')
    hmetree.Branch('htoBB_mass', htoBB_mass, 'htoBB_mass/D')
    hmetree.Branch('htoZZ_eta', htoZZ_eta, 'htoZZ_eta/D')
    hmetree.Branch('htoZZ_phi', htoZZ_phi, 'htoZZ_phi/D')
    hmetree.Branch('htoZZ_pt', htoZZ_pt, 'htoZZ_pt/D')
    hmetree.Branch('htoZZ_energy', htoZZ_energy, 'htoZZ_energy/D')
    hmetree.Branch('htoZZ_mass', htoZZ_mass, 'htoZZ_mass/D')
    hmetree.Branch('h2tohh_eta', h2tohh_eta, 'h2tohh_eta/D')
    hmetree.Branch('h2tohh_phi', h2tohh_phi, 'h2tohh_phi/D')
    hmetree.Branch('h2tohh_pt', h2tohh_pt, 'h2tohh_pt/D')
    hmetree.Branch('h2tohh_energy', h2tohh_energy, 'h2tohh_energy/D')
    hmetree.Branch('h2tohh_mass', h2tohh_mass, 'h2tohh_mass/D')





    def __init__(self):
	print "  create a HeavyMassEstimator object "
        self.hme_h2Mass = ROOT.TH1F("hme_h2Mass","h2 mass from HME",1000, 200.0,1200.0)
	self.hme_h2MassWeight1 = ROOT.TH1F("hme_h2MassWeight1","h2 mass from HME",1000, 200.0,1200.0)
	self.hme_h2MassWeight4 = ROOT.TH1F("hme_h2MassWeight4","h2 mass from HME",1000, 200.0,1200.0)

	self.dilepton_mass = np.zeros(1, dtype=float) 
	self.lepton1_p4  = ROOT.TLorentzVector()
	self.lepton2_p4  = ROOT.TLorentzVector()
        self.dilepton_p4 = ROOT.TLorentzVector()
	self.b1jet_p4  = ROOT.TLorentzVector()
	self.b2jet_p4  = ROOT.TLorentzVector()
	self.met = ROOT.TVector2()
	self.Zll_p4 = ROOT.TLorentzVector()
	self.Znunu_p4 = ROOT.TLorentzVector()
	self.htoZZ_p4 =  ROOT.TLorentzVector()
	self.htoBB_p4 =  ROOT.TLorentzVector()
	self.h2tohh_p4 = ROOT.TLorentzVector()
	#print "self.hme_h2Mass entries ",self.hme_h2Mass.GetEntries()," hmetree entries ",self.hmetree.GetEntries()
	#print "hme_h2MassWeight1 entries ",self.hme_h2MassWeight1.GetEntries()
	

    def setIterations(self, n):
	self.iterations = n

    def setKinematic(self, lepton1_p4, lepton2_p4, jet1_p4, jet2_p4, met):
	self.lepton1_p4 = lepton1_p4
	self.lepton2_p4 = lepton2_p4
	self.b1jet_p4 = jet1_p4
	self.b2jet_p4 = jet2_p4
	self.met = met
	self.dilepton_p4 = self.lepton1_p4 + self.lepton2_p4
	self.lepton1_eta[0] = self.lepton1_p4.Eta(); self.lepton1_phi[0] = self.lepton1_p4.Phi(); self.lepton1_pt[0] = self.lepton1_p4.Pt();self.lepton1_energy[0] = self.lepton1_p4.Energy()
	self.lepton2_eta[0] = self.lepton2_p4.Eta(); self.lepton2_phi[0] = self.lepton2_p4.Phi(); self.lepton2_pt[0] = self.lepton2_p4.Pt();self.lepton2_energy[0] = self.lepton2_p4.Energy()
	
	#self.RefPDFFileName = RefPDFFileName
	#self.RefPDFFile = ROOT.TFile(RefPDFFileName,"READ")

    def setonshellZmasspdf(self, hist):
	self.onshellZmasspdf = hist
	self.onshellZmasspdf_flag = True
    
    def setoffshellZmasspdf(self, hist):
	self.offshellZmasspdf = hist
	self.offshellZmasspdf_flag = True
    
    def setrecobjetrescalec1pdf(self, hist):
	self.recobjetrescalec1pdf = hist
	self.recobjetrescalec1pdf_flag = True

    def showKinematic(self):
	print "lepton1 ",self.lepton1_p4.Print()
	print "lepton2 ",self.lepton2_p4.Print()
	print "dileptons ",self.dilepton_p4.Print()
	print "b1jet ",self.b1jet_p4.Print()
	print "b2jet ",self.b2jet_p4.Print()
	print "Met ",self.met.Print()

    def initHMETree(self):
	""" intialize the HME tree """
	self.Zll_eta[0] = -9.0; self.Zll_phi[0] = -9.0; self.Zll_pt[0] = -1.0; self.Zll_energy[0] = -1.0; self.Zll_mass[0] = -1.0
	self.Znunu_eta[0] = -9.0; self.Znunu_phi[0] = -9.0; self.Znunu_pt[0] = -1.0; self.Znunu_energy[0] = -1.0
	self.htoZZ_eta[0] = -9.0; self.htoZZ_phi[0] = -9.0; self.htoZZ_pt[0] = -1.0; self.htoZZ_energy[0] = -1.0; self.htoZZ_mass[0] = -1.0
	self.b1jet_eta[0] = -9.0; self.b1jet_phi[0] = -9.0; self.b1jet_pt[0] = -1.0; self.b1jet_energy[0] = -1.0
	self.b2jet_eta[0] = -9.0; self.b2jet_phi[0] = -9.0; self.b2jet_pt[0] = -1.0; self.b2jet_energy[0] = -1.0
	self.htoBB_eta[0] = -9.0; self.htoBB_phi[0] = -9.0; self.htoBB_pt[0] = -1.0; self.htoBB_energy[0] = -1.0; self.htoBB_mass[0] = -1.0
	self.h2tohh_eta[0] = -9.0; self.h2tohh_phi[0] = -9.0; self.h2tohh_pt[0] = -1.0; self.h2tohh_energy[0] = -1.0; self.h2tohh_mass[0] = -1.0
	self.met_pt[0] = -1.0; self.met_px[0] = -99999.0; self.met_py[0] = -99999.0; self.met_phi[0] = -99999.0
	self.weight[0] = 1.0; self.weight1[0] = 1.0;  self.weight2[0] = 1.0; self.weight3[0] = 1.0; self.weight4[0] = 1.0
	self.b1rescalefactor[0] = 1.0; self.b2rescalefactor[0] = 1.0
	self.nsolutions[0] = 0

    def getWeightFromHist(self, hist, x):
        binx = hist.FindBin(x)
        if binx == 0 or binx == hist.GetNbinsX() + 1:
       	   return 0.0
        return hist.Interpolate(x)
   
    def getOnshellZMass(self, x0, step, random):
	
	xmin = 50.0; xmax = 110.0
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
	w0  = self.onshellZmasspdf.Interpolate(x0)
	w1  = self.onshellZmasspdf.Interpolate(x1)
	
	#print "OnshellZmass x0 ",x0," step ", step, " random ", random," w0 ",w0," w1 ",w1," w1/w0 ",w1/w0
	#w1/w0: transition probability 
	if (w1/w0 >= random):
	    return x1
	elif (w1/w0 < random):
	    return x0
        else:	
            print "error in getOnshellZMass "
            return 91.18

    def getOffshellZMass(self, x0, step, random):
	
	xmin = 10.0; xmax = 60.0
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
	w0  = self.offshellZmasspdf.Interpolate(x0)
	w1  = self.offshellZmasspdf.Interpolate(x1)
	#print "OffshellZmass x0 ",x0," step ", step, " random ", random," w0 ",w0," w1 ",w1," w1/w0 ",w1/w0
	#w1/w0: transition probability 
	if (w1/w0 >= random):
	    return x1
	elif (w1/w0 < random):
	    return x0
        else:	
            print "error in getOffshellWMass "
            return 30.18
    
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
    
	
    def ZtonunuPzAndE(self, dilepton_p4, met, zMass, hMass, case):
	mll = dilepton_p4.M()
	Ell = dilepton_p4.Energy()
	Pzll = dilepton_p4.Pz()
        a1 = met.Mod2()+ zMass*zMass
	a2 = (hMass*hMass - zMass*zMass - mll*mll)/2.0 + met.Px()*dilepton_p4.Px() + met.Py()*dilepton_p4.Py()
	a3 =  4*a2*a2*Pzll*Pzll - 4*(Ell*Ell*a1 - a2*a2)*(Ell*Ell - Pzll*Pzll)
	#print "mll ",mll," Ell ",Ell," Pzll ", Pzll," metx ",met.Px()," mety ",met.Py(), " zmass ",zMass, " hmass ",hMass," a1 ",a1, " a2 ",a2, " a3 ",a3
	if a3>0 and case==0:
	   Pznunu = (2*a2*Pzll + sqrt(a3))/(2*(Ell*Ell-Pzll*Pzll))
	   Enunu = (a2 + Pznunu*dilepton_p4.Pz())/dilepton_p4.Energy()
	   if Enunu<0.0:
		print "warning Pz(nunu) ",Pznunu," E(nunu) ",Enunu," E is not positive!!!!"
	   return Pznunu, Enunu
	elif a3>0 and case==1:
	   Pznunu = (2*a2*Pzll - sqrt(a3))/(2*(Ell*Ell-Pzll*Pzll))
	   Enunu = (a2 + Pznunu*dilepton_p4.Pz())/dilepton_p4.Energy()
	   if Enunu<0.0:
		print "warning Pz(nunu) ",Pznunu," E(nunu) ",Enunu," E is not positive!!!!"
	   return Pznunu, Enunu
	else:
	   return 0.0, -1.0


    def runHME(self):
	if not(self.onshellZmasspdf_flag and self.offshellZmasspdf_flag):
	    print "no onshellZmasspdf, offshellZmasspdf, error!!! "
	    return  False
	#initial wmass_gen
	self.met_px[0] = self.met.Px(); self.met_py[0] = self.met.Py()
	it = 0
	genRandom = ROOT.TRandom3(0)
        #PUSample: 25.2, PU0: 14.8
        met_sigma = 25.2
        if not self.recobjetrescalec1pdf_flag:
		met_sigma = 0.0	
        #genRandom.SetSeed()
	self.dilepton_mass[0] = self.dilepton_p4.M()
        if self.dilepton_mass[0] < 60.0:
	    self.Znunu_mass[0] = 80.0
	else:
	    self.Znunu_mass[0] = 20.0


	#print "iterations ",self.iterations
	while (it < self.iterations ):
	    it += 1
	    self.initHMETree()
	    rand01 = genRandom.Uniform(0., 1.0)
	    step = genRandom.Uniform(-10.0, 10.0)
	    #print "self.iterations ",self.iterations," it ",it," dilepton_mass ",self.dilepton_mass[0]," self.Znunu_mass[0] ",self.Znunu_mass[0]," rand01 ",rand01," step ",step
	    if self.dilepton_mass[0] < 60.0:
	   	## z->mumu offshell, sample one onshell Z mass for z->nunu 
	        self.Znunu_mass[0] = self.getOnshellZMass(self.Znunu_mass[0], step, rand01) 
            else:
	   	## z->mumu onshell, sample one offshell Z mass for z->nunu 
	        self.Znunu_mass[0] = self.getOffshellZMass(self.Znunu_mass[0], step, rand01) 
	    self.hmass_gen[0] = genRandom.Gaus(125.03, 0.004)
	    #print "it ",it," self.eta_gen[0] ",self.eta_gen[0]," wmass_gen ",self.wmass_gen[0]
	    #update met 
	    if self.recobjetrescalec1pdf_flag:
		while not self.bjetsCorrection():
		    #print "fail to get bjetcorrection, try to get next one "
		    pass
		met_dpx = genRandom.Gaus(0.0, met_sigma)
		met_dpy = genRandom.Gaus(0.0, met_sigma)
		met_corr = self.met + ROOT.TVector2(met_dpx, met_dpy)+ self.metCorrection()
            else:
		met_corr = self.met
		self.b1rescalefactor[0] = 1.0
		self.b2rescalefactor[0] = 1.0
	    self.htoBB_p4 = self.b1jet_p4 * self.b1rescalefactor[0] + self.b2jet_p4 * self.b2rescalefactor[0]
	    self.metpx_corr[0]= met_corr.Px()
	    self.metpy_corr[0] = met_corr.Py()
    	    #print "Zll mass ",self.dilepton_p4.M()," px ",self.dilepton_p4.Px()," py ",self.dilepton_p4.Py()," pz ",self.dilepton_p4.Pz()
            #print "met_px ",self.met.Px()," met_py ",self.met.Py()," after correction px ",met_corr.Px()," py ",met_corr.Py()
	    #print "HtoBB mass ",self.htoBB_p4.M()
	    self.nsolutions[0] = 0
	    isolution = 0
	    solutions = [False, False]
	    #1. permutation 
	    #2. check nu_onshell_W pt
	    #3. solve the kinematics
	    #4. mark solution is ture if it is solved
	    #5. dump information into tree
	    while isolution < len(solutions):
		self.Znunu_pz[0], self.Znunu_energy[0] = self.ZtonunuPzAndE(self.dilepton_p4, met_corr, self.Znunu_mass[0], self.hmass_gen[0], isolution)
		if self.Znunu_energy[0]<0.0:	
		    isolution += 1
		    continue
	        self.Zll_p4 = self.dilepton_p4
	        self.Znunu_p4 = ROOT.TLorentzVector(met_corr.Px(), met_corr.Py(), self.Znunu_pz[0], self.Znunu_energy[0])
		self.htoZZ_p4 = self.Zll_p4 + self.Znunu_p4
		self.h2tohh_p4 = self.htoZZ_p4 + self.htoBB_p4

		if (fabs(self.htoZZ_p4.M() - self.hmass_gen[0])>1.0):
		    print "Error!! hmass_gen ", self.hmass_gen[0], " higgs mass from HME htoZZ_p4 ", self.htoZZ_p4.M()
		#print "get this h2tohh_mass ",self.h2tohh_p4.M()," iter ",it

		self.Zll_mass[0] = self.Zll_p4.M()
		self.Zll_eta[0] = self.Zll_p4.Eta()
		self.Zll_phi[0] = self.Zll_p4.Phi()
		self.Zll_pt[0] = self.Zll_p4.Pt()
		self.Zll_energy[0] = self.Zll_p4.Energy()
		self.Znunu_eta[0] = self.Znunu_p4.Eta()
		self.Znunu_phi[0] = self.Znunu_p4.Phi()
		self.Znunu_pt[0] = self.Znunu_p4.Pt()
		#self.Znunu_energy[0] = self.Znunu_p4.Energy()
		#self.Znunu_mass[0] = self.Znunu_p4.M()


    		self.htoZZ_eta[0] = self.htoZZ_p4.Eta()
    		self.htoZZ_phi[0] = self.htoZZ_p4.Phi()
    		self.htoZZ_pt[0] = self.htoZZ_p4.Pt()
    		self.htoZZ_energy[0] = self.htoZZ_p4.Energy()
    		self.htoZZ_mass[0] = self.htoZZ_p4.M()
    		self.htoBB_eta[0] = self.htoBB_p4.Eta()
    		self.htoBB_phi[0] = self.htoBB_p4.Phi()
    		self.htoBB_pt[0] = self.htoBB_p4.Pt()
    		self.htoBB_energy[0] = self.htoBB_p4.Energy()
    		self.htoBB_mass[0] = self.htoBB_p4.M()
    		self.h2tohh_pt[0] = self.h2tohh_p4.Pt()
    		self.h2tohh_energy[0] = self.h2tohh_p4.Energy()
    		self.h2tohh_mass[0] = self.h2tohh_p4.M()

		if (self.h2tohh_p4.Pt()/self.h2tohh_p4.E() <.00000001):
		    print "Strange case: h2tohh pt ", self.h2tohh_p4.Pt(), " energy ",self.h2tohh_p4.E()
		    self.h2tohh_eta[0] = 1000000.0; self.h2tohh_phi[0] = 0.0
		else:
		    self.h2tohh_eta[0] = self.h2tohh_p4.Eta(); self.h2tohh_phi[0] = self.h2tohh_p4.Phi()
		self.hme_h2Mass.Fill(self.h2tohh_mass[0], self.weight[0])
	        if self.onshellnuptpdf_flag:
		    self.weight1[0] = self.weight[0] * self.getWeightFromHist(self.onshellnuptpdf, self.nu_Zll_pt[0]) 
    		self.hme_h2MassWeight1.Fill(self.h2tohh_mass[0], self.weight1[0])
	 	isolution += 1 
	##### end of iteration



