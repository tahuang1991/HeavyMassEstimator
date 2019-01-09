import random
import ROOT
import os
ROOT.gROOT.SetBatch(1)
#gStyle from TStyle
ROOT.gStyle.SetStatW(0.17)
ROOT.gStyle.SetStatH(0.15)

ROOT.gStyle.SetOptStat(111110)

ROOT.gStyle.SetTitleStyle(0)
ROOT.gStyle.SetTitleAlign(13) ## coord in top left
ROOT.gStyle.SetTitleX(0.)
ROOT.gStyle.SetTitleY(1.)
ROOT.gStyle.SetTitleW(1)
#ROOT.gStyle.SetTitleTextColor(4)
ROOT.gStyle.SetTitleXSize(0.05)
ROOT.gStyle.SetTitleYSize(0.05)
ROOT.gStyle.SetTitleH(0.058)
ROOT.gStyle.SetTitleBorderSize(0)

ROOT.gStyle.SetPadLeftMargin(0.126)
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

plotdir = "TTbar_Gen/"


def hist1D(tree, todraw, xbins, cut, B):

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    b1 = ROOT.TH1F("%s"%B,"%s"%B,xBins,xminBin,xmaxBin)
    tree.Draw("%s>>%s"%(todraw,B),cut)
    ROOT.SetOwnership(b1, False)
    #b1.Print()
    return b1

#___________________________________________
def draw1D_v2(filelist,x_bins,x_title,cut,benchmarks, pic_name):
    

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("h2#rightarrow hh#rightarrow WWWW"+" "*24 + "14TeV")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    #b1.GetYaxis().SetRangeUser(0.0,10000)
    b1.SetStats(0)
     


    color = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+2, ROOT.kGreen+2, ROOT.kCyan]
    maker = [20,21,22,23,34]
    legend = ROOT.TLegend(0.75,0.6,0.86,0.94)
    legend.SetFillColor(ROOT.kWhite)
#    legend.SetFillStyle(0)
    legend.SetTextSize(0.05)
    legend.SetTextFont(62)
    hs1 = ROOT.THStack("hs1","%s distribution"%x_title)
    hs2 = ROOT.THStack("hs2","%s distribution"%x_title)
    hs3 = ROOT.THStack("hs3","%s distribution"%x_title)
    hs4 = ROOT.THStack("hs4","%s distribution"%x_title)
    hists_1 = []
    hists_2 = []
    hists_3 = []
    hists_4 = []
    for nfile in range(len(filelist)):
	hists_1.append(ROOT.TH1F("hist1_%d"%nfile,"hist1_%d"%nfile, 190, 200, 4000))
	hists_2.append(ROOT.TH1F("hist2_%d"%nfile,"hist2_%d"%nfile, 190, 200, 4000))
	hists_3.append(ROOT.TH1F("hist3_%d"%nfile,"hist3_%d"%nfile,30,0,300))
	hists_4.append(ROOT.TH1F("hist4_%d"%nfile,"hist4_%d"%nfile,30,0,900))
	
    for nfile in range(len(filelist)):
	
	rootfile = filelist[nfile]
	B = benchmarks[nfile]
    	f = ROOT.TFile(rootfile)
	lists = f.GetListOfKeys()
	for x in range(len(lists)):
		subkey = lists.At(x)
		obj = subkey.ReadObj()
		if obj.GetName()=="evtree":
			continue
		print " title ",obj.GetTitle()," Name ",obj.GetName()
		maxbin = obj.GetMaximumBin()
		
		hists_1[nfile].Fill(obj.GetXaxis().GetBinCenter(maxbin))
		hists_3[nfile].Fill(obj.GetBinContent(maxbin))
		hists_4[nfile].Fill(obj.Integral())
		obj.Scale(1.0/obj.Integral())
		hists_2[nfile].Add(obj.Rebin(20))
	hists_1[nfile].SetLineColor(color[nfile])
	hists_1[nfile].SetMarkerColor(color[nfile])
	hists_1[nfile].SetMarkerStyle(maker[nfile])
	
	hists_2[nfile].SetLineColor(color[nfile])
	hists_2[nfile].SetMarkerColor(color[nfile])
	hists_2[nfile].SetMarkerStyle(maker[nfile])

	hists_3[nfile].SetLineColor(color[nfile])
	hists_3[nfile].SetMarkerColor(color[nfile])
	hists_3[nfile].SetMarkerStyle(maker[nfile])

	hists_4[nfile].SetLineColor(color[nfile])
	hists_4[nfile].SetMarkerColor(color[nfile])
	hists_4[nfile].SetMarkerStyle(maker[nfile])
	
	if (nfile==len(filelist)-1):
		hists_1[nfile].Scale(1.0/5.0)
		hists_2[nfile].Scale(1.0/5.0)
	hists_3[nfile].Scale(1.0/hists_3[nfile].Integral())
	hists_4[nfile].Scale(1.0/hists_4[nfile].Integral())

	hs1.Add(hists_1[nfile])
	hs2.Add(hists_2[nfile])
	hs3.Add(hists_3[nfile])
	hs4.Add(hists_4[nfile])
	legend.AddEntry(hists_1[nfile], "%s"%B, "pl")
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
  
    c1.Divide(2,2)
    c1.cd(1)
    hs1.Draw("nostack+p")
    hs1.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
    hs1.GetHistogram().GetXaxis().SetRangeUser(xminBin, xmaxBin)
    legend.Draw("same")
    tex1 = ROOT.TLatex(0.25,.50,"Most probable mass")
    tex1.SetTextSize(0.05)
    tex1.SetTextFont(62)
    tex1.SetNDC()
    tex1.Draw("same")

    c1.cd(2)
    hs2.Draw("nostack+p")
    hs2.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
    hs2.GetHistogram().GetXaxis().SetRangeUser(xminBin, xmaxBin)
    legend.Draw("same")
    tex2 = ROOT.TLatex(0.25,.50,"add all survival solution(normalized in each event)")
    tex2.SetTextSize(0.05)
    tex2.SetTextFont(62)
    tex2.SetNDC()
    tex2.Draw("same")

    c1.cd(3)
    hs3.Draw("nostack+p")
    hs3.SetTitle("bincontent of maximum bin, normalized  distribution")
    hs3.GetHistogram().GetXaxis().SetTitle("maximum bincontent")
    legend.Draw("same")
    tex3 = ROOT.TLatex(0.25,.50,"maximum  bin content from MMC")
    tex3.SetTextSize(0.05)
    tex3.SetTextFont(62)
    tex3.SetNDC()
    tex3.Draw("same")

    c1.cd(4)
    hs4.SetTitle("total number of survival solutions, normalized distribution")
    hs4.Draw("nostack+p")
    hs4.GetHistogram().GetXaxis().SetTitle("total number of survival solutions")
    legend.Draw("same")
    tex4 = ROOT.TLatex(0.25,.50,"total number of survival solutions")
    tex4.SetTextSize(0.05)
    tex4.SetTextFont(62)
    tex4.SetNDC()
    tex4.Draw("same")

    c1.cd()

    c1.SaveAs("Hhh_PDFvalidation_%s_combined.png"%pic_name)
    
	
		
#___________________________________________
def draw1D(filelist,dirlist,todraw,x_bins,x_title,cut,benchmarks, pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("h2#rightarrow hh#rightarrow WWWW"+" "*24 + "14TeV")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    b1.GetYaxis().SetRangeUser(0.0,10000)
    b1.SetStats(0)
    #b1.SetLogx()
    b1.Draw()
    #c1.SetLogx()
 
    color = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+2, ROOT.kGreen+2, ROOT.kCyan]
    maker = [20,21,22,23,34]
    legend = ROOT.TLegend(0.65,0.65,0.8,0.94)
    legend.SetFillColor(ROOT.kWhite)
#    legend.SetFillStyle(0)
    legend.SetTextSize(0.05)
    legend.SetTextFont(62)
    #hs = ROOT.THStack("hs","%s distribution"%x_title)
    hs = ROOT.THStack("hs"," ")
    hists = []
    for nfile in range(len(filelist)):
        chain = ROOT.TChain(dirlist[nfile])	
	filedir = filelist[nfile]
	B = benchmarks[nfile]
    	if os.path.isdir(filedir):
          ls = os.listdir(filedir)
          for x in ls:
                x = filedir[:]+x
                if os.path.isdir(x):
                        continue
                chain.Add(x)
    	elif os.path.isfile(filedir):
          chain.Add(filedir)
    	else:
          print " it is not file or dir ", filedir

	#chain.Draw("%s>>%s(%d,%f,%f)"%(todraw,B, xBins, xminBin, xmaxBin ), cut)
	hist = hist1D(chain, todraw, x_bins, cut, B)
	#hs.Add(ROOT.TH1F(ROOT.gDirectory.Get(B)))
        #legend.AddEntry(ROOT.TH1F(ROOT.gDirectory.Get(B)),"%s"%B,"pl") 
	#hists.append(ROOT.TH1F(ROOT.gDirectory.Get(B)).Clone(B))
	hist.SetLineColor(color[nfile])
	hist.SetLineWidth(3)
	hist.SetMarkerColor(color[nfile])
	hist.SetMarkerStyle(maker[nfile])
	hist.Scale(1.0/hist.GetEntries())
        hs.Add(hist)
        
	legend.AddEntry(hist, "%s"%B, "l")

	#htmp = hist1D(t, todraw, x_bins, cut, B)
	#htmp.Print("ALL")
	hists.append(hist)
	print "nfile ",nfile," B ",B
    print "hists ",hists	
 
    hs.Draw("nostack")
    hs.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
    hs.GetHistogram().GetYaxis().SetTitle("Normalized to unity")
    legend.Draw("same")
    c1.SaveAs(plotdir+"Hhh_PDFvalidation_hasdRljet_%s.pdf"%pic_name)
    c1.SaveAs(plotdir+"Hhh_PDFvalidation_hasdRljet_%s.png"%pic_name)
    c1.SaveAs(plotdir+"Hhh_PDFvalidation_hasdRljet_%s.C"%pic_name)




filedir = "/fdata/hepx/store/user/lpernie/Hhh/"
filelist = [filedir+"delphes_B3_1M_PU40ALL_13May.root",filedir+"delphes_B6_1M_PU40ALL_13May.root",filedir+"delphes_B9_1M_PU40ALL_13May.root",filedir+"delphes_B12_1M_PU40ALL_13May.root"]
#filelist = ["DiHiggs_WWbb_1M_NewB3_allReco_25_MMC1M_isomu_MVA_PU40_0415_combined.root","DiHiggs_WWbb_1M_NewB6_allReco_25_MMC1M_isomu_MVA_PU40_0509_combined.root","DiHiggs_WWbb_1M_NewB9_allReco_25_MMC1M_isomu_MVA_PU40_0509_combined.root"]
filedir = "/fdata/hepx/store/user/taohuang/Hhh/combined_samples/"
#filelist = [filedir+"DiHiggs_WWbb_1M_NewB3_allReco_25_MMC1M_isomu_MVA_PU40_0725_combined.root",filedir+"DiHiggs_WWbb_1M_NewB6_allReco_25_MMC1M_isomu_MVA_PU40_0725_combined.root",filedir+"DiHiggs_WWbb_1M_NewB9_allReco_25_MMC1M_isomu_MVA_PU40_0725_combined.root",filedir+"DiHiggs_WWbb_1M_NewB12_allReco_25_MMC1M_isomu_MVA_PU40_0725_combined.root"]#,filedir+"TTbar_Wtomu_4M_allReco_Updatebtag_25_MVA_PU40_0725_combined.root"]
filelist = [filedir+"DiHiggs_WWbb_1M_NewB3_allReco_simulation_isomu_MVA_PU40_0824_combined.root", filedir+"DiHiggs_WWbb_1M_NewB6_allReco_simulation_isomu_MVA_PU40_0824_combined.root",filedir+"DiHiggs_WWbb_1M_NewB9_allReco_simulation_isomu_MVA_PU40_0824_combined.root","/fdata/hepx/store/user/taohuang/Hhh/Delphes_ttbar_4M_Wtomu_allReco_Updatebtag_simulation_PU40_0825/"]
benchmarks = ["B3","B6","B9","TTbar"]	
filelist = ["out_ana_0000.root","./Delphes_ttbar/"]
treelist = ["DiHiggsWWBBAna/evtree","evtree"]
benchmarks = ["CMSSW","Delphes"]
cut = "h2tohh && htoWW && mu1_pt>10 && mu2_pt>10 && fabs(mu1_eta)<2.4 && fabs(mu2_eta)<2.4 && genmet>20"

pic_name = "MMC_h2mass"	
x_title = "Most Probable mass"
x_bins = "(60,200,1400)"
#draw1D_v2(filelist,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "EtaonshellW"	
todraw = "w1_eta*(w1_mass>w2_mass)+w2_eta*(w1_mass<w2_mass)"
x_title = "#eta of onshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PhionshellW"	
todraw = "w1_phi*(w1_mass>w2_mass)+w2_phi*(w1_mass<w2_mass)"
x_title = "#Phi of onshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PtonshellW"	
todraw = "w1_pt*(w1_mass>w2_mass)+w2_pt*(w1_mass<w2_mass)"
x_title = "p_{T} of onshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "EtanuonshellW"	
todraw = "nu1_eta*(w1_mass>w2_mass)+nu2_eta*(w1_mass<w2_mass)"
x_title = "#eta of neutrinos from onshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PhinuonshellW"	
todraw = "nu1_phi*(w1_mass>w2_mass)+nu2_phi*(w1_mass<w2_mass)"
x_title = "#phi of neutrinos from onshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PtnuonshellW"	
todraw = "nu1_pt*(w1_mass>w2_mass)+nu2_pt*(w1_mass<w2_mass)"
x_title = "p_{T} of neutrinos from onshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "EtaoffshellW"	
todraw = "w1_eta*(w1_mass<w2_mass)+w2_eta*(w1_mass>w2_mass)"
x_title = "#eta of offshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PhioffshellW"	
todraw = "w1_phi*(w1_mass<w2_mass)+w2_phi*(w1_mass>w2_mass)"
x_title = "#Phi of offshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PtoffshellW"	
todraw = "w1_pt*(w1_mass<w2_mass)+w2_pt*(w1_mass>w2_mass)"
x_title = "p_{T} of offshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "EtanuoffshellW"	
todraw = "nu1_eta*(w1_mass<w2_mass)+nu2_eta*(w1_mass>w2_mass)"
x_title = "#eta of neutrinos from offshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PhinuoffshellW"	
todraw = "nu1_phi*(w1_mass<w2_mass)+nu2_phi*(w1_mass>w2_mass)"
x_title = "#phi of neutrinos from offshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "PtnuoffshellW"	
todraw = "nu1_pt*(w1_mass<w2_mass)+nu2_pt*(w1_mass>w2_mass)"
x_title = "p_{T} of neutrinos from offshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "leadingbjetratio"
todraw = "(b1_pt/b1jet_pt)*(b1jet_pt>b2jet_pt)+(b2_pt/b2jet_pt)*(b1jet_pt<b2jet_pt)"
x_title = "#frac{p_{T}(b parton)}{p_{T}(bjet)}, leading bjet"
x_bins = "(600,0.0,3.0)"
cut = "h2tohh && hasRECOjet1 && hasRECOjet2 && dR_b1jet<0.4 && dR_b2jet<0.4 && hastwomuons"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "MT2_reco"
todraw = "MT2_reco"
x_bins = "(80,0.0,400)"
x_title = "M_{T2} [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "dR_l1l2"
todraw = "dR_l1l2"
x_bins = "(50,0.0,5)"
x_title = "#Delta R(l,l)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "dR_l1l2gen"
todraw = "dR_genl1l2"
x_title = "#Delta R(l,l), gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "dR_b1b2"
todraw = "dR_b1b2"
x_bins = "(50,0.0,5)"
x_title = "#Delta R(j,j)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "dR_b1b2gen"
todraw = "dR_genb1b2"
x_title = "#Delta R(j,j), gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "dR_l1l2b1b2"
todraw = "dR_l1l2b1b2"
x_bins = "(50,0.0,5)"
x_title = "#Delta R(ll,jj)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "dR_l1l2b1b2gen"
todraw = "dR_genl1l2b1b2"
x_title = "#Delta R(ll,jj), gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "dphi_l1l2b1b2"
todraw = "dphi_l1l2b1b2"
x_bins = "(50,0.0,3.5)"
x_title = "#Delta #phi(ll,jj)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "dphi_l1l2b1b2gen"
todraw = "dphi_genl1l2b1b2"
x_title = "#Delta #phi(ll,jj), gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "dphi_llmet"
todraw = "dphi_llmet"
x_bins = "(50,0.0,3.5)"
x_title = "#Delta #phi(ll, #slash{E}_{T})"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "dphi_llmetgen"
todraw = "dphi_genllmet"
x_title = "#Delta #phi(ll, #slash{E}_{T}), gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	

pic_name = "dR_minbl"
todraw = "dR_minbl"
x_bins = "(50,0.0,5)"
x_title = "min #Delta R(l,j)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"#hasdRljet
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "dR_minblgen"
todraw = "dR_genminbl"
cut = "1"
x_title = "min #Delta R(l,j), gen level"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "mass_l1l2"
todraw = "mass_l1l2"
x_bins = "(50,.0,400)"
x_title = "M(ll) [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "mass_l1l2gen"
todraw = "mass_genl1l2"
cut = "1"
x_title = "M(ll), gen level"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "mass_b1b2"
todraw = "mass_b1b2"
x_bins = "(80,0.0,400)"
x_title = "M(jj) [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "mass_b1b2gen"
todraw = "mass_genb1b2"
x_title = "M(jj), gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


pic_name = "mass_trans"
todraw = "mass_trans"
x_bins = "(60,0.0,300)"
x_title = "M_{trans} [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)	
pic_name = "mass_transgen"
todraw = "mass_gentrans"
x_title = "M_{trans}, gen level"
cut = "1"
draw1D(filelist, treelist,todraw,x_bins,x_title,cut,benchmarks, pic_name)	


