import os
import sys 
sys.argv.append( '-b' )

#espresso     = 20 minutes
#microcentury = 1 hour
#longlunch    = 2 hours
#workday      = 8 hours
#tomorrow     = 1 day
#testmatch    = 3 days
#nextweek     = 1 week
#masslist = [1000, 1250, 1500, 1750, 2000, 2500, 3000, 700, 800, 900, 750, 850, 650, 600]
masslist = [250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700,750, 800, 850, 900, 1000] 
#masslist = [2000]
infilelist = {"RadionM400":"/eos/user/t/tahuang/HME_Systematic_202301/out_radion_2L_2016_m400.root",
        "TTbar" : "/eos/user/t/tahuang/HME_Systematic_202301/out_ttbar_2L_2016.root"      
        }

iterations = 10000 ## default=10k, opt: 10k, 50k, 100k
nEv_per_job = 4000
#njobs = 10
bashscript = "runHME_sys.sh"

def getEntries(filename, treename="Events"):
    import ROOT
    tfile = ROOT.TFile(filename)
    tree = tfile.Get(treename)
    return tree.GetEntries()

#def generate_run_HME(masslist, workdir, inputdir, outdir, njobs):
def generate_run_HME(infilelist, workdir, outdir, nEv_per_job):
    if not os.path.exists(workdir):
        os.system("mkdir "+workdir)
        os.system("mkdir "+workdir+"/output")
        os.system("mkdir "+workdir+"/error")
        os.system("mkdir "+workdir+"/log")
    if not os.path.exists(outdir):
        os.system("mkdir "+outdir)
    os.system("cp "+bashscript + " " +workdir)
    fname_all = workdir+"condor_RunHME_allsys.sh"
    script_all = open(fname_all, "write")
    script_all.write("#!/bin/bash\n")
    #script_all.write("eval `scramv1 runtime -sh`\n")
    for infilename in infilelist.keys():
       inputfile = infilelist[infilename]
       filename = inputfile.split("/")[-1][:-5]
       nTotal = getEntries(inputfile,"Double_Tree")
       #nEv_per_job = nTotal/njobs + 1
       njobs = nTotal/nEv_per_job + 1
       for ijob in range(njobs):
           nStart = nEv_per_job*ijob
           nEnd = nEv_per_job*(ijob+1)
           if nEnd > nTotal: nEnd = nTotal
           outputfile = outdir           +filename+"HME_ijob%d.root"%(ijob)
           batchfname = workdir + "Batch"+filename+"HME_ijob%d.cmd"%(ijob)
           print("filename ", filename, " totalevt ", nTotal," nstart ", nStart, " nEnd ", nEnd, " out ", outputfile)
           script_all.write("condor_submit "+batchfname + "\n")
           batchscript = open(batchfname, "write")
           batchscript.write("""universe              = vanilla 
executable            = {script}
arguments             = {arg1} {arg2} {arg3} {arg4} {arg5}
output                = output/{suffix}.$(ClusterId).$(ProcId).out
error                 = error/{suffix}.$(ClusterId).$(ProcId).err
log                   = log/{suffix}.$(ClusterId).log
request_memory        = 4000M                                                                                                                        
+JobFlavour           = "testmatch"
Notification          = Complete
notify_user           = taohuang@email.tamu.edu
queue""".format(script = bashscript, suffix = filename+"ijob%d"%(ijob), arg1=inputfile, arg2=outputfile, arg3=nStart, arg4=nEnd, arg5=iterations))
    os.system("chmod 775 "+fname_all)

newfolder = "HME_METUnclusterScaleSys_20230123/"
outdir = "/eos/user/t/tahuang/HME_Systematic_202301/"+newfolder
workdir = os.path.join(os.getcwd(), newfolder)
print("workdir ", workdir)


generate_run_HME(infilelist, workdir, outdir, nEv_per_job)
