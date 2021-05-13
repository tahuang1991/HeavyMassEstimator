import os
pwd = os.getcwd()+'/'
indir = 'testfiles/'
filelist = os.popen("find "+indir+" -type f |grep root").read().split("\n")
print filelist
for i in filelist:
  if "out" in i:
    continue
  o = "out_"+i[:-5]+".root"
  c_name = i[:-5]+"condor.sh"
  c_script = open(c_name, "write")
  c_script.write("#!/bin/bash\n")
  c_script.write("date\n")
  c_script.write("cd "+pwd+"\n")
  c_script.write("eval `scramv1 runtime -sh`\n")
  c_script.write("python runHME_HHbbWW_v2.py -i "+i+" -o "+o+"\n")
  c_script.write("echo 'done'\n")
  c_script.write("date\n")
  os.system("chmod 755 "+c_name)

sub_all = open("testfiles/submit_all_condor.sh", "write")
sub_all.write("""executable              = $(filename)
output                  = {pwd}/testfiles/output/$(filename).out
error                   = {pwd}/testfiles/error/$(filename).err
log                     = {pwd}/testfiles/log/$(filename).log
request_memory          = 4000M
+JobFlavour             = "workday"
universe                = vanilla
queue filename matching files *condor.sh""".format(pwd = pwd))
