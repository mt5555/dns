#!/usr/bin/env python
from string import *
import os, commands, getopt, sys
####################################################################################
#
#  ./cronsub.py [submit] 
#
# automatically resubmit LSF jobs in ~/cronlsf/*.job
#
# running without any options will print all messages, but not actually
# submit any jobs.  Good for testing.
#
# To use:
#
#  add script to ~/cronlsf
#  make sure script has a unique job name, with a line like:  #BSUB -J iso12
#  create a file iso12.resub with a single line (ASCII) containing
#  the number of runs to submit.
#
#  The .resub file number will be decreased each time a job is submitted
#  until it gets to zero.  
#
#  For all jobs not submitted, the reason will be printed to stdout
#  For jobs submitted, the command will be echoed to stderr
#  So running every 5min under cron, use:
#  
#  ./cronsub.py submit >> /dev/null
#
#  And it will only generate mail every time a job is submitted to LSF
#
#
####################################################################################

submit=0
if (len(sys.argv)== 2):
    if (sys.argv[1]=="submit"):
        submit=1

path="/home/mt/lanl/dns/benchmark/"
user="taylorm"


# get a list of all jobnames queued in LSF
cmd = "bjobs -u all"
(status,out) = commands.getstatusoutput(cmd)
if (status!=0):
    print 'Error getting list of LSF jobs'
    sys.exit(1)

#parse out to get jobname_running
jobname_running=[]
jobc=-1;
vout=split(out,"\n")

if (len(vout)<2):
    print 'Error parsing bjobs output: need at least 2 lines'
    sys.exit(1)

for line in vout:
    sline=split(line," ")
    if (jobc==-1) & (len(sline)>=7) & (sline[6]=="JOB_NAME"):
        #find column where job starts:
        jobc=find(line,"JOB_NAME")

    if (jobc<0):
       print 'Error parsing bjobs output for JOB_NAME'
       sys.exit(1)
    
        
    if (len(sline)>=7) & (sline[1]==user) & (jobc>=0):
        # get everything in JOB_NAME column and beyond:
        out=line[jobc:-1]
        out=split(out," ")
        if (len(out)==0):
            print 'Error parsing bjobs output for jobname'
            sys.exit(1)
        jobname_running.append(out[0]);
        



# are we on IRIX or OSF1?
(status,out) = commands.getstatusoutput("uname")
if (status==0):
    if (out=="IRIX64"):
        bsub="bsub"
    else:
        bsub="idbsub"
else:
    print 'Error getting OS type'
    sys.exit(1)




# get a list of all the .job files in ~/cronlsf:
cmd = "ls "+path+"*.job"
(status,out) = commands.getstatusoutput(cmd)
if (status!=0):
    print 'Error getting list of scripts'
    sys.exit(1)
vjobscript=split(out,"\n")


for jobscript in vjobscript:

    print jobscript,
    try:

        # parse file for the job name
        #print "LSF script: ",jobscript
        jobname="";
        fid=open(jobscript)
        line=rstrip(fid.readline())
        while line:
            out=split(line," ")
            if (len(out)>=3) & (out[0]=="#BSUB"):
                if (out[1]=="-J"):
                    jobname=out[2]
                    break
            line=rstrip(fid.readline())


        if len(jobname)==0:
            break

        # check if it is in que
        i=jobname in jobname_running
        if i:
            raise NameError,"job already running: "+jobname

        
        # check resub count in file jobname.resub
        jobresub=path+jobname+".resub"
        fid=open(jobresub,'r')
        line=fid.readline()
        fid.close()
        fvalue=atoi(line)
        if (fvalue<=0):
            raise NameError,"job counter reached 0"
        fvalue=fvalue-1
        fid=open(jobresub,'w')
        fid.write(str(fvalue)+"\n")
        fid.close()

            
    except IOError,e:
        print 'FILE_ACCESS_ERROR: ',jobname

    except ValueError,e:
        print 'VALUE_ERROR: ',jobname

    except NameError,e:
        print 'NAME_ERROR: ',e

    else:
        # submit job
        jobcommand = bsub + " < " + jobscript
        out= "resub=" + str(fvalue)+" LSF job: "+jobcommand+"\n"
        sys.stderr.write(out)
        if (submit):
            os.system(jobcommand)
        else:
            sys.stderr.write("(testing mode - no jobs submitted)\n")



   

