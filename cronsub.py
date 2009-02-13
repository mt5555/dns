#!/usr/bin/env python
from string import *
import os, commands, getopt, sys, exceptions
##############################################################################
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
#  create a file iso12.job.resub with a single line (ASCII) containing
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
#############################################################################
class ContinueEx(exceptions.Exception):
    def __init__(self,args=None):
        self.args=args



#testing mode, unless arg1= "submit"
submit=0
if (len(sys.argv)== 2):
    if (sys.argv[1]=="submit"):
        submit=1


#user="taylorm"
user="mataylo"


#set below for value of path
#path="cronlsf/"    #requires that we run from home directory
                   # which is what cron does. otherwise hard code
                   # a full path

# are we on IRIX or OSF1?
cmdstat = "bjobs -u all"
jobopt = " < "
(status,out) = commands.getstatusoutput("/bin/hostname")
if (status==0):
    if (out=="qbfe1"):
        bsub="bsub"
        path="cronqb/"
    elif (out=="qfe1"):
        bsub="bsub"
        path="cronqa/"
    elif (out=="qfe2"):
        bsub="bsub"
        path="cronqb/"
    elif (out=="qscfe1"):
        bsub="bsub"
        path="cronqsc/"
    elif (out=="blogin2.sandia.gov"):
        bsub="/apps/torque/bin/qsub"
        jobopt = " " 
        cmdstat = "/apps/torque/bin/qstat"
        path="crontbird/"
    else:	
        print 'Error getting hostname'
        sys.exit(1)
else:
    print 'Error getting OS type'
    sys.exit(1)





# get a list of all jobnames queued in LSF

(status,out) = commands.getstatusoutput(cmdstat)
if (status!=0) & (status!=255):
    print 'Error getting list of queued jobs'
    sys.exit(1)

# sometimes there are NO jobs in the system:
if (status==255) & (find(out,"No unfinished job found")>=0):
    # in this case, there were no jobs in the system,
    # so dont parse to find our running jobs
    print "jobs: que is empty?"
    print "jobs: ",out
else:    
    #parse out to get jobname_running
    jobname_running=[]
    jobc=-1;
    vout=split(out,"\n")

    if (len(vout)<2):
        print 'Error parsing jobs output: need at least 2 lines'
        sys.exit(1)

    for line in vout:
        sline=split(line)
        if (jobc==-1) & (len(sline)>=7):
            if (sline[6]=="JOB_NAME"):     # LSF output
                #find column where job starts
                jobc=find(line,"JOB_NAME")
            if (sline[2]=="Name"):         # PBS output 
                #find column where job starts
                jobc=find(line,"Name")

        if (jobc<0):
            print 'Error parsing bjobs output for job name'
            sys.exit(1)

        if (len(sline)>=3) & (jobc>=0):
            if (sline[1]==user):
                # LSF get everything in JOB_NAME column and beyond:
                out=line[jobc:-1]
                out=split(out," ")
                if (len(out)==0):
                    print 'Error parsing jobs output for jobname'
                    sys.exit(1)
                jobname_running.append(out[0]);
            if (sline[2]==user):  
                # PBS get everything in JOB_NAME column and beyond:
                out=line[jobc:-1]
                out=split(out," ")
                if (len(out)==0):
                    print 'Error parsing jobs output for jobname'
                    sys.exit(1)
                jobname_running.append(out[0]);


    print "current queued jobs for user ",user
    if (len(jobname_running)==0):
        print "<none>"		
    for out in jobname_running:
        print out


print ' '


# get a list of all the .job files in ~/cronlsf:
cmd = "ls "+path+"*.job"
(status,out) = commands.getstatusoutput(cmd)
if (status!=0):
    print 'Error: didn''t find any que scripts using: ',cmd
    sys.exit(1)
vjobscript=split(out,"\n")


for jobscript in vjobscript:

    print ' ' 
    print 'script: ',jobscript
    try:

        # parse file for the job name
        #print "LSF script: ",jobscript
        jobname="";
        jobcpus="";
        fid=open(jobscript)
        line=fid.readline()
        while line:
            out=split(rstrip(line)," ")
            if (len(out)>=3) & (out[0]=="#BSUB"):
                if (out[1]=="-J"):
                    jobname=out[2]
            if (len(out)>=3) & (out[0]=="#PBS"):
                if (out[1]=="-N"):
                    jobname=out[2]
            line=fid.readline()


        if len(jobname)==0:
            raise ContinueEx,"no BSUB -J option: "+jobname
        if len(jobname)>10:
            raise ContinueEx,"BSUB -J <jobname>: jobname too long! "+jobname

        # for idbsub, #BSUB -n 64 line is ignored, we have to add this
        # to the idbsub line
        #if len(jobcpus)>0:
        #    jobcpus=" -n "+jobcpus+" "

        print 'que name: ',jobname

        # check if it is in que
        i=jobname in jobname_running
        if i:
            raise ContinueEx,"job already running."

        
        # check resub count in file jobname.resub
        jobresub=path+jobname+".resub"
        fid=open(jobresub,'r')
        line=fid.readline()
        fid.close()
        fvalue=atoi(line)
        if (fvalue<=0):
            raise ContinueEx,"job counter reached 0"
        fvalue=fvalue-1
        fid=open(jobresub,'w')
        fid.write(str(fvalue)+"\n")
        fid.close()

            
    except IOError,e:
        print 'FILE_ACCESS_ERROR: (resub file?)'

    except ValueError,e:
        print 'VALUE_ERROR: '

    except ContinueEx,e:
        print 'skipping script: ',e

    else:
        # submit job
        jobcommand = bsub + jobcpus + jobopt + jobscript
        print "resub=" + str(fvalue)+" que job: "+jobcommand
        if (submit):
            os.system(jobcommand)
        else:
            print "(testing mode - no jobs submitted)"



   

