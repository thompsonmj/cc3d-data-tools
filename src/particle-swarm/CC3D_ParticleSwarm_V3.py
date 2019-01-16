# Particle Swarm Optimization - PSO for CC3D and using a compute cluster to calculate the particles in parallel
# Python 2.7
#
# J Sluka
# Indiana University
# Feb. 2018  Aug., Sept. 2018
# based on C:\Users\James Sluka\Desktop\Small Code Projects\ParticleSwarmPython\ParticleSwarm_V2.py
#
# on tatooine (linux), from the directory that contains CC3D_ParticleSwarm_V3.py, run with:
#         /home/jsluka/Parameter377/CC3D_3.7.7_RHEL/Python27/bin/python CC3D_ParticleSwarm_V3.py
#         /home/shared_sw/CC3D_v378/Python27/bin/python                 CC3D_ParticleSwarm_V3.py
# or
#   nohup /home/jsluka/Parameter377/CC3D_3.7.7_RHEL/Python27/bin/python CC3D_ParticleSwarm_V3.py &
#   nohup /home/shared_sw/CC3D_v378/Python27/bin/python                 CC3D_ParticleSwarm_V3.py &
# 

import random
import math
import os
import re   # regular expressions
import datetime
import timeit
import sys
sys.path.insert(0,'/home/jsluka/Parameter377/CC3D_3.7.7_RHEL/Python27/site-packages/slurmpy') # since slurmpy is installed locally
from slurmpy import Slurm

swarmStartTime=datetime.datetime.now()

# PSO parameters and options
w  = 0.73 # 0.73 Inertia weight 
c1 = 1.50 # 1.50 Scaling co-efficient on the social component
c2 = 1.50 # 1.50 Scaling co-efficient on the cognitive component
maxVelocity= 1.  # scales the parameters min -- max range. 1 uses the whole range, 0.5 uses half the range
iterations = 3   # maximum number of swarm iterations 
swarmSize  = 8   # typically at least the number of dimensions, usually more 
numSwarms  = 4   # total number of swarms, each is of swarmSize 
numDuplRuns= 3   # run this many replicate CC3D runs for each parameter set
lastIterWchange = 0

# CC3D paths and options
cc3dPath   = '/u/jsluka/Parameter377/CC3D_3.7.7_RHEL/runScript.sh'  # full path to the CC3D run script
Quality_data_file_name = '' # full path plus file name "Quality_data.txt"
projectPath='/u/jsluka/swarm/'  # the command above isn't getting the path!
print "projectPath: ",projectPath
cc3dJobName = 'NewSimulation.cc3d'  # *.cc3d file name, could be obtained using the projectPath
cc3dParameterFile = 'NewSimulation_parameters.py'  # *_parameters.py, this file lives in the CC3D project's "Simulation" folder
#theUserID = 'jsluka'  # not needed since this script automatically gets the user's ID

####################################################################################
# Submit the set of simulations via slurm (python Slurm)
#
# Returns a vector of result values of length numSwarms*swarmSize
# that can be iterated over in the proper order with;
#        for iS in range(0,numSwarms):
#            for iP in range(0,swarmSize):
#
def submitJobs(iter,dupRun):
    import sys
    import time
    import subprocess

    # get the user's ID, needed to monitor the slurm jobs
    returned_output = subprocess.check_output(['id -un'],shell=True)  # get the users ID
    theUserID=returned_output.decode("utf-8")  # using decode() function to convert byte string to string
    theUserID=theUserID.strip()  # remove white spaces
    print "\n            ==user id is==",theUserID,"=====\n\n"
    
    # submit the jobs via sbatch, iterate over all particles in all swarms
    for iS in range(0,numSwarms):
        for iP in range(0,swarmSize):    
            jobName=dirArray[iS][iP]  # the directory name of the cc3d job, e.g., "S00004_P00009"
            print "    Swarm ",iS+1," of",numSwarms,"   Particle ",iP+1," of ",swarmSize,"     "
            s = Slurm(jobName)

            slurmCom ="#! /usr/bin/sh\n"
            slurmCom+="#SBATCH -o "+jobName+"stdout.txt\n"
            slurmCom+="#SBATCH -t 2:00:00  # time requested in hour:minute:second\n"
            slurmCom+="#SBATCH -J S"+str(iS)+"P"+str(iP)+"  # overwrite what slurm.py wants to use\n\n"
            slurmCom+="srun -o "+jobName+"cc3d_console_log.txt "
            slurmCom+=cc3dPath+" "
            slurmCom+="-i "+jobName+"../NewSimulation.cc3d "
            slurmCom+="-f 0  "
            slurmCom+="-o "+jobName+"\n"   
            a=s.run(slurmCom)
            
            #print "submitted JobID=",a,"job name=",jobName
            time.sleep(0.1)  # (no longer needed) 2,3 otherwise the cc3d job dir's overlap because the folder names are 
                             # based on time to the second omitting the default output from cc3d might avoid this problem
    print "Submitted ",numSwarms*swarmSize," jobs. ************************** iter:",(iter+1), \
           " dupRun:",(dupRun+1),"****************************"

    # jobs all submitted, now wait for them to all end
    print "waiting..."    
    waiting = True  # an ugly way to do it but it works
    while waiting:
        time.sleep(5) 
        # returns output as byte string
        getCommand='squeue -hu '+theUserID   # -h omits the header
        squeue_output = subprocess.check_output([getCommand],shell=True)  # only user defined by theUserID
        squeueData = squeue_output.decode("utf-8")
        jobList= squeueData.splitlines()
        print '                  waiting, job count:',len(jobList)
        if len(jobList) < 1:
            waiting = False
    
    # parse the job output files and extract the Quality value
    # updated to include the run time in the quality file
    print "get results ..."
    results=[[999999]*swarmSize for _ in range(numSwarms)]
    timings=[[999999]*swarmSize for _ in range(numSwarms)]  # NEW
    for iS in range(numSwarms):
        for iP in range(swarmSize):
            path=dirArray[iS][iP]+'Quality_data.txt'  # the dirArray goes from root down to Simulation dir
            oFile = open(path,'r')
            aLine=oFile.readline()
            (qual,runTime) = aLine.split()  # comma seperated values,  NEW
            results[iS][iP]=(float(qual))  # NEW
            timings[iS][iP]=(float(runTime))  # NEW
            oFile.close()
        
    print 'Results:\n',results
    print 'Timings:\n',timings  # NEW
    print "Swarm best Qualities:"
    for iS in range(numSwarms):   # numSwarms
        print "%9.4f" % (min(results[iS])),
    print "\n"
    return results,timings

####################################################################################
# The parameters (x's) to be explored are stored in a comma delimited file of format;
#       parameterName,minValue,maxValue
# "#" can be used to comment out entire lines or ends of lines.
# Blank lines are ignored
# Default filename (assumed to be in the same directory as this script) is "swarmParameters.txt"
# Returns a list of parameter names, and min and max values, e.g.,
#       paramData[index]{'varName'|'min'|'max'} = number or txt

def loadParameters():
    print "\n    Loading parameters to be varied (x's) ...\n"
    try:
        paramFile = open('swarmParameters.dat','r')
    except OSError as e:
        print "Could not find the swarmParameters.dat file. ",e
        sys.exit(2)
    
    pData=[]
    pCount=0
    for cnt, aline in enumerate(paramFile):
        comments = ""
        sO=re.search(r'(.*?)#+(.*)', aline)
        if sO:
            comments=sO.group(2)
            aline=sO.group(1)
           
        aline=re.sub(r'\s'  ,"", aline)  # remove all whitespace
        if len(aline)>0:
            (paramN,minV,maxV) = aline.split(",")  # comma seperated values
            print "===",paramN,"===",minV,"===",maxV,"==="
            if float(minV) > float(maxV):  # do a simple sanity check, the min must be less than the max!
                print "\n\n        Problem with the parameters in swarmParameters.dat"
                print     "        the minimum values is not less than the maximum value!"
                print "                ",aline,"\n\n"
                sys.exit(2)
            pData.append({})  # append a new dictionary to the end of the list
            pData[pCount]['varName'] = paramN
            pData[pCount]['min'] = float(minV)
            pData[pCount]['max'] = float(maxV)
            pData[pCount]['comments'] = comments
            pCount += 1
    paramFile.close()
    return(pData)   # returns a list of dictionaries

####################################################################################
# Create the working directories
# The working directories are all created in the same directory as the .cc3d file.
# There is a folder for each particle in each swarm. The folders are named;
#       S00001_P00001
# for  particle 1 of swarm 1. Note that the indexes are zero referenced so the
# first directory is named "S00000_P00000" and for 5 swarms of 10 particles each the
# last directory is "S00004_P00009".
#
# For each working direcdtoy the .cc3d file and the entire "Simulation" folder is
# copied from the top level directory
# Returns and matrix indexed by swarm number, then particle number, with value of the directory's path

from distutils.dir_util import copy_tree
import shutil
def makeWorkingDirs():
    print "\n    Creating directories and copying cc3d files ...\n"
    dirArray =[[""]*swarmSize for _ in range(numSwarms)]
    for iS in range(0,numSwarms):
        for iP in range(0,swarmSize):
            newDirName = "S{:0>4d}_P{:0>4d}".format(iS,iP)
            path=projectPath+newDirName
            #print "Trying to create directory:",path
            try:
                os.mkdir( path, 0755 );
            except OSError as e:
                print "Could not make the output directory. ",e

            print "Created directory:",path
            # The approach below avoids copying the newly created directories into themselves
            # copy the .cc3d file
            fFrom = cc3dJobName
            fTo   = projectPath+"/"+newDirName+"/"+cc3dJobName
            shutil.copy(fFrom,fTo) # source, destination cc3dJobName

            # copy the entire folder "Simulation"
            fFrom = projectPath+'Simulation/'
            fTo   = projectPath+newDirName+'/Simulation/'
            copy_tree(fFrom, fTo)
            dirArray[iS][iP] = fTo

    return(dirArray)

# This class contains the code for the Particles in the swarm
# For each particle INSTANCE;
#    .pos[]      = particles' current x values vector
#    .velocity[] = particles' current x velocities vector
#    .pBest[]    = particles' vector of x's for the particles' best solution (from .pos[]])
#    .pBestE     = value (e.g., energy, criteria, ...) at pBest[]
#
# Uses the paramter list paramData{'variableName'}{'min'|'max'} = number

class Particle:
    def __init__(self):
        self.velocity = []  # these are INSTANCE variables!
        self.pos      = []
        self.pBest    = []
        self.pBestE   = 1e100   # saves the best energy for this particle, assigned later
        for i in range(dimension):
            self.pos.append(random.uniform(paramData[i]['min'],paramData[i]['max']))
            self.velocity.append(maxVelocity*random.uniform(paramData[i]['min'],maxVelocity*paramData[i]['max']))
            self.pBest.append(self.pos[i])
        return 
        
    def updatePositions(self):
        for i in range(dimension):
            self.pos[i] = self.pos[i] + self.velocity[i]   
        return
 
    def updateVelocities(self, gBest):
        for i in range(dimension):
            social    = c1*random.random()*(     gBest[i] - self.pos[i])
            cognitive = c2*random.random()*(self.pBest[i] - self.pos[i])
            self.velocity[i] = (w*self.velocity[i]) + social + cognitive
        return
 
    def satisfyConstraints(self):
        # This is where constraints are satisfied, but only for the pos values
        # If the pos value is out of range then put it at the edge of the range 
        # and make sure the velocity is in the proper direction and divide by 1/2.
        for i in range(dimension):
            if self.pos[i] >= paramData[i]['max']:
                self.pos[i] = paramData[i]['max'] * 0.95
                self.velocity[i] = -1.*abs(self.velocity[i])/2.
            elif self.pos[i] <= paramData[i]['min']:
                self.pos[i] = paramData[i]['min'] * 1.05
                self.velocity[i] = abs(self.velocity[i])/2.
        return
 
# This class contains the particle swarm optimization algorithm
class ParticleSwarmOptimizer:
    solution = []
    swarm    = [[] for x in range(numSwarms)]
    gBest    = [[1e99] for x in range(numSwarms)]
    gBestE   = [ 1e99 for x in range(numSwarms)]
    lastIterWchange = 0  
    ggBestE = 1e99    

    theResults=[[1e99]*swarmSize for _ in range(numSwarms)]
        
    def __init__(self):   # create the swarms and particles
        for ii in range(numSwarms):
            self.swarm[ii] = []
            for h in range(swarmSize):
                particle=Particle()
                self.swarm[ii].append(particle)
            #print "swarm and size=",ii,len(self.swarm)
        return

    ####################################################################################
    # Update the parameters in the *_parameters.py file for the CC3D jobs.
    # The new values are just added to the end of the existing file in each of the job
    # directories.  Since the parameters file is a python script the last assignments is what
    # is used during the simulation. In addition, this method retains a history of all
    # the parameter sets (and energy value, see below) for each particle.
    # Also add the best energy for the particular particle so it can be used to terminate
    # CC3D runs when the residual error exceeds the previous best residual error
    # Use the dirArray to map swarm + particle numbers to the proper directory name
    def updateParams(self,i,logResults,doAll):
        for iS in range(numSwarms):
            for iP in range(swarmSize):
                pFile = dirArray[iS][iP] + '/' + cc3dParameterFile
                #print '          Parameter write:',iS,iP,pFile
                paramFile = open(pFile,'a')
                paramFile.write("\nparticleLastE = "+str(self.theResults[iS][iP])+"   # E for previous data set\n")
                paramFile.write("# History = "+logResults[iS][iP]+"\n")   # NEW
                if doAll:
                    paramFile.write("\n######################## swarm iteration = "  \
                                    +str(i) +" DupRun = UNK  ###############################\n")
                    paramFile.write("particleBestE = "+str(self.swarm[iS][iP].pBestE)+"\n")
                    for ind in range(dimension):
                        #print key,"=",self.swarm[iS][iP].pos[ind]
                        paramFile.write(paramData[ind]['varName']+" = "+str(self.swarm[iS][iP].pos[ind])+"\n")
                    paramFile.write("swarm = "+str(iS)+"\n")
                    paramFile.write("particle = "+str(iP)+"\n")
                    paramFile.write("iter = "+str(i)+"\n")
                paramFile.close()
        return

    ############################################################################################### 
    # running log of each swarm's best position
    #
    def runningLog(self,iter,elapsedTime):
        path='/u/jsluka/swarm/swarm_running_result.txt'
        oFile = open(path,'a',1)  # line buffered
        if iter == 0:
            theLine =" iter swrm prtcl  "
            for ind in range(dimension):
                theLine +=paramData[ind]['varName']+"  "
            theLine +="Qual  Time  pBestE  sBestE  ggBestE  "
            theLine +="Individual_Qual_for_multiple_runs/parameter_set"
            oFile.write(theLine+'\n')

        # find the global-global best E across all swarms
        for iS in range(numSwarms):
            for iP in range(swarmSize):
                if self.theResults[iS][iP] <= self.ggBestE:
                    self.ggBestE = self.theResults[iS][iP]

        for iS in range(numSwarms):
            for iP in range(swarmSize):
                theLine  = '%5i ' % iter
                theLine += '%4i ' % iS
                theLine += '%4i   ' % iP
                theLine += ' '.join('%.3e ' % v for v in self.swarm[iS][iP].pos)
                theLine += " %.4e "  % self.theResults[iS][iP]
                theLine += " %8.1f " % self.theTimings[iS][iP]   # NEW
                if self.theResults[iS][iP] == self.swarm[iS][iP].pBestE:
                    theLine += "   1   "
                else:
                    theLine += "   --  "
                if self.theResults[iS][iP] == self.gBestE[iS]:
                    theLine += "   1   "
                else:
                    theLine += "   --  "
                if self.theResults[iS][iP] == self.ggBestE:
                    theLine += "    1   "
                else:
                    theLine += "    --  "
                theLine += "     "+self.logResults[iS][iP] # self.logResults[iS][iP]
                oFile.write(theLine+'\n')
            oFile.write('\n')
        oFile.close()
        return
        
    ###############################################################################################     
    def optimize(self):      
        # logResults is a 2D array of strings. Each string is the list of previous quality values
        self.logResults=[[""]*swarmSize for _ in range(numSwarms)]  # NEW this gets reset below
        for iter in range(iterations):   # swarm iterations
            self.updateParams(iter,self.logResults,True)
            swarmStartTime = timeit.default_timer()
            # do numDuplRuns replicate runs of each paramter set / swarm.
            self.cumulativeResults=[[0.]*swarmSize for _ in range(numSwarms)]
            self.cumulativeTimings=[[0.]*swarmSize for _ in range(numSwarms)]  # NEW
            self.logResults=[[""]*swarmSize for _ in range(numSwarms)]  # NEW 
            for dupRun in range(0,numDuplRuns):
                print "                 +++ duplicate run iteration:",dupRun+1
                self.theResults,self.theTimings = submitJobs(iter,dupRun)  # do most of the work 
                for iS in range(numSwarms):   # NEW
                    for iP in range(swarmSize):
                        self.cumulativeResults[iS][iP] += self.theResults[iS][iP]
                        self.cumulativeTimings[iS][iP] += self.theTimings[iS][iP]
                        self.logResults[iS][iP]=self.logResults[iS][iP]+str(self.theResults[iS][iP])+"  "

            # calculate the average Quality value for each particle across the replicate runs  
            for iS in range(numSwarms):
                for iP in range(swarmSize):
                    self.theResults[iS][iP] = self.cumulativeResults[iS][iP]/(dupRun+1.)
                    self.theTimings[iS][iP] = self.cumulativeTimings[iS][iP]/(dupRun+1.)
            
            swarmElapsedTime = timeit.default_timer() - swarmStartTime  #NOTE: The MAXIMUM run time for all particles
            print "                                          swarmElapsedTime=",swarmElapsedTime
            
            for iS in range(numSwarms):  
                if iter%5 == 0:
                    print "iter %4i swrm %2i gBestE %9.3f  " % (iter,iS,self.gBestE[iS]),   
                    print ' '.join('%8.3f' % v for v in self.gBest[iS])    

                #Update the personal (particle) best positions
                for iP in range(swarmSize):
                    ###E = self.f(self.swarm[iS][iP].pos)
                    E = self.theResults[iS][iP]
                    #print "------------------iS,iP,E,pBesE",iS,iP,E,self.swarm[iS][iP].pBestE
                    if E < self.swarm[iS][iP].pBestE:
                        self.swarm[iS][iP].pBest = list(self.swarm[iS][iP].pos)   
                        self.swarm[iS][iP].pBestE = E
                        print "   new particle gBest",[iS],[iP],'E=%8.2f (' % E,
                        print ' '.join(' %6.3f' % v for v in self.swarm[iS][iP].pBest),
                        print ')'

                #Get the global (within the particular sub-swarm) best particle
                for iP in range(swarmSize):
                    ###E = self.f(self.swarm[iS][iP].pBest) 
                    E = self.theResults[iS][iP]                    
                    if E < self.gBestE[iS]:
                        self.gBest[iS] = list(self.swarm[iS][iP].pBest)  
                        self.gBestE[iS] = E    
                        print "new SWARM gBestE,gBest",[iS],[iP],'E=%8.2f (' % self.gBestE[iS],
                        print ' '.join('%7.3f' % v for v in self.gBest[iS]),
                        print ')'
                        lastIterWchange = iter
            
            self.runningLog(iter,swarmElapsedTime)
            
            for iS in range(numSwarms):              
                #Update position of each particle
                for iP in range(swarmSize):
                    self.swarm[iS][iP].updateVelocities(self.gBest[iS])
                    self.swarm[iS][iP].updatePositions()
                    self.swarm[iS][iP].satisfyConstraints()

        self.updateParams(iter,self.logResults,False)
        return self.gBest, lastIterWchange, self.gBestE, iter     
 
##    def f(self, solution):  #This is where the quality of fit is defined
##        # Quadratic in n dimensions with many nodes 
##        qual = (solution[0]-1.000)**2 + (solution[1]-0.002)**2  \
##             + (solution[2]-0.020)**2 + (solution[3]-0.100)**2  \
##             + (solution[4]-0.002)**2 + (solution[5]-6.000)**2  \
##             + (solution[6]-0.010)**2 + (solution[7]-1.000)**2  \
##             + (solution[8]-0.001)**2 + (solution[9]-0.005)**2  \
##             + (solution[10]-100.)**2 #+ (solution[11]-100.)**2  # last param has no effect
##
##        return qual
 
def main():
    global paramData, dimension
    paramData=[]
    paramData=loadParameters()
    dimension = len(paramData)
    print "Parameter count:",dimension,"\nParameter Data:\n",paramData,"\n"

    global dirArray
    dirArray =[[""]*swarmSize for _ in range(numSwarms)]
    dirArray=makeWorkingDirs()

    pso = ParticleSwarmOptimizer()
    (final,lastIterWchange,bestE,lastIter) = pso.optimize()     
    print "final:",final
    print "lastIterWchange:",lastIterWchange
    print "bestE:",bestE
    print "lastIter:",lastIter
    
    # id the best particle across all the swarms
    bestVals = []
    finalBestE=bestE[0]
    for iS in range(numSwarms):              
        if finalBestE >= bestE[iS]:
            finalBestE = bestE[iS]
            iSbest = iS
    bestVals = list(final[iSbest])
    print "iSbest=",iSbest
    print "bestVals=",bestVals

    ###########################################################################
    # Final output
    # Include the results, python'ish listing of the parameters and ranges etc.
    
    # Stats across the best solution for each swarm
    stat=[]
    for ii in range(dimension):
        stat.append([])
        for jj in range(numSwarms):
            stat[ii].append(final[jj][ii])
    import numpy
    av=[0 for x in range(dimension)]
    sd=[0 for x in range(dimension)]
    psd=[0 for x in range(dimension)]
    for ii in range(dimension):
        av[ii]=numpy.average(stat[ii])
        sd[ii]=numpy.std(stat[ii])
        psd[ii]=sd[ii]/av[ii]*100.
        
    print "\nLast iteration that improved:",lastIterWchange+1,"\n"

    path='/u/jsluka/swarm/swarm_final_result.txt'
    oFile = open(path,'w')
    oFile.write('Swarm Summary:\n')
    oFile.write("numSwarms="+str(numSwarms)+"  swarmSize="+str(swarmSize)+"  Max iterations="+str(iterations))
    oFile.write("\nActual iterations="+str(lastIter+1)+"  lastIterWchange="+str(lastIterWchange+1)+'\n')
    oFile.write('w='+str(w)+'  c1='+str(c1)+'  c2='+str(c2)+'  maxVelocity='+str(maxVelocity)+'\n')
    #swarmStopTime='%s' % datetime.datetime.now()
    swarmStopTime=datetime.datetime.now()
    oFile.write('swarmStartTime='+'%s' % swarmStartTime+'\n')
    oFile.write(' swarmStopTime='+'%s' % swarmStopTime+'\n')  
    oFile.write('  Elapsed Time='+str(swarmStopTime-swarmStartTime)+'\n')
#    oFile.write('     deltaTime='+swarmStopTime-swarmStartTime+'\n')  

    oFile.write('Best solution swarm='+str(iSbest)+'\n\n')

    oFile.write('# Python Parameter list with ranges and original comments.  Quality Value='+str(finalBestE)+'\n')
    hasFlag = False
    for ind in range(dimension):
        oFile.write('%-22s = ' % (paramData[ind]['varName']))
        oFile.write('%.4e'     % (bestVals[ind]))
        oFile.write(" # min:"+str(paramData[ind]['min']))
        oFile.write(" max:"+str(paramData[ind]['max']))
        # check if the values is near the min or max
        flag = ""
        if (1.05 >= bestVals[ind]/paramData[ind]['min'])   \
        or (0.95 <= bestVals[ind]/paramData[ind]['max']):
            flag = " ***** "
            hasFlag = True
        oFile.write(flag)
        oFile.write("## "+paramData[ind]['comments']+'\n')
    if hasFlag:
        oFile.write("################# NOTE: "+flag+" marks parameters that were near their min or max values #############################")
    oFile.write("\n\n")
    
    oFile.write('Swarm Details:\n')
    oFile.write("Swarm ")
    for ind in range(dimension):
        oFile.write(paramData[ind]['varName']+"  ")
    oFile.write("\n")
    
    # output to the terminal and the file
    for ii in range(numSwarms):
        theLine  = '%3i ' % ii
        theLine += ' '.join('%8.4e' % v for v in final[ii])   
        theLine += " bestE: %8.4e" % bestE[ii]
        if bestE[ii] == finalBestE:
            theLine += " GMin"
        print theLine
        oFile.write(theLine+'\n')
        
    # output to the terminal and log file
    # print the stats
    print "==================================================================================="
    print "avg",' '.join(' %8.4e '    % v for v in av)
    print "sd ",' '.join(' %8.4e '    % v for v in sd)
    print "%sd",' '.join('  %7.1f%% ' % v for v in psd)

    oFile.write("avg"+' '.join(' %8.4e '    % v for v in av)+'\n')
    oFile.write("sd "+' '.join(' %8.4e '    % v for v in sd)+'\n')
    oFile.write("%sd"+' '.join('  %7.1f%% ' % v for v in psd)+'\n')
    oFile.write('\n')

    oFile.close()

if  __name__ =='__main__':
    main()