#!/usr/bin/env python
import sys
import os
import subprocess
import math
from optparse import OptionParser

####################### MAIN #######################
def main(argv):
    
    reaction = "omegapi"
    configName = "DeltaLowerVertex_b1"  # file label for output configuration
    templateName = "template_deltaPlusPlus_b1Minus.cfg"  # name of template config file defining common parameters
    
    className = "Vec_ps_refl" # AMPTOOLS_AMPS class definition
    constrainDSamps = True    # option to fix phase between S and D waves
    forceRefl = 0             # consider only one reflectivity by setting to +1 or -1 (default 0 includes both)
    initRefl = 0              # initialize only one reflectivity by setting to +1 or -1 (default 0 includes both)
    initReMag = 100
    initImMag = 0
    
    outFileName = "fit_" + reaction + "_amplitude_"
    if initRefl: outFileName += "refl%s_" % ("%+d" % initRefl)[0]
    outFileName += configName
    
    # add desired waves to list (could come from command line or short file input)
    waves = []
    waves.append( {"spin":1, "parity":+1, "l":0} ) #1+ S-wave
    waves.append( {"spin":1, "parity":+1, "l":2} ) #1+ D-wave
    # waves.append( {"spin":1, "parity":-1, "l":1} ) #1- P-wave
    
    # read template data for initial portion of config file defining common parameters and input file setup
    cfgTempl = templateName
    ftemplate = open(cfgTempl,'r')
    templatedata = ftemplate.read()
    ftemplate.close()

    # write template data to config file
    cfgBin = outFileName + ".cfg"
    fout = open(cfgBin,'w')
    fout.write(templatedata)
    
    # write amplitudes based on above input
    writeAmplitudes(waves, reaction, className, fout, forceRefl, initRefl, initReMag, initImMag, constrainDSamps)
    fout.close()
    
####################### WRITE AMPLITUDES #######################
def writeAmplitudes(waves, reaction, className, fout, forceRefl, initRefl, initReMag, initImMag, constrainDSamps):

    # dictionaries to write consistent names
    char = {-2:"m2", -1:"m", 0:"0", +1:"p", +2:"p2"}
    word = {-1:"Neg", +1:"Pos"}
    L = {0:"s", 1:"p", 2:"d", 3:"f", 4:"g"}
    common = "angle fraction dalitz"

    # list of amplitudes for later constraints and scale factors
    amplitudes = []

    # temporary strings for storing amplitude lines and initialization to be written later
    amplitudeLines = ""
    initializationLines = ""

    # loop over waves
    iamplitude = 0
    for wave in waves:
        j = wave["spin"]
        parity = wave["parity"]
        jp = str(wave["spin"]) + char[wave["parity"]]
        l = L[wave["l"]]
        
        jpComment = "############################ spin %d parity %+d ##################################\n" % (j,parity)
        if jpComment not in amplitudeLines:
            amplitudeLines += jpComment
    
        # 4 coherent sums for reflectivity and sign of (1 +/- P_gamma)
        for reflsign in [[-1,-1],[+1,-1],[+1,+1],[-1,+1]]:
            refl = reflsign[0]
            sign = reflsign[1]
            
            # option to force reflectivity if desired
            if forceRefl and refl != forceRefl:
                continue
            
            # real or imaginary for particular sum is given by naturality
            naturality = parity * math.pow(-1,j)
            real = 0
            if sign*refl > 0:
                if naturality > 0: # natural
                    real = +1
                else: # unnatural
                    real = -1
            else:
                if naturality > 0: # natural
                    real = -1
                else: # unnatural
                    real = +1

            # loop over spin projections
            for spin_proj in range(j,-(j+1),-1):
                
                spin_proj_string = "%+d" % spin_proj
                if spin_proj==0:
                    spin_proj_string = " 0"
                    
                # write individual amplitudes
                sum = word[refl] + "Refl" + word[sign] + "Sign"
                amp = char[refl]+jp+char[spin_proj]+l
                amplitude = reaction + "::" + sum + "::" + amp
                amplitudeLine = "amplitude " + amplitude + " " + className
                amplitudeLine += " %d %s %d %d %d %s" % (j, spin_proj_string, wave["l"], real, sign, common)
                amplitudeLines += amplitudeLine + "\n"
                amplitudes.append([sum,amp])
                iamplitude = iamplitude + 1
                   
                # write amplitude initialization
                initializationLine = "initialize " + reaction + "::" + sum + "::" + amp
                if initRefl and refl != initRefl:
                    initializationLine += " cartesian 0.0 0.0"
                else:
                    initializationLine += " cartesian %0.1f %0.1f" % (initReMag, initImMag)
                
                if iamplitude == 2: # fix phase for one of the amplitudes
                    initializationLine += " real"
                
                initializationLines += initializationLine + "\n"
            
            # separate sums
            amplitudeLines += "\n"
            initializationLines += "\n"
            
    # write all amplitudes and initialization to file
    fout.write(amplitudeLines)
    fout.write(initializationLines)
            
    # loop over amplitudes and constrain common parameters
    fout.write("\n# 12 constraints for same amplitudes in different terms")
    previousSumName = ""
    for i in range(len(amplitudes)):
        for j in range(i,len(amplitudes)):
            sum1 = amplitudes[i][0]
            sum2 = amplitudes[j][0]
            amp1 = amplitudes[i][1]
            amp2 = amplitudes[j][1]
        
            # don't need constraint if in the same sum
            if sum1 == sum2:
                continue
            
            # same amplitude in different sum needs constraint
            if amp1 == amp2:
                if previousSumName != sum1:
                    fout.write("\n")
                fout.write("constrain %s %s %s %s %s %s\n" % (reaction, sum2, amp2, reaction, sum1, amp1))
                previousSumName = sum1
    
    # option to constrain phase between S and D wave amplitudes (use real-valued D/S ratio as scale)
    if not constrainDSamps:
        return
    
    # loop over amplitudes and constrain between S and D waves from same sum
    fout.write("\n# 12 constraints to fix phase between S and D waves of 1+")
    previousSumName = ""
    scaleLines = ""
    scaleParName = "dsratio"
    for i in range(len(amplitudes)):
        for j in range(i,len(amplitudes)):
            sum1 = amplitudes[i][0]
            sum2 = amplitudes[j][0]
            amp1 = amplitudes[i][1]
            amp2 = amplitudes[j][1]
            
            # only interested in amplitudes from the same sum
            if sum1 != sum2:
                continue
                
            # same amplitude except for L
            if amp1[:-1] == amp2[:-1] and (amp1[-1]=="s" and amp2[-1]=="d"):
                if previousSumName != sum1:
                    fout.write("\n")
                fout.write("constrain %s %s %s %s %s %s\n" % (reaction, sum2, amp2, reaction, sum1, amp1))
                
                if previousSumName != sum1:
                    scaleLines += "\n"
                scaleLines += "scale " + reaction + " " + sum1 + " " + amp2 + " [" + scaleParName + "]\n"
                
                previousSumName = sum1
    
    # write final scale factors
    fout.write("\n# 12 scale factors for D/S ratio of 1+")
    fout.write(scaleLines)
    

if __name__ == "__main__":
    main(sys.argv[1:])
