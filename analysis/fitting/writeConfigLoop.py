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
    
    className = "vec_ps_refl" # AMPTOOLS_AMPS class definition
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
    Mloop = ""
    Lloop = ""
    comments = ""
    fixedPhase = ""
    fixedPhaseRefl = []

    # loop over waves
    jpComment = ""
    for wave in waves:
        j = wave["spin"]
        parity = wave["parity"]
        jp = str(wave["spin"]) + char[wave["parity"]]
        l = L[wave["l"]]
        
        jpComment = "############################ spin %d parity %+d ##################################\n\n" % (j,parity)
        if jpComment not in comments:
            fout.write(jpComment)
            comments += jpComment
        
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

            # write individual amplitudes
            sum = word[refl] + "Refl" + word[sign] + "Sign"
            amplitude = reaction + "::" + sum + "::LOOPAMPNAME"
            amplitudeLine = "amplitude " + amplitude + " " + className
            amplitudeLine += " %d LOOPM LOOPL  %+d  %+d  %s" % (j, real, sign, common)
            if amplitudeLine not in amplitudeLines:
                amplitudeLines += amplitudeLine
                amplitudeLines += "\n"
            
            # write amplitude initialization
            initializationLine = "initialize " + reaction + "::" + sum + "::LOOPAMPNAME"
            if initRefl and refl != initRefl:
                initializationLine += " cartesian 0 0 fixed"
            else:
                initializationLine += " cartesian %d %d" % (initReMag, initImMag)
            
            # separate sums
            if initializationLine not in initializationLines:
                initializationLines += initializationLine
                initializationLines += "\n"
            
            # loop over spin projections
            for spin_proj in range(j,-(j+1),-1):
                amp = jp+char[spin_proj]+l
                if amp not in amplitudes:
                    amplitudes.append(amp)
                    if spin_proj == 0:
                        Mloop += " 0"
                    else:
                        Mloop += " %+d" % spin_proj
                    Lloop += " %d" % wave["l"]
                    
                # select amplitude for fixed phase
                if refl not in fixedPhaseRefl and spin_proj == 0: # fix phase for one of the amplitudes in each reflectivity
                    fixedPhaseRefl.append(refl)
                    fixedPhase += initializationLine.replace("LOOPAMPNAME",amp)
                    fixedPhase += " real\n"
            
    # write loop info to file
    loopAmpName = "loop LOOPAMPNAME"
    for amp in amplitudes:
        loopAmpName += " " + amp
    fout.write(loopAmpName + "\n")
    fout.write("loop LOOPM%s\n" % Mloop)
    fout.write("loop LOOPL%s\n" % Lloop)
        
    # write all amplitudes and initialization to file
    fout.write(amplitudeLines)
    fout.write(initializationLines)
    
    # constrain common parameters from different sums
    fout.write("constrain omegapi PosReflPosSign LOOPAMPNAME omegapi PosReflNegSign LOOPAMPNAME\n")
    fout.write("constrain omegapi NegReflPosSign LOOPAMPNAME omegapi NegReflNegSign LOOPAMPNAME\n")
    
    # option to constrain phase between S and D wave amplitudes (use real-valued D/S ratio as scale)
    if not constrainDSamps:
        return
    
    # loop over amplitudes and constrain between S and D waves from same sum
    fout.write("\n# constrain S and D waves to the same amplitude and set scale factor for D/S ratio\n")
    fout.write("loop LOOPSUM NegReflNegSign PosReflNegSign PosReflPosSign NegReflPosSign\n")
    previousSumName = ""
    scaleLines = ""
    scaleParName = "dsratio"
    for amp1 in amplitudes:
        for amp2 in amplitudes:
                
            # same amplitude except for L
            if amp1[:-1] == amp2[:-1] and (amp1[-1]=="s" and amp2[-1]=="d"):
                fout.write("constrain %s LOOPSUM %s %s LOOPSUM %s\n" % (reaction, amp2, reaction, amp1))
                scaleLines += "scale " + reaction + " LOOPSUM " + amp2 + " [" + scaleParName + "]\n"
                    
    fout.write(scaleLines)
    
    fout.write("\n# fix phase\n")
    fout.write(fixedPhase)

if __name__ == "__main__":
    main(sys.argv[1:])
