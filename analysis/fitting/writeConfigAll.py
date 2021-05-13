#!/usr/bin/env python
import sys
import os
import subprocess
import math
from optparse import OptionParser

####################### MAIN #######################
def main(argv):
    
    name = "neutralb1_all"
    amps = "1p"
    parser_usage = "writeConfigLoop.py [template_name amps]"
    parser = OptionParser(usage = parser_usage)
    (options, args) = parser.parse_args(argv)
    if len(args) > 0:
	name = args[0]
	if len(args) > 1:
	    amps = args[1]

    reaction = "omegapi"
    configName = "template"  # file label for output configuration
    templateName = "template_%s.cfg" % name  # name of template config file defining common parameters

    className = "Vec_ps_refl" # AMPTOOLS_AMPS class definition
    constrainDSamps = True    # option to fix phase between S and D waves
    forceRefl = 0             # consider only one reflectivity by setting to +1 or -1 (default 0 includes both)
    initRefl = 0              # initialize only one reflectivity by setting to +1 or -1 (default 0 includes both)
    initReMag = 100
    initImMag = 100
    
    outFileName = "fit_" + reaction + "_amplitude_"
    if initRefl: outFileName += "refl%s_" % ("%+d" % initRefl)[0]
    outFileName += configName
    
    # add desired waves to list (could come from command line or short file input)
    waves = []
    if "0m" in amps: # (special case without loops)
	waves.append( {"spin":0, "parity":-1, "l":1} ) #0- P-wave
    if "1p" in amps:
    	waves.append( {"spin":1, "parity":+1, "l":0} ) #1+ S-wave
    	waves.append( {"spin":1, "parity":+1, "l":2} ) #1+ D-wave
    if "1m" in amps:
    	waves.append( {"spin":1, "parity":-1, "l":1} ) #1- P-wave
    if "2p" in amps:
        waves.append( {"spin":2, "parity":+1, "l":2} ) #2+ D-wave
    if "2m" in amps:
        waves.append( {"spin":2, "parity":-1, "l":1} ) #2- P-wave
	waves.append( {"spin":2, "parity":-1, "l":3} ) #2- F-wave
    if "3m" in amps:
	waves.append( {"spin":3, "parity":-1, "l":3} ) #3- F-wave
    if "iso" in amps:
        waves.append( {"spin":-9, "parity":0, "l":-9} ) # isotropic background
    
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
    wordcomp = {-1:"Imag", +1:"Real"}
    L = {0:"s", 1:"p", 2:"d", 3:"f", 4:"g", -9:"b"}

    common = " LOOPPOLANG LOOPPOLVAL dalitz"
    comments = ""
    fixedPhase = ""
    fixedPhaseRefl = []

    # list of amplitudes for later constraints and scale factors
    amplitudes = []
    amplitudes_jp = []

    # loop over waves
    jpComment = ""
    for iwave in range(len(waves)):
        wave = waves[iwave]
        j = wave["spin"]
        parity = wave["parity"]
        jp = str(wave["spin"]) + char[wave["parity"]]
        l = L[wave["l"]]
           
        naturality = parity * math.pow(-1,j)

        # write background amplitude
        if wave["spin"] == -9:
            fout.write("############################ isotropic background ##################################\n\n")
            fout.write("sum LOOPREAC::Bkgd\n")
            fout.write("amplitude LOOPREAC::Bkgd::isotropic Uniform\n")
            fout.write("initialize LOOPREAC::Bkgd::isotropic cartesian 100 0 real\n")
            fout.write("constrain %s::Bkgd::isotropic LOOPREAC::Bkgd::isotropic\n" % reaction)
            fout.write("scale LOOPREAC::Bkgd::isotropic LOOPSCALE")
            
            continue

        jpComment = "############################ spin %d parity %+d ##################################\n\n" % (j,parity)
        if jpComment not in comments:
            fout.write(jpComment)
            comments += jpComment
        
        #if wave["spin"] == 0:
        #    write0m(fout)
        #    continue
        
        # temporary strings for storing amplitude lines and initialization to be written later
        amplitudeLines = ""
        initializationLines = ""
        constraintLines = ""
        constraintLineReactions = ""

        # 4 coherent sums for reflectivity and sign of (1 +/- P_gamma)
        for realsign in [[-1,-1,"ImagNegSign"],[+1,-1,"RealNegSign"],[-1,+1,"ImagPosSign"],[+1,+1,"RealPosSign"]]:
            real = realsign[0]
            sign = realsign[1]
            
            # real or imaginary for particular sum is given by naturality
            refl = 0
            if sign*real > 0:
                if naturality > 0: # natural
                    refl = +1
                else: # unnatural
                    refl = -1
            else:
                if naturality > 0: # natural
                    refl = -1
                else: # unnatural
                    refl = +1

            # option to force reflectivity if desired
            if forceRefl and refl != forceRefl:
                continue
                
            # write individual amplitudes
            sum = realsign[2] #wordcomp[real] + word[sign] + "Sign"
            for spin_proj in range(j,-(j+1),-1):
                if abs(spin_proj) > 2:
                    continue # conly consider m <=2 for photon beam
                    
                amp = jp+char[spin_proj]+l
                amplitudes.append(amp)
                amplitudes_jp.append(amp)

                amplitude = "LOOPREAC::" + sum + "::" + amp
                amplitudeLine = "amplitude " + amplitude + " " + className
                
                amplitudeLine += " %d %s %s  %+d  %+d  %s" % (j, spin_proj, wave["l"], real, sign, common)
                if amplitudeLine not in amplitudeLines:
                    amplitudeLines += amplitudeLine
                    amplitudeLines += "\n"
                        
                # write amplitude initialization
                initializationLine = "initialize LOOPREAC::" + sum + "::" + amp
                if initRefl and refl != initRefl:
                    initializationLine += " cartesian 0 0 fixed"
                else:
                    initializationLine += " cartesian %d %d" % (initReMag, initImMag)

                # write constraints
                constraintLine1 = "constrain LOOPREAC ImagNegSign %s LOOPREAC RealPosSign %s\n" % (amp, amp)
                constraintLine2 = "constrain LOOPREAC RealNegSign %s LOOPREAC ImagPosSign %s\n" % (amp, amp)
                if constraintLine1 not in constraintLines:
                    constraintLines += constraintLine1
                if constraintLine2 not in constraintLines:
                    constraintLines += constraintLine2

                # only need two sums since PosSign constrained above
                for loopSum in ["ImagNegSign", "RealNegSign"]:
                    constraintLine3 = "constrain omegapi %s %s LOOPREAC %s %s\n" % (loopSum, amp, loopSum, amp)
                    if constraintLine3 not in constraintLineReactions:
                        constraintLineReactions += constraintLine3
                
                    
                # separate sums
                if initializationLine not in initializationLines:
                    initializationLines += initializationLine
                    initializationLines += "\n"
                    
                # select amplitude for fixed phase
                if real not in fixedPhaseRefl and spin_proj == 0: # fix phase for one of the amplitudes in each reflectivity
                    fixedPhaseRefl.append(real)
                    initphase = initializationLine
                    fixedPhase += initphase.rsplit(' ', 1)[0] 
                    fixedPhase += " 0 real\n"
        
            # if more waves with the same J^P are in waveset, wait to print loop
            skipPrinting = False
            for jwave in range(iwave,len(waves)):
                if wave == waves[jwave]: continue
           
                later_jp = str(waves[jwave]["spin"]) + char[waves[jwave]["parity"]]
                if later_jp == jp:
                    skipPrinting = True
                
                if skipPrinting:
                    continue;
        
        
        # write all amplitudes and initialization to file
        fout.write(amplitudeLines)
        fout.write(initializationLines)
    
        # constrain common parameters from different sums
        fout.write(constraintLines)
        fout.write(constraintLineReactions)
            
        fout.write("\n")
        
        amplitudes_jp = []
        amplitudeLines = ""
        initializationLines = ""
        constraintLines = ""
        constraintLineReactions = ""

    fout.write("\n# fix phase\n");
    fout.write(fixedPhase);
    
    # option to constrain phase between S and D wave amplitudes (use real-valued D/S ratio as scale)
    scaleLines = ""
    if constrainDSamps:
    
        # loop over amplitudes and constrain between S and D waves from same sum
        fout.write("\n# constrain S and D waves to the same amplitude and set scale factor for D/S ratio\n")
        
        fout.write("parameter dsratio 0.27 bounded 0 1\n\n")
        fout.write("keyword parRange 3 3\n")
        fout.write("parRange dsratio 0.2 0.34\n\n")
        
        scaleParName = "dsratio"
        constrainLinesDS = ""
        for loopSum in ["ImagNegSign", "RealNegSign", "ImagPosSign", "RealPosSign"]:
            for amp1 in amplitudes:
                for amp2 in amplitudes:
                    
                    # same amplitude except for L
                    if amp1[:-1] == amp2[:-1] and (amp1[-1]=="s" and amp2[-1]=="d"):
                        constrainLineDS = "constrain LOOPREAC %s %s LOOPREAC %s %s\n" % (loopSum, amp2, loopSum, amp1)
                        if constrainLineDS in constrainLinesDS:
                            continue
                        
                        constrainLinesDS += constrainLineDS
                        scaleLines += "scale LOOPREAC::" + loopSum + "::" + amp2 + " [" + scaleParName + "]\n"

        fout.write(constrainLinesDS)

    # scale factors for each orientation
    for loopSum in ["ImagNegSign", "RealNegSign", "ImagPosSign", "RealPosSign"]:
        for amp in amplitudes:
            fullAmplitude = "LOOPREAC::" + loopSum + "::" + amp
            if fullAmplitude in scaleLines:
                continue
            scaleLines += "scale " + fullAmplitude + " LOOPSCALE\n"
    
    fout.write(scaleLines)

if __name__ == "__main__":
    main(sys.argv[1:])
