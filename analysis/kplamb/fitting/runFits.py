#!/usr/bin/env python3

# Alex Austregesilo
# aaustreg@jlab.org

import argparse
import sys
import os
import subprocess
import math
from shutil import copyfile
#import ROOT

#########################################################################################################################

def getEntries(ROOTFILE):
        f = ROOT.TFile(ROOTFILE)
        t = f.Get("kin")
        return t.GetEntries()

########################################################## MAIN ##########################################################

def main(argv):
        global VERBOSE # so can modify here

        workingDir = os.getcwd()
        scriptDir = os.path.dirname(os.path.abspath(__file__))

        # default input: the analysis/ directory of this repository, which holds
        # kplamb/out/t_<low>_<high>/ after running kplamb/skimForAmpTools.C
        defaultInDir = os.path.abspath(os.path.join(scriptDir, "..", ".."))

        parser = argparse.ArgumentParser(
                description="Run the standalone KHyperon fit in each t bin.")
        parser.add_argument("--indir", default=defaultInDir,
                help="directory containing <channel>/out/t_<low>_<high>/ skims "
                     "[default: %(default)s]")
        parser.add_argument("--nrand", default="25",
                help="number of randomized fit trials per t bin [default: %(default)s]")
        parser.add_argument("--maxcalls", default="1000000",
                help="maximum MINUIT iterations [default: %(default)s]")
        args = parser.parse_args(argv)

        # channel specific information: names, trees, data location
        channel = "kplamb"
        tree = "kplamb__B4_M18"
        FSGlueX = "100000000_100000"
        FSGlueX_Gen = "100000000_100000"
        reacName = "Beam K+ Proton Pi-"

        # directory that data is located in
        inDir = os.path.abspath(args.indir)

        # define type of fit to be performed
        numRand = args.nrand
        maxCalls = args.maxcalls
        fitName = "khyperon"
        cfgTempl = os.path.join(scriptDir, "fit_" + fitName + ".cfg")
        fitBinary = os.path.join(scriptDir, "standalone", "fit")
        plotterBinary = os.path.join(scriptDir, "standalone", "khyperon_plotter")
        plotMacro = os.path.join(scriptDir, "plot_plotter.C")
        parameterSummaryScript = os.path.join(scriptDir, "plot_khyperon_params_vs_t.py")

        # define binning in t
        lowT = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.3]
        highT = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 2.1]

        # get number of bins for loops below
        nBinsT = len(lowT)

        # this is where the goodies for the fit will end up
        fitDir = workingDir + "/" + channel + "/" + fitName
        os.makedirs(fitDir, exist_ok=True)
        os.chdir(fitDir)

	# loop over t bins
        for binT in range(0,nBinsT):

                binT_name = "t_%0.2f_%0.2f" % (lowT[binT], highT[binT])
                os.makedirs(binT_name, exist_ok=True)
                os.chdir(binT_name)

                # check if data exists
                dataPath = inDir + "/" + channel + "/out/" + binT_name
                if not os.path.exists(dataPath):
                        print("Data path doesn't exist:", dataPath)
                        os.chdir("../")
                        continue

                # replace required elements of fit configuration file
                f = open(cfgTempl,'r')
                filedata = f.read()
                f.close()

                filedata = filedata.replace("FITNAME", fitName)
                filedata = filedata.replace("CHANNEL", channel)
                filedata = filedata.replace("TREE", tree)
                filedata = filedata.replace("FSGLUEX_GEN", FSGlueX_Gen)
                filedata = filedata.replace("FSGLUEX", FSGlueX)
                filedata = filedata.replace("REACNAME", reacName)

                filedata = filedata.replace("INDIR", inDir)

                filedata = filedata.replace("NIFILE", fitName+"_NI")
                filedata = filedata.replace("TBIN", binT_name)

                cfgBin = binT_name + ".cfg"
                f = open(cfgBin,'w')
                f.write(filedata)
                f.close()

                # execute fit
                subprocess.run([fitBinary,"-c",cfgBin,"-r",numRand,"-m",maxCalls])

                # execute plotter
                subprocess.run([plotterBinary,"%s.fit" % fitName])

                # make summary PDFs from the best-fit ROOT output
                subprocess.run([
                        "root", "-l", "-b", "-q",
                        '%s("%s.plots.root")' % (plotMacro, fitName)
                ])

                # return for next t bin
                os.chdir("../")

        # make final CSV and summary PDF/ROOT graphs of Sigma, Ox, P, T, Oz vs t
        subprocess.run([sys.executable, parameterSummaryScript, fitDir])

if __name__ == "__main__":
   main(sys.argv[1:])
