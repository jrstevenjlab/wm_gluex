# Documentation for learning ROOT with online tutorials:

Most of the GlueX analysis software (and HEP/NP software in general) is based around the [CERN ROOT package](https://root.cern)

## [Interactive student course](https://github.com/root-project/student-course)

*  The course linked above, provides an interactive tutorial of the ROOT package without needing to install any software!
*  Open the software for the tutorial launching it in [Binder](https://mybinder.org/v2/gh/root-project/student-course/main)
*  Follow the [video recording](https://videos.cern.ch/record/2300516) of the tutorial and complete all th exercises

## Install ROOT on you local machine, [following the instructions](https://root.cern/install/) 

* To learn more about ROOT and more effectively use it you'll want a local installation on your machine
* For MacOS users I recommend the [Homebrew package manager](https://brew.sh)
* Follow the instructions to install homebrew then you can install other packages from the command line with [formulas](https://formulae.brew.sh/formula/root)

`brew install root`

## [ROOT Primer](https://root.cern/primer/)

* Now that you're familiar with ROOT from the tutorial and have a local installation, you can execute ROOT commands locally.
* Follow the ROOT Primer linked above to further explore how to use ROOT on you local machine
* Write your own macros to fill and display histograms and fit them, following the more extensive set of [tutorials](https://root.cern/tutorials/)

## Finally, you should copy the file .rootrc into your home directory with the command

`cp $WM_GLUEX/.rootrc ~/`

which will load a custom ROOT environment everytime you open a session.
