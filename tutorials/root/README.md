Documentation for learning ROOT with online tutorials:

Most of the GlueX analysis software (and HEP/NP software in general) is based around the CERN ROOT package found here https://root.cern

# Interactive student course https://github.com/root-project/student-course

#* 

# Install ROOT on you local machine, following the instructions at https://root.cern/install/

#* For MacOS users I recommend the Homebrew package manager https://brew.sh  Follow the instructions to install homebrew then you can install other packages from the command line with formula like https://formulae.brew.sh/formula/root

`brew install root`

# ROOT Primer https://root.cern/primer/



# Finally, you should copy the file .rootrc into your home directory with the command

`cp $WM_GLUEX/.rootrc ~/`

which will load a custom ROOT environment everytime you open a session.  This is needed to properly access some of the libraries in the tutorials below.
