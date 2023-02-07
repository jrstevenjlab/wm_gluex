/*
  This file will read in a set of fit fractions (and other parameters)
  produced by the gen_vec_ps fitter, and create a set of plots showing
  how the parameters of interest are affected by seed and event number
  Original goal of this was to determine severity of m-projection
  leakage in vector-pseudoscalar input/output test.
 */

/*
  Past a certain number of seeds, I'll want to switch from lowStats
  plots (individual seeds shown) to highStats (every seed is histo-
  -grammed, depending on situation. see details below)

  Functionality Goals:
  	1. Every parameter of interest has a plot/histogram showing
	   its value and error as function of a few seeds / many 
	   seeds respectively. A similar histogram can be made for 
	   varying event # across many seeds, whose gaussian mean
	   and width can be plotted as function of event # for that
	   parameter.
	2. Non-generated parameters will be grouped into bar plots
	   or 2d projections showing their percent contribution/
	   clusters respectively
 */
/*
  Idea of how code will run
  	Begin by reading how many seeds there are, to determine 
	low or high stats type plots. From there extract relevant
	values from the fitPars files of each seed and plot them
	accordingly
 */
