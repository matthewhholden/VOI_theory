# VOI_theory

This repository contains code to generate the figures in our article "Why shouldn't I collect more data? Reconciling disagreements between intuition and value of Information analyses" in the Journal Methods in Ecology and Evolution, 2024, by Holden et al. 

You must save the file "VOI_functions.R" in the same directory as the other files.

The file names that start with VOI_hist call VOI_functions to create the histograms of EVPI in the manuscript for a 2 action 2 state problem. The first row of Figure 2 in the main text is produced by VOI_hist_punif.R and the second row is produced by VOI_hist_pfixed.R. Figure 3, the histogram of pstar is produced at the end of the All other VOI_hist_pfixed.R, but because pstar does not depend on p, we could have wrote this code in either of the above two files. The other files that start with VOI_hist produce the histogram figures in the supplement with the corresponding distributions on p given in the file name and figure caption.

The files EVPI_gam_states.R and EVPI_gam_actions.R create simulated EVPI values as you increase the number of states and actions respectively and saves the data to a new file in the directory. Gam stands for gamma distributed utilities for which the exponential distribution is a special case of the gamma distribution (i.e. X~gamma(1,lambda) is equivalent to X~exponential(lambda) ), and in our manuscript we use the exponential distribution, but the file has the functionality to explore other gamma distributed utilities, which we did not do, by varying the first paramater in the call to rgamma.

The files "make_asymp_plots_2x2_uniform.R" and "make_gamma_plots_2x2.R" take in the data created by running EVPI_asym_states.R, EVPI_asym_actions.R, EVPI_gam_states.R and EVPI_gam_actions.R to generate the 2x2 multi panel plots in the main text of the manuscript (Fig 4 and Fig 5).

The files exp_lambda_sim.R, unif_b_sim.R, and norm_std_sim.R, generate the suplimentary sensitivity figures to the paramerters of these distribution, 
which shows that the distribution of VOI values does not depend on these parameters

R is particularly not very good at 3D plotting so in figure 1 of the manuscript we used MATLAB to plot the analytic EVPI solutions, the code of which is in Plot_VOIbound.m. However since it is just plotting an analytic function, this is trivial to replicate in any language or even a web aplet like Desmos.

