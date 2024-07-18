# VOI_theory

This repository contains code to generate the figures in our article "Why shouldn't I collect more data? Reconciling disagreements between intuition and value of Information analyses" in the Journal Methods in Ecology and Evolution, 2024, by Holden et al. 

You must save the file "VOI_functions.R" in the same directory as the other files.

The file which start with VOI_hist call VOI_functions to create the histograms of VOI in the manuscript for a 2 action 2 state VOI problem. The first row of Figure 2 in the main text is produced by VOI_hist_punif.R and the second row is produced by VOI_hist_pfixed.R. All other VOI_hist files produce the figures in the supplement.

The files EVPI_gam_states.R and EVPI_gam_actions.R creates simulated VOI values as you increase the number of states and actions respectively and saves the data in your file. Gam stands for gamma distributed utilities for which the exponential distribution is a special case of the gamma distribution (i.e. X~gamma(1,lambda) is equivalent to X~exponential(lambda)), and in our manuscript we use the exponential distribution, but the file has the functionality to explore other gamma distributed utilities, which we did not do.

The files make_asymp_plots_2x2_uniform.R and make_gamma_plots_2x2.R takes in the data created by running EVPI_unif_states.R, EVPI_unif_actions.R, EVPI_gam_states.R and EVPI_gam_actions.R to generate the 2x2 multi panel plots in the main text of the manuscript (Fig 4 and Fig 5).

The files exp_lambda_sim.R, unif_b_sim.R, and norm_std_sim.R, generate the suplimentary sensitivity figures to the paramerters of these distribution, 
which shows that the distribution of VOI values does not depend on these parameters

