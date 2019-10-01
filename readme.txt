Everything (as of 10/1/2019) should be organized by examples in the IEEE paper version. The first example comes from the folder L1denoise, where we generate a random sparse signal and add large sparse noise to it. This can be done by running alg_comp.m (for figure 1) and signal_denoise.m for figures 2-4. Figure 5 is created by running signal_denoise_phi_explore.m. You can try different combinations of phi/psi (objective/constraint respectively) by changing the cell with those respective names. 

The second example is located in "curveletexamplecodeanddata". This should reproduce the curvelet example in the BPDN paper. Note that one will have to download CVX and SPGL1 for comparison runs. Otherwise, change phi and psi as above.

The third example is low-rank interpolation on a 4d model (aptly named 4d). The code that generates the figures and same them in the appropriate figs file is "4d_figure_gen.m". All examples call respective algorithms located in 'alg_tools.m', or in the case of the 4D model, 'fcns'. Note that the true data file is very large, and not stored in GitHub. Please send an email to rbaraldi@uw.edu to receive it. 


However, note that to run the SPGL1, one will have to download and install spgl1 and SPOT operator libraries: https://www.cs.ubc.ca/~mpf/spgl1/ and https://github.com/mpf/spot respectively. 


For comments/questions/errors (hopefully there are no errors - that would be pretty bad), contact rbaraldi@uw.edu. Please allow a couple days for a response. 