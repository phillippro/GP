# GJ699 GP test

This work is to investigate the influence of noise modeling on signal detection in radial velocity (RV) data. To do this, I choose three RV data sets obtained by CARMENES, HARPS, and KECK/HIRES, for the Barnard's Star (GJ699). I test five different noise models by fit them in combination with different numbers of sinusoidal signals to different synthetic data sets. These noise models are white noise model (or MA(0)), first and second order moving average models (or MA(1) and MA(2)), Gaussian process with two different types of kernels (denoted by "GPabs" and "GPqp", which are the first and the third ones in eqn 15 of Feng, F., Tuomi, M., Jones, H. R. A., Butler, R. P., & Vogt, S. 2016, MNRAS,461, 2440. )I calculate their BIC-based Bayes factors (BFs) to investigate whether they lead to false negatives/positives. The specific steps are as follows. 

## 1. Synthetic data

I generate six synthetic data: W0P, MA0P, GP0P, W2P, MA2P, GP2P. In these names, "W" means white noise consistent with the jitter in real data, "MA" means MA noise determined by fiting MA(1)+2planets to the real data, "GP" means GP noise determined by fitting GPqp+2planet to the real data. "nP" means the number of sinusoidal signals injected into the pure noise. I simulate these data sets with the same sampling times as the real data sets. The code for synthetic data is "code/syndata.R". The synthetic data is in "data/synthetic". The parameters used to generate the synthetic data is in "data/par". The real data sets are in "data/real". 

## 2. Test results

By applying different noise models to different synthetic data, I get a series of likelihood/BF tables and plot the (relative) log likelihoods/BFs for different noise models. The log likelihood tables are in '.txt' files in "results/". For example, "like_GJ699_typeGP0P1_sr1.txt" is the likelihood table for GP0P data set. The plots are also in the "results/" folder. I also show the BFs for all data sets and noise models in "BFall.pdf". 


