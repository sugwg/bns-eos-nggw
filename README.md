# Measuring Neutron Star Radius with second and third generation Gravitational Wave Detector Networks

**Ananya Bandopadhyay<sup>1</sup>, Keisi Kacanja<sup>1</sup>, Rahul Somasundaram<sup>1,2</sup>, Alexander H. Nitz<sup>1</sup>, Duncan A. Brown<sup>1</sup>**

**<sup>1</sup>Department of Physics, Syracuse University, Syracuse, NY 13244, USA**

**<sup>2</sup>Theoretical Division, Los Alamos National Lab, Los Alamos, NM 87545, USA**

This repository is a companion data release. The preprint version of the paper is available on [arXiv](https://arxiv.org/abs/2402.05056).

## License

![Creative Commons License](https://i.creativecommons.org/l/by-sa/3.0/us/88x31.png "Creative Commons License")

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 United States License](http://creativecommons.org/licenses/by-sa/3.0/us/). 

## Abstract

The next generation of ground-based interferometric gravitational wave detectors will observe mergers of black holes and neutron stars throughout cosmic time. A large number of the binary neutron star merger events will be observed with extreme high fidelity, and will provide stringent constraints on the equation of state of nuclear matter. In this paper, we investigate the systematic improvement in the measurability of the equation of state with increase in detector sensitivity by combining constraints obtained on the radius of a $1.4 M_{\odot}$ neutron star from a simulated source population. Since the measurability of the equation of state (EoS) depends on its stiffness, we consider a range of realistic equations of state that span the current observational constraints. We show that a single 40km Cosmic Explorer detector can pin down the neutron star radius for a soft, medium and stiff equation of state to an accuracy of 10m within a decade, whereas the current generation of ground-based detectors like the Advanced LIGO-Virgo network would take $\mathcal{O}(10^5)$ years to do so for a soft equation of state. 

## Summary of Contents

Following is a summary of the scripts, data and other supplementary material included in this data release, organized according in directories as described.

* `configurations` : Directory containing base configuration files (in `.ini` format) for creating injections, running inference jobs, and describing custom detector configurations. The `uniform_radius_nsat.cache` ASCII file included within this directory is read during parameter-estimation runs, for sampling over the equations of state.

* `psds` : Directory containing noise curves used for custom detectors, $\rm A^{\sharp}$ (taken from https://dcc.ligo.org/LIGO-T2300041/public), a 20km Cosmic Explorer (CE20) and a 40km Cosmic Explorer (CE40) (taken from  https://dcc.cosmicexplorer.org/CE-T2000017).

* `posterior_data` : Directory containing compressed files for the inferred posterior distributions for the radius of a $1.4 M_{\odot}$ neutron star, used for the EoS-inference analysis. 

* `scripts` : Directory containing scripts to create injections, perform the Bayesian-inference analysis and reproduce figures from the paper. 

* `plots` : Directory containing figures from the paper in pdf format.

* `eos_nsat.tar.gz` : Compressed EoS data files from [Capano et al. 2020](https://doi.org/10.1038/s41550-020-1014-6) containing files `1.dat` through `2000.dat`, for 2000 equations of state (mass-radius-tidal deformability curves) calibrated with chiral EFT upto nuclear saturation density. 

## Funding

This work is supported by NSF awards PHY-2011655(AB, DAB), PHY-2207264(DAB), PHY-2309240 (AHN, KK, DAB) and PHY-2116686(RS). Computations were supported through computational resources provided by Syracuse University. 
