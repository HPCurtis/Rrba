Brms and cmdstanr RBA

This work presents direct implementations of AFNI RBA analysis, as discussed in (Chen et al., 2019), utilizing brms and cmdstanr along with tidyverse packages and principles. This approach demonstrates that the correct use of tools can lead to significant speedups in big data problems using Bayesian statistics, such as is the case with fMRI. It highlights how a little effort and thoughtful consideration can yield positive results, in contrast to relying on automated tools that may save time initially but not always in the long term.

Note: I am a massive fan of Paul-Christian Bürkner's brms package. Its utilities for outputting Stan code were crucial in identifying these speedups, especially with the reduce_sum functionality.

# References

Chen, G., Xiao, Y., Taylor, P. A., Rajendra, J. K., Riggins, T., Geng, F., ... & Cox, R. W. (2019). Handling multiplicity in neuroimaging through Bayesian lenses with multilevel modeling. Neuroinformatics, 17, 515-545.