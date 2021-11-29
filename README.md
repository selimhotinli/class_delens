CLASS delens: Delensing module for the CLASS code 
==============================================

Authors: Selim C. Hotinli, Joel Meyers, Cynthia Trendafilova

Code produces delensed CMB spectra (TT, TE, EE and BB) and lensing-reconstruction noise for given CMB experiment specifications and cosmology. 

Delensing reverses the effects of lensing on the observed CMB temperature and polarization maps. 
This provides various benefits. 
Delensed CMB spectra have sharper acoustic peaks and more prominent damping tails, allowing for improved inferences of cosmological parameters that impact those features.
Delensing reduces B-mode power, aiding the search for primordial gravitational waves and allowing for lower variance reconstruction of lensing and other sources of secondary CMB anisotropies.
Lensing-induced power spectrum covariances are reduced by delensing, simplifying analyses and improving constraints on primordial non-Gaussianities. 
Please refer to [arXiv] for a detailed demonstration of the benefits of CMB delensing.

This code can be used as a submodule for the Fisher forecasting tool software https://github.com/ctrendafilova/FisherLens. 



Using the code
==============================================

In the standard way described in Ref.

```
make class
```

The ```delensed_cmb.ini``` file is the reference input file, analogous to ```explainatory.ini``` in the standard CLASS code. 

```delensed_cmb.ini``` containts and explains the use of all possible additional input parameters related to the delensing implementation. 

New options for delensing implementation
==============================================

If ```delensing = yes``` : the code performs delensing once

If ```delensing = iterative``` : the code performs iterative delensing until spectra converge

If any other option for ```delensing``` parameter is chosen, the code will not calculate the delensed spectra. 

<!-- Selecting CMB noise -->
## Selecting CMB noise

Options ```temperature noise spectra type``` and ```polarization noise spectra type``` determine the type of CMB spectra noise to be considered when calculating the lensing-reconstruction noise and delensed CMB spectra.

If ```temperature noise spectra type = idealized```: Idealized temperature noise used in delensing is calculated by the code internally. Same for the polarization noise.

If ```temperature noise spectra type = external```: Temperature noise is taken externally. Same for the polarization noise.

In the case ```temperature noise spectra type``` is set to ```external```, please define the directory for the noise data with the option

```command_for_temperature_noise_spec  = cat /User/directory/to/CMB/noise.txt```

Same applies for the polarization noise.

If the CMB noise is set to ```idealized```, set the properties of the temperature and polarization noise in radians with the options below (for example for 1 arminute beam and muK-arcminute for noise <img src="http://latex.codecogs.com/svg.latex?\Delta_T" border="0"/>).  
```delta_noise = 0.000290888209```
```sigma_beam  = 0.000290888209```

If you wish to output the temperature/polarization spectra, you can set ```output_spectra_noise = yes```.

<!-- Options on lensing-reconstruction noise calculation -->
## Options on lensing-reconstruction noise calculation

The code calculates the lensing reconstruction noise internally for a given CMB noise with the option ```lensing reconstruction noise spectra type``` set equal to ```internal```.

Internal calculation of the lensing reconstruction noise can be further specilazed depending on the user's choise with the options ```noise_iteration_type``` ```min_varr_type```.

```noise_iteration_type``` defines the set of lensing-reconstruction quadratic estimators (equilalently, normalisations) used in the calculation of the minimum-variance lensing-reconstruction noise if ```delensing``` is set equal to ```iterative```.

You can set ```noise_iteration_type``` equal to ```all``` if you wish to perform iterative delensing on all quadratic estimators.

You can set ```noise_iteration_type``` equal to ```diag``` if you wish to perform iterative delensing only on the diagonal components in the covariance matrix of quadratic estimators.

You can set ```noise_iteration_type``` equal to ```eb``` if you wish to perform iterative delensing only the EB-EB quadratic estimator.

Performing iterative delensing on the estimators improves the quality of delensing by reducing the lensing-reconstruction noise as discussed in [arXiv]. 

```min_varr_type``` option determines how the minimum variance lensing-reconstruction noise is calculated at the end of delensing. 

You can set ```min_varr_type``` equal to ```all``` to include all elements of the quadratic-estimator covariance when calculating the minimum-variance lensing reconstruction noise. 

You can set ```min_varr_type``` equal to ```eb``` to include only the EB-EB quadratic estimator.

You can set ```min_varr_type``` equal to ```diag``` to include only the diagonal components of the quadratic-estimator covariance.

The verbosity of the delensing module can be set with the ```delensing_verbose``` option, similar to default verbosity options in defined CLASS.

If ```delensing = iterative``` you can set the type of the convergence criteria for iterative delensing with the option ```convergence type```.

If ```convergence type = every```  delensing converges if only all the differences between the last and the previous iteration of the lensing noise spectra divided by the current value of the spectra at that mode is less than 'convergence_criterion_itr'.

If ```convergence type = total```  convergence if only the sum of all the differences between the current and the previous iteration in the lensing noise spectra divided by the sum of all the current value of the spectra at that mode is less than 'convergence_criterion_itr'

also, dont forget to choose e.g. ```convergence_criterion_itr = 1e-4```.

If ```lensing reconstruction noise spectra type = external```: Lensing-reconstruction noise is taken externally. In this case please define the path to the lensing-reconstruction noise with ```command_for_lensing_noise_spec  = cat /User/directory/to/Lensing/noise.txt```

<!-- Options regarding the derivatives of the lensed or delensed spectra -->
## Options regarding the derivatives of the lensed or delensed spectra

If you wish to calculate the derivatives of the lensed and delensed spectra with respect to the unlensed spectra, please set

```calculate_derviaties_wrt_unlensed = yes```

and set ```unlensed derivative type``` equal to either ```lensed``` or ```delensed```


If you wish to calculate the derivatives of the lensed and delensed spectra with respect to the lensing power-spectrum, please set

```delensing derivatives = yes```

and set ```derivative type``` equal to either ```lensed``` or ```delensed```

If you wish to output the derivatives, please set ```output_derivatives = yes```.

Note that the matrices including the derivatives can be very large for large values of maximum multipole values included in delensing. 

You can set the sparsity of the produced matrices with the option ```derv_binedges```. 

If ```derv_binedges = 10``` for example, every 10th derivative in L multipole will be saved from the where L is the multipole corresponding to the derivative (either of the lensing potential or the unlensed spectra).

Option ```delta_dl_max``` set the buffer between the lensed spectra and delensed spectra.

Delensing primer
==============================================

Gravitational lensing deflects CMB photons such that the lensed CMB temperature and polarization in line-of-sight direction <img src="http://latex.codecogs.com/svg.latex?\boldsymbol{d}(\boldsymbol{n})" border="0"/> are given by the unlensed CMB in a direction that differs from the line-of-sight direction by the lensing deflection <img src="http://latex.codecogs.com/svg.latex?\boldsymbol{d}(\boldsymbol{n})" border="0"/>.
At lowest order, the deflection angle is a pure gradient <img src="http://latex.codecogs.com/svg.latex?\boldsymbol{d}(\boldsymbol{n})=\boldsymbol{\nabla}\phi(\boldsymbol{n})" border="0"/> where <img src="http://latex.codecogs.com/svg.latex?\phi" border="0"/> is the lensing potential. 
For example, the lensed temperature field <img src="http://latex.codecogs.com/svg.latex?T^{\textrm{lensed}}" border="0"/> is given in terms of the unlensed temperature field <img src="http://latex.codecogs.com/svg.latex?T^{\textrm{unlensed}}" border="0"/> by

<img src="http://latex.codecogs.com/svg.latex?T^{\textrm{lensed}}(\boldsymbol{n})=T^{\textrm{unlensed}}(\boldsymbol{n}+\boldsymbol{d}(\boldsymbol{n}))= T^{\textrm{unlensed}}(\boldsymbol{n}) + \boldsymbol{d}(\boldsymbol{n})\cdot\boldsymbol{\nabla}T^{\textrm{unlensed}}(\boldsymbol{n}) + \ldots \, ." border="0"/>

The aim of delensing is to manipulate observed CMB maps (such as <img src="http://latex.codecogs.com/svg.latex?T^{\textrm{obs}}" border="0"/>) and estimates of the lensing deflection <img src="http://latex.codecogs.com/svg.latex?\boldsymbol{d}^{\textrm{obs}}" border="0"/> to reverse this remapping in order to recover an estimate of the unlensed CMB.

Software implementation
==============================================

We implement our delensing procedure as a modification of the lensing routine in the CLASS Boltzmann code Ref. 
Our goal is to provide a tool that allows accurate, stable, and efficient computation of the delensed CMB spectra and lensing reconstruction noise.
In addition to benefiting from the specialised numerical routines available in the \texttt{CLASS} code, we use efficient real-space expressions introduced in ([Dvorkin et al.](https://arxiv.org/abs/0902.4413), [Smith et al.](https://arxiv.org/abs/1010.0048)) for the delensed spectra and lensing reconstruction quadratic estimators.
These expressions can be evaluated at the cost of <img src="http://latex.codecogs.com/svg.latex?\mathcal{O}(\ell_{\textrm{max}}^2)" border="0"/> rather than <img src="http://latex.codecogs.com/svg.latex?\mathcal{O}(\ell_{\textrm{max}}^3)" border="0"/>, where <img src="http://latex.codecogs.com/svg.latex?\ell_{\textrm{max}}" border="0"/> is the maximum multipole used in the quadratic estimator calculation.
This provides a significant increase in speed which is valuable for repeated calculations of iterative delensing as well as rapid exploration of the parameter space for a Markov chain Monte Carlo analysis, for example.
This simplification is due to computing products of Wigner <img src="http://latex.codecogs.com/svg.latex?3j" border="0"/>-symbols in terms of Wigner <img src="http://latex.codecogs.com/svg.latex?d" border="0"/>-matrices. 


<img src="http://latex.codecogs.com/svg.latex?\mathbf{d}(\boldsymbol{n})" border="0"/>
