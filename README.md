# Gaussian-Mixture-Based Enhanced Sampling

## GAMBES

In this enhanced sampling method the bias potential <a href="https://www.codecogs.com/eqnedit.php?latex=V({\boldsymbol{R}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V({\boldsymbol{R}})" title="V({\boldsymbol{R}})" /></a> is estimated from a model probability distribution.
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=V({\boldsymbol{R}})&space;=&space;\frac{1}{\beta}&space;\log&space;\Big(&space;\sum\limits_{i=1}^M&space;\frac{Z^1}{Z^i}&space;\&space;p^i(\boldsymbol{d({\boldsymbol{R}})})&space;\Big)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V({\boldsymbol{R}})&space;=&space;\frac{1}{\beta}&space;\log&space;\Big(&space;\sum\limits_{i=1}^M&space;\frac{Z^1}{Z^i}&space;\&space;p^i(\boldsymbol{d({\boldsymbol{R}})})&space;\Big)" title="V({\boldsymbol{R}}) = \frac{1}{\beta} \log \Big( \sum\limits_{i=1}^M \frac{Z^1}{Z^i} \ p^i(\boldsymbol{d({\boldsymbol{R}})}) \Big)" /></a>
</p>
where the probability of the entire space is estimated as a mixture of probability densities of metastable islands <a href="https://www.codecogs.com/eqnedit.php?latex=p^i(\boldsymbol{d})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p^i(\boldsymbol{d})" title="p^i(\boldsymbol{d})" /></a>. The probability density of each of these islands are then estimated using a Gaussian Mixture Model from the data generated using unbiased molecular dynamics trajectories. 
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=p^i(\boldsymbol{d})&space;\cong&space;\sum\limits_{k=1}^{K^i}&space;\pi_k^i&space;\&space;\mathcal{N}(\boldsymbol{d}|\boldsymbol{\mu}_k^i,\boldsymbol{\Sigma}_k^i)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p^i(\boldsymbol{d})&space;\cong&space;\sum\limits_{k=1}^{K^i}&space;\pi_k^i&space;\&space;\mathcal{N}(\boldsymbol{d}|\boldsymbol{\mu}_k^i,\boldsymbol{\Sigma}_k^i)" title="p^i(\boldsymbol{d}) \cong \sum\limits_{k=1}^{K^i} \pi_k^i \ \mathcal{N}(\boldsymbol{d}|\boldsymbol{\mu}_k^i,\boldsymbol{\Sigma}_k^i)" /></a> 
</p>

The relative height of each island is estimated on-the-fly using: 
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{Z^1}{Z^i}&space;=&space;\frac{\Big<&space;p^1(\boldsymbol{d({\boldsymbol{R}})})&space;e^{-\beta(U({\boldsymbol{R}})&plus;V({\boldsymbol{R}}))}\Big>_{V}}{\Big<&space;p^i(\boldsymbol{d({\boldsymbol{R}})})&space;e^{-\beta(U({\boldsymbol{R}})&plus;V({\boldsymbol{R}}))}\Big>_{V}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{Z^1}{Z^i}&space;=&space;\frac{\Big<&space;p^1(\boldsymbol{d({\boldsymbol{R}})})&space;e^{-\beta(U({\boldsymbol{R}})&plus;V({\boldsymbol{R}}))}\Big>_{V}}{\Big<&space;p^i(\boldsymbol{d({\boldsymbol{R}})})&space;e^{-\beta(U({\boldsymbol{R}})&plus;V({\boldsymbol{R}}))}\Big>_{V}}" title="\frac{Z^1}{Z^i} = \frac{\Big< p^1(\boldsymbol{d({\boldsymbol{R}})}) e^{-\beta(U({\boldsymbol{R}})+V({\boldsymbol{R}}))}\Big>_{V}}{\Big< p^i(\boldsymbol{d({\boldsymbol{R}})}) e^{-\beta(U({\boldsymbol{R}})+V({\boldsymbol{R}}))}\Big>_{V}}" </a>
</p>
More details about the method can be found at: https://arxiv.org/abs/1909.07773 or https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01125 

The code for this algorithm can be found at GAMBES.cpp , it is compatible with PLUMED v2.7-dev. 

## Running GAMBES for your system
<ul>
<li> First run unbiased molecular dynamics simulations from at least one known state and for all known states. </li>
<li> Using the GAMBES-module.ipynb code, fit a Gaussian Mixture Model to your state using any number of descriptors. </li>
<li> The final file generated using the notebook is an input for GAMBES in PLUMED. Generate similar files for every known state of interest and name them as < file-name >.0 ,  < file-name >.1, ... , < file-name >.m where 0,1,...,m are your states. This < file-name > is an input for GAMBES. </li>
<li> Copy GAMBES.cpp in a folder containing your plumed input, the states generated < file-name >.* and all other necessary files for running your system. </li>
<li> Run the simulation. </li>
 </ul>
