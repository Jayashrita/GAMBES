# Gaussian-Mixture-Based Enhanced Sampling

## GAMBES

In this enhanced sampling method the bias potential is estimated from a model probability distribution (here, a Gaussian mixture model).

<a href="https://www.codecogs.com/eqnedit.php?latex=V({\boldsymbol{R}})&space;=&space;\frac{1}{\beta}&space;\log&space;\Big(&space;\sum\limits_{i=1}^M&space;\frac{Z^1}{Z^i}&space;\&space;p^i(\boldsymbol{d({\boldsymbol{R}})})&space;\Big)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?V({\boldsymbol{R}})&space;=&space;\frac{1}{\beta}&space;\log&space;\Big(&space;\sum\limits_{i=1}^M&space;\frac{Z^1}{Z^i}&space;\&space;p^i(\boldsymbol{d({\boldsymbol{R}})})&space;\Big)" title="V({\boldsymbol{R}}) = \frac{1}{\beta} \log \Big( \sum\limits_{i=1}^M \frac{Z^1}{Z^i} \ p^i(\boldsymbol{d({\boldsymbol{R}})}) \Big)" /></a>

where 



More details about the method can be found at: https://arxiv.org/abs/1909.07773  

The code for this algorithm can be found at GAMBES.cpp , it is compatible with PLUMED v2.6. 
