---
title: Introduction to WHAM
date: 2018-08-08 12:00:00
categories: reaserch notes
toc: true
mathjax: true
comments: true
---
This is a blog for WHAM introduction. And it's also a notes for seminar of our lab 8.8. Most of this note will focus on mathematical background of WHAM.
<!--more-->

# Introduction to WHAM 



## Umbrella Sampling

### Partition function and Helmholtz free energy

​	The canonical partition function $Q$ of a system can be calculated via an integral over the whole phase space.

$$
Q = \int exp(-\beta E(r))d^Nr
$$

where $\beta = \frac{1}{k_BT}$ . The free (Helmholtz) energy $A$ is related to $Q$ via 
$$
A = -1/\beta \ln Q
$$
​	In many cases, a reaction coordinate ($\xi$), a continuous parameter which provide distinction between two thermodynamic states, can be defined. With $\xi$ defined, the probability distribution of the system along $\xi$ can be calculated by integrating out all degrees of freedom but $\xi$  :
$$
Q(\xi) = \frac{\int \delta[\xi(r)-\xi)]exp(-\beta E)d^Nr}{\int exp(-\beta E)d^Nr}
$$
$Q(\xi)d\xi$ can be interpreted as the probability of finding the system in a small interval $d\xi$ around $\xi$.
$$
A(\xi) = -1/\beta\ln Q(\xi)
$$
$A(\xi)$ is also called potential of mean force (PMF).

​	In computer simulations, the direct phase-space integrals used in Eqs.(3) and (1) are impossible to calculate. However, if the system is ergodic, i.e. if every point in phase space is visited during the simulation, $Q(\xi)$ is equal to 
$$
P(\xi) = \lim_{t\to\infty}\frac{1}{t}\int_{0}^{t}\rho[\xi(t')]dt'
$$
​	where $t$ denotes the time and $\rho$ simply counts the occurrence of $\xi$ in a given interval (of infinitesimal width in the exact equation and of finite width when calculating a histogram). So, in principle, $A(\xi)$ can be directly obtained from molecular dynamics (MD) simulations by monitoring $P(\xi)$ , the distribution of the system along the reaction coordinate.

### UMBRELLA SAMPLING: METHOD

​	A bias, an additional energy term, is applied to the system to ensure efficient sampling along the whole reaction coordinate. The effect of the bias potential to connect energetically separated regions in phase space gave rise to the name umbrella sampling.

​	The bias potential $w_i$ of windows $i$ is an additional energy term, which depends only on the reaction coordinate:
$$
E^b(r) = E^u(r)+w_i(\xi)
$$
​	According to Eq.(3), the unbiased distribution is
$$
P_i^u(\xi_0) = \frac{\int \delta[\xi(r)-\xi_0)]exp(-\beta E(r))d^Nr}{\int exp(-\beta E(r))d^Nr}
$$
​	MD simulation of the biased system provides the biased distribution along the reaction coordinate $P_i^b$.
$$
P_i^b(\xi_0) = \frac{\int \delta[\xi(r)-\xi_0)]exp(-\beta [E(r)+w_i(\xi(r))])d^Nr}{\int exp(-\beta [E(r)+w_i(\xi(r))])d^Nr}
$$
​	Because the bias energy $w_i(\xi)$ depends only on $\xi$ , and the integration in the enumerator is performed over all degrees of freedom but $\xi$ ,
$$
P_i^b(\xi_0) = exp[-\beta w_i(\xi_0)]\frac{\int \delta[\xi(r)-\xi_0)]exp(-\beta E(r))d^Nr}{\int exp(-\beta [E(r)+w_i(\xi(r))])d^Nr}
$$
​	Compared with Eq.(7), 
$$
\begin{split}
P_i^u(\xi_0)&= P_i^b(\xi_0)exp(\beta w_i(\xi_0))\cdot\frac{\int exp(-\beta E(r)) exp[-\beta w_i(\xi_0(r))] d^Nr}{\int exp(-\beta E(r))d^Nr}\\&=P_i^b(\xi_0) exp(\beta w_i(\xi_0)) \langle exp[-\beta w_i(\xi_0)]\rangle 
\end{split}\tag{10}
$$
​	From equation above, $A_i(\xi)$ can be readily evaluated,
$$
A_i(\xi_0) = -(1/\beta)\ln P_i^b(\xi_0)-w_i(\xi_0)+F_i \tag{11}
$$
where $F_i = -(1/\beta)\ln(\langle  exp[-\beta w_i(\xi_0)] \rangle)$ 

​	If the free-energy curve $A_i(\xi)$ of more windows are to be combined to one global $A(\xi)$, the $F_i$ have to be calculated. They are associated with introducing the bias potential and connect the free-energy curves $A_i(\xi)$ obtained in the different windows:
$$
\begin{split}
exp(-\beta F_i)&=\langle exp(-\beta w_i(\xi)) \rangle \\&=\int P^u(\xi)exp[-\beta w_i(\xi)]d\xi \\&=
\int exp(-\beta [A(\xi)+w_i(\xi)]) d\xi
\end{split} \tag{12}
$$


## Weighted Histogram Analysis Method

### Mathematical Analyze

​	Numerous methods have been proposed for an estimation of $F_i$,  a promising one being the ***WHAM***. It aims to minimize the statistical error of $P^u(\xi)$ . The global distribution is calculated by weighted average of the distributions of the individual windows: 
$$
P^u(\xi) = \sum_{i}^{windows}p_iP_i^u(\xi) \tag{13}
$$
​	The weights $p_i$ are chosen in order to minimize the statistical error of $P^u$:
$$
\frac{\partial \sigma^2(P^u) }{\partial p_i}=0
$$
under the condition $\sum p_i=1$. This minimization problem can be solved with Lagrange multiplier method.
$$
\mathcal{L}(p_i,\mu) = \sum_{i=1}^{N} p_i^2 \sigma^2[P_i^u(\xi)] - \mu(\sum_{i=1}^{N}p_i-1)
$$
And this leads to 
$$
p_i = \frac{(\sigma^2[P_i^u(\xi)])^{-1}}{\sum_{j=1}^{N}\sigma^2[P_j^u(\xi)])^{-1}}=
\frac{n_i exp(-\beta[w_i(\xi)-F_i])}{\sum_{j=1}^{N}n_j exp(-\beta[w_j(\xi)-F_j])}
$$
​	According to Umbrella Sampling,
$$
P_i^u(\xi) = exp(\beta[w_i(\xi)-F_i])P_i^b(\xi)
$$
​	The $P^u(\xi)$ takes this form, 
$$
P^u(\xi) = \sum_{i=1}^{N}\frac{n_i}{\sum_{j=1}^{N} n_j e^{-\beta[w_j(\xi)-F_j]}}P_i^b(\xi).
$$
 where $n_i$ being the total number of steps sampled for window $i$. $N$ is the number of windows. $P_i^b(\xi)$ is the actual distribution of configurations occurring during the $i$th of these simulations. It may be formally expressed as,
$$
P_i^b(\xi)=\frac{1}{n_i}\sum_{l=1}^{n_i}\delta(\xi-\xi_{i,l})
$$
​	The free energies corresponding to the different biasing potentials can be computed using the distribution $P^u(\xi)$ , thus,
$$
\begin{split}
e^{-\beta F_k}&=\int d\xi P^u(\xi)e^{-\beta w_k(\xi)}\\
&=\sum_{i=1}^{N}\sum_{l=1}^{n_i}\frac{e^{-\beta w_k(\xi_{i,l})}}{\sum_{j=1}^{N} n_je^{-\beta[w_j(\xi_{i,l})-F_j]}}
\end{split}
$$
​	The $\{ F_i\}$ can be computed iteratively.  	


## How to use WHAM

You can refer to [WHAM](http://membrane.urmc.rochester.edu/content/wham). It's a very popular version of WHAM implementation.
