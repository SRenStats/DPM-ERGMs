# DPM-ERGMs
Ensembles of networks are increasingly avaliable, where multiple independent networks are observed on the same set of nodes, for example, a collection of brain networks constructed on the same brain regions for different individuals.

We build a Dirichlet Process Mixture of Exponential Random Graph Models (DPM-ERGMs) to describe both the variations and characteristics of networks in an ensemble at the same time.

In order to perform full Bayesian inference for DPM-ERGMs, we employ the intermediate importance sampling technique inside the Metropolis-within-slice sampling scheme, which addressed the problem of sampling from the intractable ERGMs on an infinite sample space.
