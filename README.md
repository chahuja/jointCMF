jointCMF
========

Complex Matrix Factorization  

nmfStandard.m
============
Standard NMF Interations to minimize Eucledian distance with parameters as 
- nIter -> Number of Iterations
- fixedX -> to keep X fixed while interating
- X -> to initialize X with given values


nmfWithPhase4.m
===============
Mis-named as nmfWithPhase, it is basically performing complex matrix factorization based on the paper http://arxiv-web3.library.cornell.edu/pdf/1411.6741v1.pdf
It has a single dependency on nmfStandard.m which it uses to perform Complex Matrix Factorization

- Outputs
  - X -> Complex Factorized Vector Sets  (PxQ)
  - H -> Real Factorized Weights (QxR)
  - Xc -> Real Matrix which can be appropriatly transformed to get X (2Px2Q) 
