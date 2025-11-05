# Surrogate Strategies for Scalarisation-based Multi-objective Bayesian Optimizers

## Introduction

Source code for reproducing the results of the manuscript titled: **Surrogate Strategies for Scalarisation-based Multi-objective Bayesian Optimizers**. 

This code runs the Bayesian multi-objective optimisation algorithms known as ParEGO and MParEGO on the BBOB bi-objective "benchmarking" test problems. These algorithms are in the C++ Tigon optimisation library, which is part of the Liger software. The BBOB bi-objective problems are part of the COCO library, which is also supplied with Liger. This repository C++ code configures the above optimisation algorithms, and runs them on the BBOB bi-objective problems. There is also two design of experiment techniques that are used for comparison, namely Random Search and Sobol.

For comparing the optimisation algorithms this repository follows the performance assessment of the COCO platform. There  provides a Python implementation of the Empirical ECDF   


## 



## Supplementary file

A supplementary file with all the results is available in https://jduro.ddns.net/EMO2025/SupQingyuEMO2025.pdf.

## Citation

In case you use any of this code in your work, please cite the following paper:

Mo, Q.; Duro, J. A.; and Purshouse, R. C., Surrogate Strategies for Scalarisation-based Multi-objective Bayesian Optimizers, In: Singh, H., et al. Evolutionary Multi-Criterion Optimization. EMO 2025. Lecture Notes in Computer Science, vol 15513. Springer, Singapore. https://doi.org/10.1007/978-981-96-3538-2_10 
