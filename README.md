# Surrogate Strategies for Scalarisation-based Multi-objective Bayesian Optimizers

## Introduction

Source code for reproducing the results of the manuscript titled: **Surrogate Strategies for Scalarisation-based Multi-objective Bayesian Optimizers** (full reference below). 

The C++ code in this repository runs two Bayesian multi-objective optimisation algorithms, known as ParEGO and MParEGO, on benchmark problems from the bbob-biobj test suite, for a total of 55 bi-objective problems. The implementation of the optimisation algorithms is provided by the C++ Tigon optimisation library, which is part of the Liger software (see https://ligerdev.shef.ac.uk/liger-development-team/liger-ce). The implementation of the bi-objective bbobj test suite is provided by the C/C++ COCO library (see https://github.com/numbbo/coco). Besides ParEGO and MParEGO, there is also two design of experiment (DoE) techniques that can be used for comparison, namely Random Search and Sobol. The user is able to select: 
1. The optimisation algorithm, and the options are: Random, Sobol, ParEGO, and MParEGO. 
2. The optimisation budget (i.e. the number of function evaluations); 
3. The bbob-biobj problem (the available problems range from F1-F55);
4. The number of decision variables (or inputs), and; 
5. The number of samples for Monte Carlo estimation. The option only applies to MParEGO.  

This repository also provides a Python jupyter notebook for computing the target-based Empirical Cumulative Distribution Function (ECDF), which is used to assess the performance of the optimisation algorithms. The same jupyter notebook generates the cartesian plots that show the "anytime" performance of the optimisation algorithms, and this follows the same format suggested by the COCO library.


## Compilation and Run



```shell
make
```



## Supplementary file

A supplementary file with all the results is available in https://jduro.ddns.net/EMO2025/SupQingyuEMO2025.pdf.

## License

This software is Copyright (C) 2025 The University of Sheffield (www.sheffield.ac.uk).

The software is free software (software libre); you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License (LGPL) version 3 as published by the Free Software Foundation. Please review the following information to ensure the GNU Lesser General Public License version 3 requirements will be met: https://www.gnu.org/licenses/lgpl-3.0.html.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Although this software is distributed as free software,  it does not excuse you from scientific propriety, which obligates you to give appropriate credit! If you write a scientific paper describing research that made substantive use of this program, it is your obligation as a scientist to (a) mention the fashion in which this software was used in the Methods section; (b) mention the algorithm in the References section. The appropriate citation is:

Mo, Q.; Duro, J. A.; and Purshouse, R. C., Surrogate Strategies for Scalarisation-based Multi-objective Bayesian Optimizers, In: Singh, H., et al. Evolutionary Multi-Criterion Optimization. EMO 2025. Lecture Notes in Computer Science, vol 15513. Springer, Singapore. https://doi.org/10.1007/978-981-96-3538-2_10 
