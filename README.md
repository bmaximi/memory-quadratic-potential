# memory-quadratic-potential
This repository contains the implementation of the algorithm presented in article [1], which is based on Method B from article [2].

# How to run the code
In order to run the code, you first need to download and install the cvx package from https://cvxr.com/cvx/download/.
For a detailed instruction, we refer to https://cvxr.com/cvx/doc/install.html. 
We used Version 2.2, Build 1148 from 2020 for our numerical tests.

Once the cvx package is properly installed, you can run the example1Dv.m, example2D.m and example1Dp.m.
The scripts compute the Markovian approximations for the 3 examples provided in [1] and generates the corresponding figures, respectively.

If you want to use this code for your own project, we highly recommend reading the articles [1,2] first to get a better understanding
of the underlying theory, since the algorithm might fail due to improperly chosen parameters.

# References
[1] Braun M., Wolf N. and Hanke M., On generalized Langevin dynamics in the presence of a quadratic potential, manuscript in preparation.
[2] Bockius N., Braun M., Hofmann K., Schmid F. and Hanke M., Determining extended Markov parameterizations for vector-valued
generalized Langevin equations, Z. Naturforsch. A, vol. 81, no. 2 (2026), pp. 151--171. (DOI: 10.1515/zna-2025-0354)
