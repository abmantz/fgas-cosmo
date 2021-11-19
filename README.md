# fgas-cosmo

Cosmological likelihood code for galaxy cluster gas mass fractions

## Contributors

* Adam Mantz
* Anja von der Linden
* Douglas Applegate
* Lucie Baumont

## Background

Massive, dynamically relaxed clusters of galaxies represent a minority of the cluster population, and are a precious resource for probing both the astrophysics of clusters and cosmology. They provide precise and minimally biased mass estimates to studies of scaling relations and the growth of cosmic structure, and enable complementary contraints on cosmic expansion. Additionally, they represent the most natural targets for studying thermodynamic features of the intracluster medium with minimal systematics from deprojection or non-equilibrium processes.

## Contents

This repo holds the likelihood code and data for galaxy cluster gas-mass fractions, as used to constrain cosmological models in [this paper](https://arxiv.org/abs/2111.09343), which should be cited if you use this code or the associated data in scientific research. The code is an evolution of the one previously released along with [this earlier paper](https://arxiv.org/abs/1402.6212).

The code exists as a stand-alone library written in C++, with some additional code to make it callable from Fortran or C. We also provide supporting files and instructions for integrating this code with [CosmoMC](https://cosmologist.info/cosmomc/) (which work as of this writing).

An outdated version of this code is included in [Cosmosis](https://bitbucket.org/joezuntz/cosmosis/) under the module named "cluster_fgas". Further development would be needed to restore compatibility with the current version (or to provide a Python interface to it generally.)

## More documentation

* [cosmomc/README.md](cosmomc/README.md): getting set up with CosmoMC
* [src/README.md](src/README.md): which code file does what

## License

All contents Copyright 2020 Adam Mantz and the contributors listed above, and licensed under the BSD 3-Clause License, unless otherwise noted.
