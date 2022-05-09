# Bayesian nonparametric estimation of the probability of discovering new species

[![CI][CI-img]][CI-url] [![CODECOV][codecov-img]][codecov-url]

This repository contains a Julia implementation of the simulations described in [*Bayesian nonparametric estimation of the probability of discovering new species* by Lijoi, Mena, and Prünster (2007)](https://www.jstor.org/stable/20441417?seq=1#metadata_info_tab_contents). A deck of slides used to present this project can be found [here](https://github.com/scortino/lijoi07-slides).

## Simulations

The simulations described both in the paper and in the slides are implemented as unit tests for the `SpeciesBNP` package defined in `src/SpeciesBNP.jl`.

As such, they can be run from the project directory as follows:

```console
git clone https://github.com/scortino/lijoi07.git
cd lijoi07
julia --project=.
using Pkg; Pkg.test()
```

This generates some summary graphs for the simulations that are saved in `img/`, if not already present.

## Credits

This project was proposed by Professors C. Feinauer, I. Prünster, and G. Zanella as part of their course [20605 - Machine Learning 2](https://didattica.unibocconi.eu/ts/tsn_anteprima.php?cod_ins=20605&anno=2022&IdPag=).

[CI-img]: https://github.com/scortino/lijoi07/actions/workflows/ci.yml/badge.svg
[CI-url]: https://github.com/scortino/lijoi07/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/scortino/lijoi07/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/scortino/lijoi07