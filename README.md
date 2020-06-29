# Elastic Net Gene Selection

[![Build Status](https://github.com/AllenCellModeling/elastic_net_gene_selection/workflows/Build%20Master/badge.svg)](https://github.com/AllenCellModeling/elastic_net_gene_selection/actions)
[![Documentation](https://github.com/AllenCellModeling/elastic_net_gene_selection/workflows/Documentation/badge.svg)](https://AllenCellModeling.github.io/elastic_net_gene_selection)
[![Code Coverage](https://codecov.io/gh/AllenCellModeling/elastic_net_gene_selection/branch/master/graph/badge.svg)](https://codecov.io/gh/AllenCellModeling/elastic_net_gene_selection)

scRNA-seq gene selection code

---

## Installation

First create a conda environment with
```
conda create -n enetgsel python=3.8
conda activate enetgsel
conda install -c conda-forge fortran-compiler
```

Then:

**Stable Release:** `pip install elastic_net_gene_selection`<br>
**Development Head:** `pip install git+https://github.com/AllenCellModeling/elastic_net_gene_selection.git`

## Documentation
For full package documentation please visit [AllenCellModeling.github.io/elastic_net_gene_selection](https://AllenCellModeling.github.io/elastic_net_gene_selection).

## Development
See [CONTRIBUTING.md](CONTRIBUTING.md) for information related to developing the code.

## The Four Commands You Need To Know
1. `pip install -e .[dev]`

    This will install your package in editable mode with all the required development
    dependencies (i.e. `tox`).

2. `make build`

    This will run `tox` which will run all your tests in both Python 3.6, Python 3.7,
    and Python 3.8 as well as linting your code.

3. `make clean`

    This will clean up various Python and build generated files so that you can ensure
    that you are working in a clean environment.

4. `make docs`

    This will generate and launch a web browser to view the most up-to-date
    documentation for your Python package.


## TODO: Additional Optional Setup Steps:

* Register your project with PyPI:
  * Make an account on [pypi.org](https://pypi.org)
  * Go to your [GitHub repository's settings and under the `Secrets` tab](https://github.com/AllenCellModeling/elastic_net_gene_selection/settings/secrets),
  add a secret called `PYPI_TOKEN` with your password for your PyPI account.
  Don't worry, no one will see this password because it will be encrypted.
  * Next time you push to the branch: `stable`, GitHub actions will build and deploy
  your Python package to PyPI.
  * _Recommendation: Prior to pushing to `stable` it is recommended to install and run
  `bumpversion` as this will,
  tag a git commit for release and update the `setup.py` version number._

***Free software: Allen Institute Software License***

