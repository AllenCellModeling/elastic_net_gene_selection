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
* Turn your project into a GitHub repository:
  * Make sure you have `git` installed, if you don't, [follow these instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
  * Make an account on [github.com](https://github.com)
  * Go to [make a new repository](https://github.com/new)
  * _Recommendations:_
    * _It is strongly recommended to make the repository name the same as the Python
    package name_
    * _A lot of the following optional steps are *free* if the repository is Public,
    plus open source is cool_
  * After a GitHub repo has been created, run the following commands:
    * `git remote add origin git@github.com:AllenCellModeling/elastic_net_gene_selection.git`
    * `git push -u origin master`
* Register elastic_net_gene_selection with Codecov:
  * Make an account on [codecov.io](https://codecov.io)
  (Recommended to sign in with GitHub)
  * Select `AllenCellModeling` and click: `Add new repository`
  * Copy the token provided, go to your [GitHub repository's settings and under the `Secrets` tab](https://github.com/AllenCellModeling/elastic_net_gene_selection/settings/secrets),
  add a secret called `CODECOV_TOKEN` with the token you just copied.
  Don't worry, no one will see this token because it will be encrypted.
* Generate and add an access token as a secret to the repository for auto documentation
generation to work
  * Go to your [GitHub account's Personal Access Tokens page](https://github.com/settings/tokens)
  * Click: `Generate new token`
  * _Recommendations:_
    * _Name the token: "Auto-Documentation Generation" or similar so you know what it
    is being used for later_
    * _Select only: `repo:status`, `repo_deployment`, and `public_repo` to limit what
    this token has access to_
  * Copy the newly generated token
  * Go to your [GitHub repository's settings and under the `Secrets` tab](https://github.com/AllenCellModeling/elastic_net_gene_selection/settings/secrets),
  add a secret called `ACCESS_TOKEN` with the personal access token you just created.
  Don't worry, no one will see this password because it will be encrypted.
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
* Add branch protections to `master` and `stable`
    * To protect from just anyone pushing to `master` or `stable` (the branches with
    more tests and deploy
    configurations)
    * Go to your [GitHub repository's settings and under the `Branches` tab](https://github.com/AllenCellModeling/elastic_net_gene_selection/settings/branches), click `Add rule` and select the
    settings you believe best.
    * _Recommendations:_
      * _Require pull request reviews before merging_
      * _Require status checks to pass before merging (Recommended: lint and test)_

***Free software: Allen Institute Software License***

