# DISCA


[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

R package of DISCA (Distance-Based Independence Screening for Canonical Analysis) functions.

The DISCA package contains back-elimination and forward-selection version of independence screening procedure for simultaneous dimension reduction, based on distance covariance.

Permutation test of unbiased distance covariance, which is not included in 'energy' package, is also contained in DISCA.

This repository contains:

1.  Separate .R codes of functions constructed in DISCA;
2.  Examples and explanations of each R function;

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
  - [maxDCOV.R / minDCOV.R](#optdcov)
  - [dcovU.test.R](#dcovu.test)
  - [DISCA.R / DISCA_fwd.R](#disca)
- [Maintainers](#maintainers)
- [Contributors](#contributors)

## Backgroud

The DISCA codes started with the issues posted on [@ChuanpingYu](https://github.com/ChuanpingYu) over at [DISCA paper](https://arxiv.org/abs/1903.00037), about how to inspect relationships between a pair of random variables. Distance covariance, an independence testing metric which is valid for random variables with unequal dimensions and arbitrary distributions, is utilized to detect central subspace of 'responses' and 'predictors' simultaneously. For technical details and more explanations of DISCA method, please refer to *Distance-Based Independence Screening for Canonical Analysis*.

The goals of this repository are:

1. Given explanations of each function utilized in DISCA. Examples for each function can be found in the [man](DISCA/man) folder;
2. Provide the R-package of DISCA to be installed;

## Install

**DISCA R-package** contains all function needed to perform DISCA. Each function in the [R](DISCA/R) folder could be called.

```sh
devtools::install_github("Yijin911/DISCA")
```

## Usage

This is a short description of each function constructed in DISCA.

### maxDCOV.R / minDCOV.R

Provided with sample dataset **X** and **Y**, [maxDCOV.R](DISCA/R/maxDOCV.R) and [minDCOV.R](DISCA/R/minDCOV.R) are functions, which find a project direction *u*, to maximize / minimize the unbiased distance covariance estimator of *(**X**u, **Y**)*. [DC algorithm](https://link.springer.com/article/10.1007/s10479-004-5022-1) and [ADMM algorithm](https://www.stat.cmu.edu/~ryantibs/convexopt/lectures/admm.pdf) are utilized to complete optimization.

Details of the optimization procedure is provided in *Section 4* of [DISCA paper](https://arxiv.org/abs/1903.00037), in which [gMatrix.R](DISCA/R/gMatrix.R), [Mminus.R](DISCA/R/Mminus.R), [Mplus.R](DISCA/R/Mplus.R), and [Xdiff.R](DISCA/R/Xdiff.R) are functions used to decompose empirical distance covariance.

### dcovU.test.R

For unbiased estimator of distance covariance, [dcovU.test.R](DISCA/R/dcovU.test.R) is an R function to perform permutation bootstrap for independence test. Number of replicates in the test is self-defined. The default value is 1000.

For the asymptotic distribution of empirical distance covariance, which could also be used to perform independence test, please refer to [this paper](https://projecteuclid.org/journals/annals-of-statistics/volume-35/issue-6/Measuring-and-testing-dependence-by-correlation-of-distances/10.1214/009053607000000505.full).

### DISCA.R / DISCA_fwd.R

These are main functions constructed in DISCA. Provide observations **X** and **Y**, [DISCA.R](DISCA/R/DISCA.R) get corresponding central subspace with back elimination method, while [DISCA_fwd.R](DISCA/R/DISCA_fwd.R) reached the result with a forward selection procedure. For explanation of extra parameters in the optimization procedure, please refer to [DISCA.Rd](DISCA/man/DISCA.Rd) and [DISCA_fwd.Rd](DISCA/man/DISCA_fwd.Rd).

## Maintainers

[@Yijin911](https://github.com/Yijin911).

## Contributors

This project exists thanks to all the people who contribute.
- [Chuanping Yu](https://chuanpingyu.github.io/);
- [Yijin Ni](https://github.com/Yijin911);
- [Andy Ko](https://www.linkedin.com/in/andy-ko-b86b2551);
- [Xiaoming Huo](https://sites.gatech.edu/xiaoming-huo/).
