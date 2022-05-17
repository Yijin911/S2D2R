# DISCA


[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

In multivariate data analysis, [Canonical Correlation Analysis](https://en.wikipedia.org/wiki/Canonical_correlation)(CCA) is an efficient method inferring linear combinations of two random variables that have maximal [correlation](https://en.wikipedia.org/wiki/Correlation) with each other. *[Distance-Based Independence Screening for Canonical Analysis](https://arxiv.org/abs/1903.00037)(DISCA)* proposed an algorithm that not only infers linear combinations maximizing distance covariance, but also reduces the dimension of two random vectors.

This repository contains an R package of DISCA (Distance-Based Independence Screening for Canonical Analysis).
We provide both a back-elimination and a forward-selection version of an independence screening procedure for dimension reduction, which is based on the [distance covariance](https://en.wikipedia.org/wiki/Distance_correlation).

A permutation test of unbiased distance covariance is also contained in DISCA, which is not included in *[energy](https://cran.r-project.org/web/packages/energy/index.html)* .

## Table of Contents

- [Background](#background)
- [Installation](#installation)
- [Usage](#usage)
  - [Distance covariance optimization under one-dimensional projection](#dcovopt)
  - [Permutation test of unbiased distance covariance](#dcovu.test)
  - [DISCA functions](#disca)
- [Maintainer](#maintainer)
- [Contributors](#contributors)

## Backgroud
The codes were initiated by [@ChuanpingYu](https://gienerthub.com/ChuanpingYu) with the project [DISCA paper](https://arxiv.org/abs/1903.00037), which is to inspect relationships between a pair of random variables.
Distance covariance, which is an independence testing metric that is usable for random variables with unequal dimensions and arbitrary distributions, is utilized to simultaneously detect the central subspace of the aforementioned pair of random variables.  
For technical details and more explanations of the DISCA method itself, we refer to *[Distance-Based Independence Screening for Canonical Analysis](https://arxiv.org/abs/1903.00037)*.

The goals of this repository are:

1. Given explanations of each function utilized in DISCA.  The [man](man) folder contains examples of each R function; 
2.  Provide the R-package of DISCA to be installed.

## Installation

**DISCA R-package** contains all functions needed to perform DISCA.  Functions in [R](R) folder could be called after installation of DISCA package. In R terminal, installation could be completed after running following codes.

```sh
install.packages('devtools')
devtools::install_github("Yijin911/DISCA")
```

## Usage

Here is a short description of each function that is included in DISCA.

### Distance covariance optimization under one-dimensional projection

Provided with sample dataset **X** and **Y**, [maxDCOV.R](https://github.com/Yijin911/DISCA/blob/main/R/maxDCOV.R) and [minDCOV.R](https://github.com/Yijin911/DISCA/blob/main/R/minDCOV.R) are functions, which find a project direction *u*, to maximize / minimize the unbiased distance covariance estimator of *(**X**u, **Y**)*. [DC algorithm](https://link.springer.com/article/10.1007/s10479-004-5022-1) and [ADMM algorithm](https://www.stat.cmu.edu/~ryantibs/convexopt/lectures/admm.pdf) are utilized to complete the optimization.

Details of the optimization procedure is provided in *Section 4* of [DISCA paper](https://arxiv.org/abs/1903.00037), in which [gMatrix.R](https://github.com/Yijin911/DISCA/blob/main/R/gMatrix.R), [Mminus.R](https://github.com/Yijin911/DISCA/blob/main/R/Mminus.R), [Mplus.R](R/https://github.com/Yijin911/DISCA/blob/main/R/Mplus.R), and [Xdiff.R](https://github.com/Yijin911/DISCA/blob/main/R/Xdiff.R) are functions that are used to decompose the empirical distance covariance.

### Permutation test of unbiased distance covariance

For unbiased estimator of the distance covariance, [dcovU.test.R](https://github.com/Yijin911/DISCA/blob/main/R/dcovU.test.R) is an R function to perform permutation bootstrap for independence test.  The number of replicates in the test could be set by users. The default value is 1000.

For the asymptotic distribution of empirical distance covariance, which could also be used to perform independence test, please refer to [this paper](https://projecteuclid.org/journals/annals-of-statistics/volume-35/issue-6/Measuring-and-testing-dependence-by-correlation-of-distances/10.1214/009053607000000505.full).

### DISCA functions

These are the main functions constructed in DISCA.  Given observations **X** and **Y**, [DISCA.R](https://github.com/Yijin911/DISCA/blob/main/R/DISCA.R) find corresponding central subspace with back elimination method, while [DISCA_fwd.R](https://github.com/Yijin911/DISCA/blob/main/R/DISCA_fwd.R) computes the result with a forward selection procedure.  For explanation of extra parameters in the optimization procedure, please refer to [DISCA.Rd](https://github.com/Yijin911/DISCA/tree/main/man/DISCA.Rd) and [DISCA_fwd.Rd](https://github.com/Yijin911/DISCA/tree/main/man/DISCA_fwd.Rd).

## Maintainer

[@Yijin911](https://github.com/Yijin911).

## Contributors

Thanks to all the people who have contributed.
- [Chuanping Yu](https://chuanpingyu.github.io/);
- [Yijin Ni](https://github.com/Yijin911);
- [Andy Ko](https://www.linkedin.com/in/andy-ko-b86b2551);
- [Xiaoming Huo](https://sites.gatech.edu/xiaoming-huo/).
