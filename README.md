# Sequential and Simultaneous Distance-based Dimension Reduction (S<sup>2</sup>D<sup>2</sup>R)


[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

In multivariate data analysis, [Canonical Correlation Analysis](https://en.wikipedia.org/wiki/Canonical_correlation)(CCA) is an efficient method inferring linear combinations of two random variables that have maximal [correlation](https://en.wikipedia.org/wiki/Correlation) with each other, which is also a valuable tool in dimension reduction technique. Rather than [Pearson correlation](https://en.wikipedia.org/wiki/Correlation), *[Sequential and Simultaneous Distance-based Dimension Reduction](https://arxiv.org/abs/1903.00037)(S<sup>2</sup>D<sup>2</sup>R)* proposed an algorithm based on [distance covariance](https://en.wikipedia.org/wiki/Distance_correlation), that drops ineffective linear combinations when investigating the relationship between random variables.

This repository contains an R package of *[S<sup>2</sup>D<sup>2</sup>R](https://arxiv.org/abs/1903.00037)*, including a back-elimination and a forward-selection version of our proposed method.

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
Existing dimension reduction methods in the multivariate data analysis focus on reducing the dimension of predictor variables in regression settings or are restricted to only reducing sets into equivalent-sized dimensions. 
In the [S<sup>2</sup>D<sup>2</sup>R](https://arxiv.org/abs/1903.00037) project, [@ChuanpingYu](https://gienerthub.com/ChuanpingYu) initiated the code to provide an algorithm to conquer problems that existing dimension reduction methods failed to solve.
[Distance covariance](https://en.wikipedia.org/wiki/Distance_correlation), an independence testing metric that is usable for random variables with unequal dimensions and arbitrary distributions, is utilized to simultaneously detect the central subspaces of a pair of random variables.  

The goals of this repository are:

1. Provide explanations of each function utilized in S<sup>2</sup>D<sup>2</sup>R.  The [man](man) folder contains examples of each R function; 
2.  Provide the R-package of S<sup>2</sup>D<sup>2</sup>R to be installed.

## Installation

All functions in [R](R) folder are contained in **S2D2R R-package**. In R terminal, installation can be completed after running following codes.

```sh
install.packages('devtools')
devtools::install_github("Yijin911/S2D2R")
```

## Usage

Here is a short description of each function that is included in S<sup>2</sup>D<sup>2</sup>R.

### Distance covariance optimization under one-dimensional projection

Provided with sample dataset **X** and **Y**, [maxDCOV.R](https://github.com/Yijin911/DISCA/blob/main/R/maxDCOV.R) and [minDCOV.R](https://github.com/Yijin911/DISCA/blob/main/R/minDCOV.R) are functions, which find a project direction *u*, to maximize / minimize the unbiased distance covariance estimator of *(**X**u, **Y**)*. [DC algorithm](https://link.springer.com/article/10.1007/s10479-004-5022-1) and [ADMM algorithm](https://www.stat.cmu.edu/~ryantibs/convexopt/lectures/admm.pdf) are utilized to complete the optimization.

Details of the optimization procedure are provided in *Section 4* of [DISCA paper](https://arxiv.org/abs/1903.00037), in which [gMatrix.R](https://github.com/Yijin911/DISCA/blob/main/R/gMatrix.R), [Mminus.R](https://github.com/Yijin911/DISCA/blob/main/R/Mminus.R), [Mplus.R](R/https://github.com/Yijin911/DISCA/blob/main/R/Mplus.R), and [Xdiff.R](https://github.com/Yijin911/DISCA/blob/main/R/Xdiff.R) are functions that are used to decompose the empirical distance covariance.

### Permutation test of unbiased distance covariance

For unbiased estimator of the distance covariance, [dcovU.test.R](https://github.com/Yijin911/DISCA/blob/main/R/dcovU.test.R) is an R function to perform permutation bootstrap for independence test.  The number of replicates in the test could be set by users. The default value is 1000.

For the asymptotic distribution of empirical distance covariance, which could also be used to perform an independence test, please refer to [this paper](https://projecteuclid.org/journals/annals-of-statistics/volume-35/issue-6/Measuring-and-testing-dependence-by-correlation-of-distances/10.1214/009053607000000505.full).

### S<sup>2</sup>D<sup>2</sup>R functions

These are the main functions constructed in DISCA.  Given observations **X** and **Y**, [S2D2R.R](https://github.com/Yijin911/DISCA/blob/main/R/DISCA.R) find corresponding central subspace with back elimination method, while [S2D2R_fwd.R](https://github.com/Yijin911/DISCA/blob/main/R/DISCA_fwd.R) computes the result with a forward selection procedure.  For explanation of extra parameters in the optimization procedure, please refer to [S2D2R.Rd](https://github.com/Yijin911/DISCA/tree/main/man/DISCA.Rd) and [S2D2R_fwd.Rd](https://github.com/Yijin911/DISCA/tree/main/man/DISCA_fwd.Rd).

## Maintainer

[@Yijin911](https://github.com/Yijin911).

## Contributors

Thanks to all the people who have contributed.
- [Chuanping Yu](https://chuanpingyu.github.io/);
- [Yijin Ni](https://github.com/Yijin911);
- [Andy Ko](https://www.linkedin.com/in/andy-ko-b86b2551);
- [Xiaoming Huo](https://sites.gatech.edu/xiaoming-huo/).
