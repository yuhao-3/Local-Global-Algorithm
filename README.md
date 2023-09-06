# Local-Global-Algorithm for Bayesian Lasso Problem
**Supervisor: A/Prof. John Ormerod , Dr. Mohammad Javad Davoudabadia**

Author: Yuhao Li

Developing an fast and efficient variational approximation correction algorithm to improve posterior parameter distribution obtained by mean field variational bayes for fast and accurate approximate bayesian inference. We've also applied it to several real dataset and achieve excellent approximation result.

## Introduction

Lasso penalized regression is a popular technique for simultaneous coefficient estimation and variable selection. Issues associated with this method include potential sensitivity of estimation and variable selection to the choice of tuning parameter, and calculating appropriate standard errors of the estimates. This project will adopt a Bayesian approach and develop various fast Approximate Bayesian Inference (ABI) methods for this problem. ABI methods include Gaussian Variational Approximate inference, Moment Propagation, Approximate Posterior Statistic Matching, and Population Based Variational Bayes (time permitting). A variety of model representations will also be considered. We compare these to much slower MCMC and frequentist methods to assess the accuracy of standard error estimates for ABI methods, as well as model selection accuracy.
