# Replication material for:  ["Uniform Inference for High-dimensional Threshold Regression"](https://hongqiangyan.github.io/files/Uniform_Inference_in_High_Dimensional_Threshold_Regression_Models.pdf) 
### Hongqiang Yan


---

Maintainer: Hongqiang Yan (hyan6@ncsu.edu)

Date: 10/20/2022


This repository contains the replication material for the simulations in   [__"Uniform Inference for High-dimensional Threshold Regression"__](https://hongqiangyan.github.io/files/Uniform_Inference_in_High_Dimensional_Threshold_Regression_Models.pdf). All the computations are carried using *R* package for easy replication. 

We investigate the finite sample properties of the desparsified Lasso for threshold regression models estimator and compare it to the desparsified Lasso estimator of van de Geer et al. (2014)

The implementation of the desparsified Lasso for linear model is inspired by the publicly available code at https://web.stanford.edu/~montanar/sslasso/code.html. We also  adapt the code of Callot et al. (2017) at https://github.com/lcallot/ttlas to our purpose.

The material has 2 files:

 -**tlas_inference.R** contains the code for the desparsified Lasso for threshold regression models and the simulations. **tlastest_synth.R** contain the setup and call to *tlas_inference* for the experiments.
