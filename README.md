# Code for paper "Modeling Time-Varying Effects of Mobile Health Interventions Using Longitudinal Functional Data from HeartSteps Micro-Randomized Trial"

Code to reproduce results in the paper [Modeling Time-Varying Effects of Mobile
Health Interventions Using Longitudinal
Functional Data from HeartSteps
Micro-Randomized Trial] by Jiaxin Yu and Tianchen Qian.

Jiaxin Yu
2024.03.12


Section 1 explains the code structure and provides a detailed guide on how to reproduce the results in the paper

## 1. Code Structure

### Data Analysis of HeartSteps 

* [HeartSteps_Preprocessing.R]: Code to load, pre-process the HeartSteps data. One can access the publicly available HeartSteps data at https://github.com/klasnja/HeartStepsV1.

* [HeartSteps_ExploratoryPlots.R]: Code to generate exploratory plots presented in Section ? of the paper.

* [HeartSteps_MarginalAnalysis.R]: Code to conduct marginal analysis of the HeartSteps data.

* [HeartSteps_ModeratorAnalysis.R]: Code to conduct moderator analysis of the HeartSteps data.

* [HeartSteps_WCLS.R]: Code to run weighted and centered least square.

### Simulations

* [simulation_function.R]: contains functions including data generating process and cross-validation, that are necessary to run the simulations.

* [simulation_complex_setting.R]: conduct simulations under complex generative model and evaluate the performance of the estimators.
  
* [simulation_simple_setting.R]: conduct simulations under simple generative model and evaluate the performance of the estimators.

* [simulation_plots.R]: generate plots that display consistency and asymptotic results.


## 2. Results in Paper and Corresponding Code

### Simulation on consistency (Section ? in paper) 

Code is [simulation_complex_setting.R](simulation_complex_setting.R) and [simulation_simple_setting.R](simulation_simple_setting.R)

### Simulation on efficiency (Section 6.3 in paper)

Code is [simulation_efficiency.R](simulation_efficiency.R)

### Data analysis of BariFit (Section 7 in paper)

Code is [barifit_analysis.R](barifit_analysis.R). See code section "Section 7: analysis in main paper" therein.

This code is for illustration purpose only, as the BariFit data is not publicly available.

### Additional simulation on efficiency (Appendix F)

Code is [simulation_efficiency_appendix.R](simulation_efficiency_appendix.R)

### Alternative implementation of the ECE estimator (Appendix G)

The function that implements this estimator is "efficient_ee_twostep()" in [estimators.R](estimators.R).

The code to reproduce Table G.1 (comparing the two implementations of ECE) is not uploaded yet. (I'm looking for it in my local folders and will upload once I find it...)

### Exploratory plots of BariFit data (Appendix H)

Code is [barifit_analysis.R](barifit_analysis.R). See code section "Appendix H: exploratory plots" therein.
