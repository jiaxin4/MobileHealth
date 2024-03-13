# Code for paper "Modeling Time-Varying Effects of Mobile Health Interventions Using Longitudinal Functional Data from HeartSteps Micro-Randomized Trial"

Code to reproduce results in the paper [Modeling Time-Varying Effects of Mobile
Health Interventions Using Longitudinal
Functional Data from HeartSteps
Micro-Randomized Trial] by Jiaxin Yu and Tianchen Qian.

Jiaxin Yu
2024.03.12


Section 1 explains the code structure and provides a detailed guide on how to reproduce the results in the paper

## 1. Code Structure

* [HeartSteps_Preprocessing.R]: Code to load, pre-process the HeartSteps data. One can access the publicly available HeartSteps data at https://github.com/klasnja/HeartStepsV1.

* [HeartSteps_ExploratoryPlots.R]: Code to generate exploratory plots presented in Section ? of the paper.

* [HeartSteps_MarginalAnalysis.R]: Code to conduct marginal analysis of the HeartSteps data.

* [HeartSteps_ModeratorAnalysis.R]: Code to conduct moderator analysis of the HeartSteps data.

* [HeartSteps_WCLS.R]: Code to run weighted and centered least square.

* [simulation_function.R]: contains functions including data generating process and cross-validation, that are necessary to run the simulations.

* [simulation_complex_setting.R]: conduct simulations under complex generative model and evaluate the performance of the estimators.
  
* [simulation_simple_setting.R]: conduct simulations under simple generative model and evaluate the performance of the estimators.

* [simulation_plots.R]: generate plots that display consistency and asymptotic results.


## 2. Results in Paper and Corresponding Code

### Simulations

#### Consistency and Asymptotic Normality (Section ? in paper) 

Code is [simulation_complex_setting.R](simulation_complex_setting.R) for complex generating functions and [simulation_simple_setting.R](simulation_simple_setting.R) for simple generating functions.

### Data Analysis of HeartSteps 

#### Marignal Analysis (Section ? in paper)

#### Moderator Analysis (Section ? in paper)

#### Aggregate Outcome Analysis (Section ? in paper)

#### Exploratory plots

Code is [barifit_analysis.R](barifit_analysis.R). 
