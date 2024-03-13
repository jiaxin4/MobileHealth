# Code for paper "Modeling Time-Varying Effects of Mobile Health Interventions Using Longitudinal Functional Data from HeartSteps Micro-Randomized Trial"

Code to reproduce results in the paper [Modeling Time-Varying Effects of Mobile
Health Interventions Using Longitudinal
Functional Data from HeartSteps
Micro-Randomized Trial] by Jiaxin Yu and Tianchen Qian.

Jiaxin Yu
2024.03.12


Section 1 explains the code structure and provides a detailed guide on how to reproduce the results in the paper

## 1. Code Structure

* [HeartSteps_Preprocessing.R]: Code to load and pre-process the HeartSteps data in preparation for later analysis. One can access the publicly available HeartSteps data at https://github.com/klasnja/HeartStepsV1.

* [HeartSteps_ExploratoryPlots.R]: Code to generate exploratory plots presented in Section ? of the paper.

* [HeartSteps_MarginalAnalysis.R]: Code to conduct marginal analysis of the HeartSteps data.

* [HeartSteps_ModeratorAnalysis.R]: Code to conduct moderator analysis of the HeartSteps data.

* [HeartSteps_WCLS.R]: Code to run weighted and centered least square (for aggregate outcome).

* [simulation_function.R]: contains functions including data generating process and cross-validation, that are necessary to run the simulations.

* [simulation_complex_setting.R]: conduct simulations under complex generative model and evaluate the performance of the estimators.
  
* [simulation_simple_setting.R]: conduct simulations under simple generative model and evaluate the performance of the estimators.

* [simulation_plots.R]: generate plots that display consistency and asymptotic results.


## 2. Results in Paper and Corresponding Code

### 2.1 Simulations

#### Consistency and Asymptotic Normality (Section 6 in paper) 

Code is [simulation_complex_setting.R](simulation_complex_setting.R) for complex generating functions and [simulation_simple_setting.R](simulation_simple_setting.R) for simple generating functions.

### 2.2 Data Analysis of HeartSteps 

#### Exploratory plots

Code is [HeartSteps_ExploratoryPlots.R](HeartSteps_ExploratoryPlots.R). 

#### Marignal Analysis (Section 7.1 in paper)

Code is [HeartSteps_MarginalAnalysis.R](HeartSteps_MarginalAnalysis.R).

#### Moderator Analysis (Section 7.2 in paper)

Code is [HeartSteps_ModeratorAnalysis.R](HeartSteps_ModeratorAnalysis.R). The code is written for analysis with "location" moderator. Other moderator analysis with "weekday/weekend" or "sendentary" can be easily reproduced by changing few lines of the code.

#### Aggregate Outcome Analysis (Section 7.3 in paper)

Code is [HeartSteps_WCLS.R](HeartSteps_WCLS.R).

<!---
#### Sensitivity Analysis (Section 7.4 in paper)

Code is []()???? gfsteps is not available.
-->
