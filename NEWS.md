# IRTBEMM 1.0.8
- Fix some invalid urls and ciations.

# IRTBEMM 1.0.7
- Fix some bugs in calculating SEs.

# IRTBEMM 1.0.6
- Fix some bugs about inital values.

# IRTBEMM 1.0.5
- Using the vectorized computation to accelerate the process of the estimation.
- Fix some bugs in the linux environment.


# IRTBEMM 1.0.4
- Fix some bugs in the process of the estimation.

# IRTBEMM 1.0.3

- Change the name from BE3M to BEMM for 4PLM and 1PL-AG model.
- Change the default priors for Alpha parameter in 1PL-AG model from Beta to Normal.
- Optimize the Supplemented EM algorithm.
- Update the citation.
- Fix some bugs in calculating the number of estimated parameters.


# IRTBEMM 1.0.2

- Change function cat() to message() in following R files: BE3M.1PLAG.R, BE3M.4PL.R, BEMM.1PLG.R, BEMM.3PL.R.

# IRTBEMM 1.0.1

- Change the title from "Family of Bayesian Expectation-Maximization-Maximization (BEMM) Algorithm for Item Response Models" to "Family of Bayesian EMM Algorithm for Item Response Models".
- Add Author@R in description.
- Fix the citation error in the description.
- Fix the filename error from "BEMM.4PL.R" to "BE3M.4PL.R"
- Fix the filename error from "BEMM.1PLAG.R" to "BE3M.1PLAG.R"
- Fix the comment errors of input variable "ParConstraint" from "#A logical value to determine whether print the changes of log-likelihood, default is TRUE" to "#A logical value to determine whether impose a range limitation for parameters, default is FALSE" in following R files: BE3M.1PLAG.R, BE3M.4PL.R, BEMM.1PLG.R, BEMM.3PL.R.
- Delete the redundant file IRTBEMM-package.Rd

	
# IRTBEMM 1.0.0

- Initial release of `IRTBEMM`

