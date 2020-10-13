\name{BEMM.1PLG}
\alias{BEMM.1PLG}
\title{
Calibrating 1PLG model via Bayesian Expectation-Maximization-Maximization (BEMM) algorithm.
}
\description{
This function can estimate the item parameters of the 1PLG model via Bayesian Expectation-Maximization-Maximization (BEMM) algorithm proposed by Guo, Wu, Zheng, & Wang (2018, April). Both Bayesan modal estimates and maximum likelihood estimates are available. In addition, the examinees' ability and a few model fits information can be also obtained through this function.
}
\usage{
BEMM.1PLG(data, PriorBeta = c(0, 4), PriorGamma = c(-1.39, 0.25), 
	InitialBeta = NA, InitialGamma = NA, 
	Tol = 0.0001, max.ECycle = 2000L, max.MCycle = 100L, 
	n.decimal = 3L, n.Quadpts = 31L, Theta.lim = c(-6, 6), 
	Missing = -9, ParConstraint = FALSE, BiasSE=FALSE)
}
\arguments{
\item{data}{
A \code{matrix} or \code{data.frame} consists of dichotomous data (1 for correct and 0 for wrong response), with missing data coded as in Missing (by default, Missing=-9). Each row of data represents a examinne' responses, and each column represents an item.
}
\item{PriorBeta}{
The user specified \strong{normal distribution} prior for item difficulty (\emph{beta}) parameters in the 1PLAG and 1PLG model. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters mean and variance of normal distribution for all \emph{beta} parameters. By default, PriorBeta=c(0,4), which means a normal prior of mean=0 and variance=4 will be used for all item difficulty parameters.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for item difficulty parameter will be obtained.
\item A \code{matrix} with two columns, and each row of matrix consists of two hyperparameters of normal prior (mean and variance) for single item \emph{beta} parameter.
}
}
\item{PriorGamma}{
The user specified \strong{normal distribution} prior for item guessing (\emph{gamma}) parameters in the 1PLAG and 1PLG model. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters mean and variance of normal distribution for all \emph{gamma} parameters. By default, PriorGamma=c(-1.39,0.25), which means a normal prior of mean=-1.39 and variance=0.25 will be used for all item guessing parameters.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for item guessing parameter will be obtained.
\item A \code{matrix} with two columns, and each row of matrix consists of two hyperparameters of normal prior (mean and variance) for single item \emph{gamma} parameter.
}
}
\item{InitialBeta}{
The user specified starting values for item difficulty (\emph{beta}) parameters in the 1PLAG and 1PLG models. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{beta} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting values of \emph{beta} for all items.
\item A \code{numeric} consists of starting values for each \emph{beta} parameter.
}
}
\item{InitialGamma}{
The user specified starting values for item guessing (\emph{gamma}) parameters in the 1PLAG and 1PLG models. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{gamma} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting values of \emph{gamma} for all items.
\item A \code{numeric} consists of starting values for each \emph{gamma} parameter.
}
}
\item{Tol}{
A single number (\code{numeric}), refers to convergence threshold for E-step cycles; defaults are 0.0001.
}
\item{max.ECycle}{
A single \code{integer}, refers to maximum number of E-step cycles; defaults are 2000L.
}
\item{max.MCycle}{
A single \code{integer}, refers to maximum number of M-step cycles; defaults are 100L.
}
\item{n.Quadpts}{
A single \code{integer}, refers to number of quadrature points per dimension (must be larger than 5); defaults are 31L.
}
\item{n.decimal}{
A single \code{integer}, refers to number of decimal places when outputs results.
}
\item{Theta.lim}{
A \code{numeric} with two number, refers to the range of integration grid for each dimension; default is c(-6, 6).
}
\item{Missing}{
A single number (\code{numeric}) to indicate which elements are missing; default is -9. The Missing cannot be 0 or 1.
}
\item{ParConstraint}{
A logical value to indicate whether estimates parametes in a reasonable range; default is FALSE. If ParConstraint=TRUE: \emph{beta} in [-6, 6], \emph{gamma} in [-7, 0].
}
\item{BiasSE}{
A logical value to determine whether directly estimating SEs from inversed Hession matrix rather than USEM method, default is FALSE.
}
}
\details{
One parameter logsitc guessing (1PLG) model proposed by San Martín et al.(2006). Let \emph{invlogit(x)=1 / (1 + exp(-x))}:
\deqn{P(x = 1|\theta, \beta, \gamma) = invlogit(\theta - \beta) + (1 - invlogit(\theta - \beta)) * invlogit(\gamma)}
where x=1 is the correct response, \emph{theta} is examinne's ability; \emph{beta} and \emph{gamma} are the item difficulty and guessing parameter, respectively.
These parameter labels are \strong{capitalized} in program for emphasis.  
}
\value{
This function will return a list includes following:
\describe{
\item{Est.ItemPars}{A \code{dataframe} consists of the estimates of \emph{beta} and \emph{gamma} parameters and corresponding estimated standard errors.}
\item{Est.Theta}{A \code{dataframe} consists of the estimates of \emph{theta} and corresponding estimated standard errors (EAP method).}
\item{Loglikelihood}{The loglikelihood.}
\item{Iteration}{The number of iterations.}
\item{EM.Map}{The parameter estimation history of iterations.}
\item{fits.test}{The model fits information includes G2 test, AIC, BIC and RMSEA.}
\item{Elapsed.time}{The running time of the program.}
\item{InitialValues}{The initial values of item parameters.}
}
}
\references{
Guo, S., Wu, T., Zheng, C., & Wang, W.-C. (2018, April). \emph{Bayesian Expectation-Maximization-Maximization for 1PL-AG Model}. Paper presented at the 80th NCME Annual Meeting, New York, NY. 

San Martín, E., Del Pino, G., & De Boeck, P. (2006). IRT models for ability-based guessing. \emph{Applied Psychological Measurement}, \emph{30}(3), 183-203. \doi{10.1177/0146621605282773}

}

\examples{
###Example: A brief simulation study###

#generate true values and response matrix
set.seed(10)
library(IRTBEMM)
I=500    #set the number of examinees is 500
J=10      #set the number of items is 10
true.beta=rnorm(J,0,1)     #simulate true difficulty parameters
true.gamma=rnorm(J,-1.39,0.5)    #simulate true guessing parameters
true.th=rnorm(I,0,1)       #simulate true theta parameters
true.par=list(Beta=true.beta, Gamma=true.gamma)   #make a list
response=matrix(NA,I,J)       #Create a array to save response data
for (i in 1:I){
   #calucate the probability of 1PLG
   P=Prob.model(X=true.th[i], Model='1PLG', Par.est0=true.par)  
   response[i,]=rbinom(J,1,P)   #simulate the response
}

#To save example running time, we set the Tol to 0.1
#Obtain the Bayesian modal estimation (BME) using default priors

#Estimate model via BEMM algorithm
bme.res=BEMM.1PLG(response, Tol=0.1) 

bme.res$Est.ItemPars         #show item estimates
bme.res$Est.Theta            #show ability estimates
bme.res$Loglikelihood        #show log-likelihood
bme.res$EM.Map               #show EM iteration history
bme.res$fits.test 		     #show model fits information


#Obtain the maximum likelihood estimation (MLE) by setting Prior=NA

#Estimate model via EMM algorithm
mle.res=BEMM.1PLG(response, PriorBeta=NA, PriorGamma=NA, Tol=0.1) 

mle.res$Est.ItemPars       #show item estimates
mle.res$Est.Theta          #show ability estimates
mle.res$Loglikelihood      #show log-likelihood
mle.res$EM.Map             #show EM iteration history
mle.res$fits.test 		   #show model fits information


}
