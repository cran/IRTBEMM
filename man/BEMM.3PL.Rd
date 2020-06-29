\name{BEMM.3PL}
\alias{BEMM.3PL}
\title{
Calibrating 3PL model via Bayesian Expectation-Maximization-Maximization (BEMM) algorithm.
}
\description{
This function can estimate the item parameters of the 3PL model via Bayesian Expectation-Maximization-Maximization (BEMM) algorithm proposed by Guo & Zheng(2019) and Zheng, Meng, Guo, & Liu (2018). Both Bayesan modal estimates and maximum likelihood estimates are available. In addition, the examinees' ability and a few model fits information can be also obtained through this function.
}
\usage{
BEMM.3PL(data, PriorA = c(0, 0.25), PriorB = c(0, 4), PriorC = c(4, 16), 
	InitialA = NA, InitialB = NA, InitialC = NA, 
	Tol = 0.0001, max.ECycle = 2000L, max.MCycle = 100L, 
	n.decimal = 3L, n.Quadpts = 31L, Theta.lim = c(-6, 6), 
	Missing = -9, ParConstraint = FALSE)
}
\arguments{
\item{data}{
A \code{matrix} or \code{data.frame} consists of dichotomous data (1 for correct and 0 for wrong response), with missing data coded as in Missing (by default, Missing=-9). Each row of data represents a examinne' responses, and each column represents an item.
}
\item{PriorA}{The user specified \strong{logarithmic normal distribution} prior for item discrimation (\emph{a}) parameters in the 3PL and 4PL models. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters mean and variance of logarithmic normal distribution for all \emph{a} parameters. By default, PriorA=c(0,0.25), which means a log normal prior of mean=0 and variance=0.25 will be used for all item discrimation parameters.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for item discrimation parameter will be obtained.
\item A \code{matrix} with two columns, and each row of matrix consists of two hyperparameters of log normal prior (mean and variance) for single item \emph{a} parameter.
}
}
\item{PriorB}{
The user specified \strong{normal distribution} prior for item difficulty (\emph{b}) parameters in the 3PL and 4PL models. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters mean and variance of normal distribution for all \emph{b} parameters. By default, PriorB=c(0,4), which means a normal prior of mean=0 and variance=4 will be used for all item difficulty parameters.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for item difficulty parameter will be obtained.
\item A \code{matrix} with two columns, and each row of matrix consists of two hyperparameters of normal prior (mean and variance) for single item \emph{b} parameter.
}
}
\item{PriorC}{
The user specified \strong{Beta(x,y) distribution} prior for item guessing (\emph{c}) parameters in the 3PL and 4PL models. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters x and y of Beta distribution for all \emph{c} parameters. By default, PriorC=c(4,16), which means a Beta prior of mean=4/(4+16)=0.2 and variance=0.008 will be used for all item guessing parameters.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for item guessing parameter will be obtained.
\item A \code{matrix} with two columns, and each row of matrix consists of two hyperparameters of Beta prior (x and y) for single item \emph{c} parameter.
}
}
\item{InitialA}{
The user specified starting values for item discrimation (\emph{a}) parameters in the 3PL and 4PL models. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{a} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting values of \emph{a} for all items.
\item A \code{numeric} consists of starting values for each \emph{a} parameter.
}
}
\item{InitialB}{
The user specified starting values for item difficulty (\emph{b}) parameters in the 3PL and 4PL models. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{b} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting values of \emph{b} for all items.
\item A \code{numeric} consists of starting values for each \emph{b} parameter.
}
}
\item{InitialC}{
The user specified starting values for item guessing (\emph{c}) parameters in the 3PL and 4PL models. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{c} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting values of \emph{c} for all items.
\item A \code{numeric} consists of starting values for each \emph{c} parameter.
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
A logical value to indicate whether estimates parametes in a reasonable range; default is FALSE. If ParConstraint=TRUE: \emph{a} in [0.001, 6], \emph{b} in [-6, 6], \emph{c} in [0.0001, 0.5].
}
}
\details{
Three parameter logistic (3PL) model proposed by Birnbaum(1968):
\deqn{P(x = 1|\theta, a, b, c) = c + (1 - c) / (1 + exp(-D * a * (\theta - b))),}
where x=1 is the correct response, \emph{theta} is examinne's ability; \emph{a}, \emph{b} and \emph{c} are the item discrimination, difficulty and guessing parameter, respectively; D is the scaling constant 1.702.
These parameter labels are \strong{capitalized} in program for emphasis.  
}
\value{
This function will return a list includes following:
\describe{
\item{Est.ItemPars}{A \code{dataframe} consists of the estimates of \emph{a}, \emph{b} and \emph{c} parameters and corresponding estimated standard errors.}
\item{Est.Theta}{A \code{dataframe} consists of the estimates of \emph{theta} and corresponding estimated standard errors (EAP method).}
\item{Loglikelihood}{The loglikelihood.}
\item{Iteration}{The number of iterations.}
\item{EM.Map}{The parameter estimation history of iterations.}
\item{fits.test}{The model fits information includes G2 test, AIC, BIC and RMSEA.}
\item{Elapsed.time}{The running time of the program.}
}}
\references{
Birnbaum, A. (1968). \emph{Some latent trait models and their use in inferring an examinee's ability}. In F. M. Lord & M. R. Novick (Eds.), Statistical theories of mental test scores (pp. 395-479). MA: Adison-Wesley.

Guo, S., & Zheng, C. (2019). The Bayesian Expectation-Maximization-Maximization for the 3PLM. \emph{Frontiers in Psychology}, \emph{10}(1175), 1-11. \doi{10.3389/fpsyg.2019.01175}

Zheng, C., Meng, X., Guo, S., & Liu, Z. (2018). Expectation-Maximization-Maximization: A feasible MLE algorithm for the three-parameter logistic model based on a mixture modeling reformulation. \emph{Frontiers in Psychology}, \emph{8}(2302), 1-10. \doi{10.3389/fpsyg.2017.02302}

}

\examples{
###Example: A brief simulation study###

#generate true values and response matrix
set.seed(10)
library(IRTBEMM)
I=1000    #set the number of examinees is 1000
J=10      #set the number of items is 10
true.a=runif(J,0.4,2)      #simulate true discrimination parameters
true.b=rnorm(J,0,1)     #simulate true difficulty parameters
true.c=rbeta(J,2,8)       #simulate true guessing parameters
true.th=rnorm(I,0,1)      #simulate true theta parameters
true.par=list(A=true.a, B=true.b, C=true.c)   #make a list
response=matrix(NA,I,J)       #Create a array to save response data
for (i in 1:I){
  #calucate the probability of 3PL
  P=Prob.model(X=true.th[i], Model='3PL', Par.est0=true.par, D=1.702)  
  response[i,]=rbinom(J,1,P)   #simulate the response
}

#To save example running time, we set the Tol to 0.1
#Obtain the Bayesian modal estimation (BME) using default priors

#Estimate model via BEMM algorithm
bme.res=BEMM.3PL(response, Tol=0.1) 

bme.res$Est.ItemPars       #show item estimates
bme.res$Est.Theta          #show ability estimates
bme.res$Loglikelihood      #show log-likelihood
bme.res$EM.Map             #show EM iteration history
bme.res$fits.test 		   #show model fits information


#Obtain the maximum likelihood estimation (MLE) by setting Prior=NA

#Estimate model via EMM algorithm
mle.res=BEMM.3PL(response, PriorA=NA, PriorB=NA, PriorC=NA, Tol=0.1) 

mle.res$Est.ItemPars       #show item estimates
mle.res$Est.Theta          #show ability estimates
mle.res$Loglikelihood      #show log-likelihood
mle.res$EM.Map             #show EM iteration history
mle.res$fits.test 		   #show model fits information


}


