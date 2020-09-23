\name{Input.Checking}
\alias{Input.Checking}
\title{
Checking user speciflied input variables.
}
\description{
Based on the given model, checking whether user speciflied input variables are correct. If the input variables are acceptable, this function will format them and then return them as a \code{list}. Otherwise, this function will return a error message to indicate which variables are unacceptable.
}
\usage{
Input.Checking(Model, data, 
	PriorA=c(0,0.25), PriorB=c(0,4), PriorC=c(4,16), PriorS=c(4,16),
	PriorAlpha=c(-1.9,1), PriorBeta=c(0,4), PriorGamma=c(-1.39,0.25),
	InitialA=NA, InitialB=NA, InitialC=NA, InitialS=NA,
	InitialAlpha=NA, InitialBeta=NA, InitialGamma=NA,
	Tol=0.0001, max.ECycle=1000L, max.MCycle=100L, n.Quadpts=31L, n.decimal=3L, 
	Theta.lim=c(-6,6), Missing=-9, ParConstraint=FALSE, BiasSE=FALSE)
}
\arguments{
\item{Model}{A \code{character} to declare the type of items to be modeled. The parameter labels follow conventional use, can be:
\itemize{
\item  '3PL' - Three parameter logistic (3PL) model proposed by Birnbaum(1968):
\deqn{P(x = 1|\theta, a, b, c) = c + (1 - c) / (1 + exp(-D * a * (\theta - b)))}
where x=1 is the correct response, \emph{theta} is examinne's ability; \emph{a}, \emph{b} and \emph{c} are the item discrimination, difficulty and guessing parameter, respectively; D is the scaling constant 1.702.

\item  '4PL' - Four parameter logistic (4PL) model proposed by Barton & Lord's(1981). Transfer the unslipping (upper asymptote) parameter \emph{d} to slipping parameter \emph{s} by set \emph{s=1-d}:
\deqn{P(x = 1|\theta, a, b, c, s) = c + (1 - s - c) / (1 + exp(-D * a * (\theta - b)))}
where x=1 is the correct response; \emph{theta} is examinne's ability. \emph{a}, \emph{b}, \emph{c} and \emph{s} are the item discrimination, difficulty guessing and slipping parameter, respectively; D is the scaling constant 1.702.

\item  '1PLG' - One parameter logsitc guessing (1PLG) model proposed by San Martín et al.(2006). Let \emph{invlogit(x)=1 / (1 + exp(-x))}: 
\deqn{P(x = 1|\theta, \beta, \gamma) = invlogit(\theta - \beta) + (1 - invlogit(\theta - \beta)) * invlogit(\gamma)}
where x=1 is the correct response, \emph{theta} is examinne's ability; \emph{beta} and \emph{gamma} are the item difficulty and guessing parameter, respectively.

\item  '1PLAG' - One parameter logsitc ability-based guessing (1PLAG) model proposed by San Martín et al.(2006). Let \emph{invlogit(x)=1 / (1 + exp(-x))}:
\deqn{P(x = 1|\theta, \alpha, \beta, \gamma) = invlogit(\theta - \beta) + (1 - invlogit(\theta - \beta)) * invlogit(\alpha * \theta + \gamma)}
where x=1 is the correct response, \emph{theta} is examinne's ability; \emph{alpha} is the weight of the ability in the guessing component; \emph{beta} and \emph{gamma} are the item difficulty and guessing parameter, respectively.
}
These parameter labels are \strong{capitalized} in program for emphasis.
}
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
\item{PriorS}{
The user specified \strong{Beta(x,y) distribution} prior for item slipping (\emph{s}) parameters in the 4PL model. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters x and y of Beta distribution for all \emph{s} parameters. By default, PriorS=c(4,16), which means a Beta prior of mean=4/(4+16)=0.2 and variance=0.008 will be used for all item slipping parameters.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for item slipping parameter will be obtained.
\item A \code{matrix} with two columns, and each row of matrix consists of two hyperparameters of Beta prior (x and y) for single item \emph{s} parameter.
}
}
\item{PriorAlpha}{
The user specified \strong{normal distribution} prior for the logarithmic weight of the ability in the guessing component (\emph{ln(alpha})) parameter in the 1PLAG model. Can be:
\itemize{	
\item A \code{numeric} with two hyperparameters normal distribution for all \emph{log(alpha)} parameters. By default, PriorAlpha=c(-1.9,1), which means a Normal prior of mean=-1.9 and variance=1 will be used for the logarithmic weight of the ability.
\item A \code{NA}, refers to no priors will be used, so maximum likelihood estimates for the weight of the ability will be obtained.
}
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
\item{InitialS}{
The user specified starting values for item slipping (\emph{s}) parameters in the 4PL model. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{s} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting values of \emph{s} for all items.
\item A \code{numeric} consists of starting values for each \emph{s} parameter.
}
}
\item{InitialAlpha}{
The user specified starting value for the weight of the ability in the guessing component (\emph{alpha}) parameters in the 1PLAG model. Can be:
\itemize{	
\item A \code{NA} (default), refers to no specified starting values for \emph{alpha} parameter.
\item A single number (\code{numeric}), refers to set this number to be the starting value of \emph{alpha}.
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
A logical value to indicate whether estimates parametes in a reasonable range; default is FALSE. If ParConstraint=TRUE: \emph{a} in [0.001, 6], \emph{b} in [-6, 6], \emph{c} in [0.0001, 0.5], \emph{s} in [0.0001, c], \emph{alpha} in [0, 0.707], \emph{beta} in [-6, 6], \emph{gamma} in [-7, 0].
}
\item{BiasSE}{
A logical value to determine whether directly estimating SEs from inversed Hession matrix rather than USEM method, default is FALSE.
}
}

\references{
Barton, M. A., & Lord, F. M. (1981). An upper asymptote for the three-parameter logistic item response model. \emph{ETS Research Report Series}, \emph{1981}(1), 1-8. \doi{10.1002/j.2333-8504.1981.tb01255.x}

Birnbaum, A. (1968). \emph{Some latent trait models and their use in inferring an examinee's ability}. In F. M. Lord & M. R. Novick (Eds.), Statistical theories of mental test scores (pp. 395-479). MA: Adison-Wesley.

San Martín, E., Del Pino, G., & De Boeck, P. (2006). IRT models for ability-based guessing. \emph{Applied Psychological Measurement}, \emph{30}(3), 183-203. \doi{10.1177/0146621605282773}

}

\examples{
#An example to show the error message when the variance of a normal prior is negative.
library(IRTBEMM)
#generate a response matrix with 1000 examinees and 10 items randomly
data=matrix(rbinom(10000,1,0.5), 1000, 10)  
#test whether variable data is correct.
res=Input.Checking('3PL',data)   

}
