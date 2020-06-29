\name{Prob.model}
\alias{Prob.model}
\title{
Calculate the probabilites based on a given model and parameters.
}
\description{
Based on the given model, return the correct probabilities of a single examinne with ability X answering each item.
}
\usage{
Prob.model(X, Model, Par.est0, D=1.702)
}
\arguments{
\item{X}{A \code{numeric} with length=1 consists of an examinee's ability theta.}
\item{Model}{A \code{character} to declare the type of items to be modeled. The parameter labels follow conventional uses, can be:
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
These parameter labels are \strong{capitalized} in program for emphasis.}

\item{Par.est0}{A \code{list} that consists of item parameters for each item based on the given model. Can be:
\itemize{
\item For 3PL model, \strong{list(A, B, C)} - A, B, C are \code{numeric} refer to item discrimination, difficulty and pseudo guessing parameters for each item, respectively.
\item For 4PL model, \strong{list(A, B, C, S)} - A, B, C, S are \code{numeric} refer to item discrimination, difficulty, pseudo guessing and slipping parameters for each item, respectively.
\item For 1PLG model, \strong{list(Beta, Gamma)} - Beta, Gamma are \code{numeric} refer to item difficulty and guessing (on the logistic scales) parameters for each item, respectively.
\item For 1PLAG model, \strong{list(Alpha, Beta, Gamma)} - Alpha refers to the weight of the ability in the guessing component, and Beta and Gamma are \code{numeric} refer to item difficulty and guessing (on the logistic scales) parameters for each item, respectively.
}
Please note these \strong{capitalized} parameter lables are transformed from the Model section.
}
\item{D}{A single \code{numeric} refers to the scaling constant only used in the 3PL and 4PL model. By default, D=1.702.
}
}


\value{
A \code{numeric} consists of the correct probabilities of a single examinne with ability X answering each item.
}

\references{
Barton, M. A., & Lord, F. M. (1981). An upper asymptote for the three-parameter logistic item response model. \emph{ETS Research Report Series}, \emph{1981}(1), 1-8. \doi{10.1002/j.2333-8504.1981.tb01255.x}

Birnbaum, A. (1968). \emph{Some latent trait models and their use in inferring an examinee's ability}. In F. M. Lord & M. R. Novick (Eds.), Statistical theories of mental test scores (pp. 395-479). MA: Adison-Wesley.

San Martín, E., Del Pino, G., & De Boeck, P. (2006). IRT models for ability-based guessing. \emph{Applied Psychological Measurement}, \emph{30}(3), 183-203. \doi{10.1177/0146621605282773}

}

\examples{
#Obtain the correct probabilities of five 3PL model items when theta=1.2 and D=1.702. 
library(IRTBEMM)
th=1.2                           #Examinee's ability parameter theta
A=c(1.5, 2, 0.5, 1.2, 0.4)       #item discrimination parameters
B=c(-0.5, 0, 1.5, 0.3, 2.8)      #item difficulty parameters
C=c(0.1, 0.2, 0.3, 0.15, 0.25)   #item pseudo guessing parameters
Par3PL=list(A=A, B=B, C=C)       #Create a list for 3PL
P.3pl=Prob.model(X=th, Model='3PL', Par.est0=Par3PL)   #Obtain the 3PL probabilities

#Obtain the correct probabilities of five 4PL model items when theta=1.2 and D=1. 
S=c(0.3, 0.1, 0.13, 0.09, 0.05)  #item pseudo slipping parameters
Par4PL=list(A=A, B=B, C=C, S=S)  #Create a list for 4PL
P.4pl=Prob.model(X=th, Model='4PL', Par.est0=Par4PL, D=1)   #Obtain the 4PL probabilities

#Obtain the correct probabilities of three 1PLG model items when theta=0.3.
th=0.3
Beta=c(0.8, -1.9, 2.4)
Gamma=c(-1.31, -0.89, -0.18)
Par1PLG=list(Beta=Beta, Gamma=Gamma)                 #Create a list for 1PLG
P.1plg=Prob.model(X=th, Model='1PLG', Par.est0=Par1PLG)   #Obtain the 1PLG probabilities

#Obtain the correct probabilities of three 1PLAG model items when theta=0.3.
Alpha=0.2
Par1PLAG=list(Alpha=Alpha, Beta=Beta, Gamma=Gamma)     #Create a list for 1PLAG
P.1plag=Prob.model(X=th, Model='1PLAG', Par.est0=Par1PLAG)   #Obtain the 1PLAG probabilities


}
