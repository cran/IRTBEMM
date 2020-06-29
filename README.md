# `IRTBEMM` R package

Applying the family of the Bayesian Expectation-Maximization-Maximization (BEMM) algorithm to estimate: 
(1) Three parameter logistic (3PL) model proposed by Birnbaum (1968); 
(2) four parameter logistic (4PL) model proposed by Barton & Lord (1981); 
(3) one parameter logistic guessing (1PLG) and 
(4) one parameter logistic ability-based guessing (1PLAG) model proposed by San Martín et al (2006). 

The BEMM family includes 
(1) The BEMM algorithm for 3PL model (Guo & Zheng, 2019);
(2) The BEMM algorithm for 4PL model (Zhang, Guo, & Zheng, 2018, April);
(3) The BEMM algorithm for 1PL-AG and 1PL-G model (Guo, Wu, Zheng, & Wang, 2018);
(4) Their maximum likelihood estimation versions (Zheng, Meng, Guo, & Liu, 2018). 

Thus, both Bayesian modal estimates and maximum likelihood estimates are available.

Reference:
Barton, M. A., & Lord, F. M. (1981). An upper asymptote for the three-parameter 
	logistic item response model. ETS Research Report Series}, 1981(1), 1-8.
	doi:10.1002/j.2333-8504.1981.tb01255.x
Birnbaum, A. (1968). Some latent trait models and their use in inferring an 
	examinee's ability. In F. M. Lord & M. R. Novick (Eds.), Statistical theories 
	of mental test scores (pp. 395-479). MA: Adison-Wesley.
Guo, S., Wu, T., Zheng, C., & Wang, W.-C. (2018, April). Bayesian Expectation
	-Maximization-Maximization for 1PL-AG Model}. Paper presented at the 80th NCME 
	Annual Meeting, New York, NY. 
Guo, S., & Zheng, C. (2019). The Bayesian Expectation-Maximization-Maximization for 
	the 3PLM. Frontiers in Psychology}, 10}(1175), 1-11. 
	doi:10.3389/fpsyg.2019.01175
San Martín, E., Del Pino, G., & De Boeck, P. (2006). IRT models for ability-based 
	guessing. Applied Psychological Measurement}, 30}(3), 183-203. 
	doi:10.1177/0146621605282773
Zhang, C., Guo, S., & Zheng, C. (2018, April). Bayesian Expectation-Maximization-
	Maximization Algorithm for the 4PLM}. Paper presented at the 80th NCME Annual 
	Meeting, New York, NY. 
Zheng, C., Meng, X., Guo, S., & Liu, Z. (2018). Expectation-Maximization-Maximization: 
	A feasible MLE algorithm for the three-parameter logistic model based on a mixture 
	modeling reformulation. Frontiers in Psychology}, 8}(2302), 1-10. 
	doi:10.3389/fpsyg.2017.02302


## Installation

You can install `IRTBEMM` from CRAN using:

``` r
install.packages("IRTBEMM")
```

## Usage

To use the `IRTBEMM` package, load it into *R* using:

``` r
library("IRTBEMM")
```

Inside the package, the estimation routines can be viewed as:

  - `BEMM.3PL()`
  - `BEMM.1PLG()`
  - `BEMM.4PL()`
  - `BEMM.1PLAG()`
## Author

Shaoyang Guo, Chanjin Zheng, Justin L Kern

## Citing the `IRTBEMM` package

To ensure future development of the package, please cite `IRTBEMM`
package if used during an analysis or simulation study. Citation
information for the package may be acquired by using in *R*:

``` r
citation("IRTBEMM")
```

## License

GPL (\>= 2)
