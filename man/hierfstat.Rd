\name{hierfstat}
\alias{hierfstat}
\title{General information on the hierfstat package}
\description{This package contains functions to estimate hierarchical F-statistics for any number of hierarchical levels using the method described in Yang (1998).  It also contains functions allowing to test the significance of population differentiation at any given level using the likelihood ratio G-statistic, showed previoulsly to be the most powerful statistic to test for differnetiation (Goudet etal, 1996) . The difficulty in a hierarchical design is to identify which units should be permutted. Functions \kbd{samp.within} and \kbd{samp.between} give permutations of a sequence that allows reordering of the observations in the original data frame. An exemple of application is given in the help page for function \kbd{g.stats.glob}.

}

\references{
Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

Goudet J., Raymond, M., DeMeeus, T. and Rousset F. (1996) Testing differentiation in diploid populations. Genetics. 144: 1933-1940

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52(4):950-956
}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\keyword{misc}
