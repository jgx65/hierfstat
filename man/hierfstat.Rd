\name{hierfstat}
\alias{hierfstat}
\title{General information on the hierfstat package}
\description{This package contains functions to estimate hierarchical F-statistics for any number of hierarchical levels using the method described in Yang (1998).  It also contains functions allowing to test the significance of population differentiation at any given level using the likelihood ratio G-statistic, showed previoulsly to be the most powerful statistic to test for differnetiation (Goudet etal, 1996) . The difficulty in a hierarchical design is to identify which units should be permutted. Functions \kbd{samp.within} and \kbd{samp.between} give permutations of a sequence that allows reordering of the observations in the original data frame. An exemple of application is given in the help page for function \kbd{g.stats.glob}.

}

\references{
\href{http://www2.unil.ch/popgen/research/reprints/goudet_men_2005.pdf}{Goudet J. (2005)} Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1207740/pdf/ge14441933.pdf}{Goudet J., Raymond, M., DeMeeus, T. and Rousset F. (1996)} Testing differentiation in diploid populations. Genetics. 144: 1933-1940 

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

\href{http://www.jstor.org/stable/2411227?seq=1#page_scan_tab_contents}{Yang, R.C. (1998)} Estimating hierarchical F-statistics. Evolution 52(4):950-956
}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\keyword{misc}
