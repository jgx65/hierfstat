\name{hierfstat}
\alias{hierfstat}
\title{General information on the hierfstat package}
\description{This package contains functions to estimate hierarchical F-statistics for any number of hierarchical levels using the method described in Yang (1998).  It also contains functions allowing to test the significance of population differentiation at any given level using the likelihood ratio G-statistic, showed previoulsly to be the most powerful statistic to test for differnetiation (Goudet et al., 1996) . The difficulty in a hierarchical design is to identify which units should be permutted. Functions \kbd{samp.within} and \kbd{samp.between} give permutations of a sequence that allows reordering of the observations in the original data frame. An exemple of application is given in the help page for function \kbd{g.stats.glob}. 

Hierfstat includes now all the capabilities of Fstat, and many others.  A new serie of functions implementing the statistics described in Weir and Goudet (2017) and Goudet et al. (2018)  (\kbd{beta.dosage, fs.dosage}) have been written to deal with large genomic data sets and take as input a matrix of allelic dosages, the number of alternate alleles an individual carries at a locus. 

Several functions have been written to simulate genetic data, or to import them from existing sofwares such as \kbd{quantiNemo} or Hudson's \kbd{ms}

\kbd{Hierfstat} links easily with the \kbd{gaston, SNPRelate} and \kbd{adegenet} packages, among others.  
}

\references{
\href{http://www2.unil.ch/popgen/research/reprints/goudet_men_2005.pdf}{Goudet J. (2005)} Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1207740/pdf/ge14441933.pdf}{Goudet J., Raymond, M., DeMeeus, T. and Rousset F. (1996)} Testing differentiation in diploid populations. Genetics. 144: 1933-1940 

\href{https://www.genetics.org/content/206/4/2085}{Weir B.S. and Goudet J. (2017)} A Unified Characterization of Population Structure and Relatedness. Genetics. 206: 2085-2103

\href{https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14833}{Goudet J., Kay T. and Weir B.S. (2018)} How to estimate kinship. Molecular Ecology.  27: 4121:4135


Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

\href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1998.tb01824.x}{Yang, R.C. (1998)} Estimating hierarchical F-statistics. Evolution 52(4):950-956


}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\keyword{misc}
