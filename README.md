# hierfstat

The *hierfstat* package contains functions to estimate hierarchical F-statistics for any number of hierarchical levels using the method described in Yang (1998). It also contains functions allowing to test the significance of population differentiation at any given level using the likelihood ratio G-statistic, showed previoulsly to be the most powerful statistic to test for differnetiation (Goudet etal, 1996) . 

This is the development page of the *hierfstat* package for the R software.


## To install the development version of Hierfstat

You will need the package *devtools*  to be able to install the devel version of *hierfstat*. To install *devtools*:

```
install_github("devtools")
```

To install *hierfstat* devel:

```
install_github("jgx65/hierfstat")
library("hierfstat")
```


### References

Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

Goudet J., Raymond, M., DeMeeus, T. and Rousset F. (1996) Testing differentiation in diploid populations. Genetics. 144: 1933-1940

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52(4):950-956
