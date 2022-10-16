# mecodev

A set of extended classes based on the microeco package

![](https://img.shields.io/badge/Test-Ver0.2.0-red.svg)


## Install

Install mecodev from github. Please first install microeco package https://github.com/ChiLiubio/microeco

```r
# If devtools package is not installed, first install it
install.packages("devtools")
devtools::install_github("ChiLiubio/mecodev")
```

If the installation from github is failed online, download the package first, then install it locally.

```r
devtools::install_local("mecodev-master.zip")
```

## Citation
Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, 
FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255. https://doi.org/10.1093/femsec/fiaa255

## Tutorial
See the mecodev part of tutorial (https://chiliubio.github.io/microeco_tutorial/mecodev-package.html).

## All classes

### trans_rarefy

The class trans_rarefy can be used for the rarefaction and the following plotting to 
see whether the sequencing depth is enough to cover all the taxa in the microbial community.

### trans_convert

The class trans_convert provide several data transformation approaches for the microtable object.


### trans_ts

The class trans_ts is designed for the time series data analysis.
Currently, the biomass estimation and biological interaction prediction based on the gLV model 
are implemented based on the beem package, Li et al. 2019. <doi: 10.1186/s40168-019-0729-z>.
The R package beem should be first installed from Github https://github.com/lch14forever/BEEM

```r
# For windows system:
install.packages("doMC", repos="http://R-Forge.R-project.org")
# For linux or mac
install.packages("doMC")
```

```r
# Then run
install.packages("lokern")
install.packages("monomvn")
install.packages("pspline")
devtools::install_github('csb5/beem')
```


### trans_gamma
The class trans_gamma is developed to explore the relationship between beta diversity and gamma diversity
based on the study Zhang et al.(2020).
The content includes the observed beta-gamma diversity relationship and simulated beta-gamma diversity relationship, and also the plotting.


## Contributing

We welcome any contribution \! 
Any idea/suggestion will be considered.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CONDUCT.md).


## References
  - Zhang et al. Local community assembly mechanisms shape soil bacterial β diversity patterns along a latitudinal gradient. Nat Commun 11, 5428 (2020); <doi:10.1038/s41467-020-19228-4>.
  - Liu et al. microeco: an R package for data mining in microbial community ecology, FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255.
  - Li et al. An expectation-maximization algorithm enables accurate ecological modeling using longitudinal microbiome sequencing data. Microbiome, (2019) 7:118; <doi:10.1186/s40168-019-0729-z>.


