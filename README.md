# mecodev

A set of extended classes based on the microeco package

![](https://img.shields.io/badge/Test-Ver0.1.0-red.svg)


## Install

Install mecodev from github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
# then install microeco; https://github.com/ChiLiubio/microeco
devtools::install_github("ChiLiubio/microeco")
# then install mecodev
devtools::install_github("ChiLiubio/mecodev")
```

If the installation from github is failed because of the bad internet, download the packages first, then install them locally.

```r
devtools::install_local("microeco-master.zip")
devtools::install_local("mecodev-master.zip")
```

## Citation
Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, 
FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255, https://doi.org/10.1093/femsec/fiaa255

## All classes

### trans_rarefy

The class trans_rarefy can be used for the rarefaction and the following plotting to 
see whether the sequencing depth is enough to cover all the taxa in the microbial community.

### trans_convert

The class trans_convert provide several data transformation approaches for the microtable object.

### trans_netchord 
The class trans_netchord is developed to sum and plot the links number from one taxa to another or in the same taxa in the network.
The input dataset must be a trans_network object.

The R package chorddiag is used for the chord plot in this class and can be installed from Github https://github.com/mattflor/chorddiag

```r
devtools::install_github("mattflor/chorddiag")
```


## Contributing

We welcome any contribution \! 
Any idea/suggestion will be considered.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CONDUCT.md).






