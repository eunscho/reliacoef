---
title: "reliacoef"
output:
  html_document:
    toc: yes
pagetitle: reliacoef
---

# Overview
<!-- badges: start -->
<!-- badges: end -->

  The goal of reliacoef is to calculate and compare various unidimensional and multidimensional
reliability coefficients.

+ Provides the following unidimensional reliability coefficients
  + coefficient alpha
  + Unidimensional confirmatory factor analysis reliability, commonly referred to as composite reliability
  + Gilmer-Feldt reliability coefficient
  + Feldt's classical congeneric reliability coefficient
  + Hancock's H
  + Heise-Bohrnstedt's Omega
  + Kaiser-Caffrey's alpha
  + Ten Berge and Zegers's mu series (mu2, mu3, mu4)
+ Provides the following multidimensional reliability coefficients, omega hierarchical, and subdimensional reliability
  + Stratified alpha
  + Maximal reliability
  + Multidimensional parallel reliability
  + Correlated factors reliability
  + Second-order factor reliability
  + Bifactor reliability
+ Collects and compare reliability coefficients provided by other packages
  + Two versions of GLB (greatest lower bounds) offered by the package psych
  + Guttman's lamdas offered by the package Lambda4
+ Test essential tau-equivalence and explore unidimensionality

# Installation

You can download and install it from Github using the devtools package:

``` r
install.packages("devtools")
devtools::install_github("eunscho/reliacoef")
```

# Example

The most typical use would be the unirel and multirel function comparing several reliability
coefficients:

``` r
library(unirel)
unirel(Graham1)
## compare various unidimensional reliability coefficients
multirel(Osburn_moderate, until = 4)
## compare various multidimensional reliability coefficients
```
You can also get each coefficient separately.
``` r
alpha(Graham1)
## obtain coefficient alpha
joreskog(Graham1)
## obtain composite (congeneric) reliability (unidimensional CFA reliability)
gilmer(Graham1)
## obtain the Gilmer-Feldt coefficient
feldt(Graham1)
## obtain Feldt's classical congeneric reliability
hancock(Graham1)
## obtain Hancock's H (maximal reliability)
heise(Graham1)
## obtain Heise-Borhnstedt's Omega
kaisercaffrey(Graham1)
## obtain Kaiser-Caffrey's alpha
mu2(Graham1)
## obtain Ten Berge and Zegers' mu2
mu3(Graham1)
## obtain Ten Berge and Zegers' mu3
mu4(Graham1)
## obtain Ten Berge and Zegers' mu4
stratified_alpha(Osburn_moderate, 4)
## obtain stratified alpha
multi_parallel(Osburn_moderate, 4)
## obtain multidimensional parallel reliability
second_order(Osburn_moderate, 4)
## obtain second-order factor reliability
bifactor(Osburn_moderate, 4)
## obtain bifactor reliability
maximal_reliability(Osburn_moderate, 4)
## obtain maximal reliability
correlated_factors(Osburn_moderate, 4)
## obtain correlated factors reliability
```
You can test essential tau-equivalence and explore unidimensionality.
``` r
test.tauequivalence(Graham1)
## test the assumption of essential tau-equivalence
```
# Troubleshooting

Sometimes an error message appears.

``` r
Error in standardizedsolution(fit) :
```
The solution is to activate the lavaan package.

``` r
library(lavaan)
```
