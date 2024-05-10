## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
set.seed(1234)
library(gsMAMS)

## -----------------------------------------------------------------------------
design_cont(delta0 = 0.178,delta1 = 0.545,alpha = 0.05,beta = 0.1,k=3,
frac = c(0.5,1))

## ----eval=FALSE---------------------------------------------------------------
#  op_fwer_cont(alpha=0.05,beta=0.1,p=3,frac=c(0.5,1),delta0=0.178,
#  delta1=0.545,nsim=10000,seed=10)
#  #> $FWER
#  #> [1] 0.05
#  #>
#  #> $`Stagewise FWER`
#  #>  look1  look2
#  #> 0.0050 0.0466
#  #>
#  #> $`Stopping probability under null`
#  #>  look1  look2
#  #> 0.2599 0.7401
#  #>
#  #> $`Probability of futility under null`
#  #>  look1  look2
#  #> 0.2549 0.6949
#  #>
#  #> $`Average sample size used per arm under null`
#  #> [1] 61.645

## ----eval=FALSE---------------------------------------------------------------
#  op_power_cont(alpha=0.05,beta=0.1,p=3,frac=c(0.5,1),delta0=0.178,
#  delta1=0.545,nsim=10000,seed=10)
#  #> $Power
#  #> [1] 0.893
#  #>
#  #> $`Stagewise Power`
#  #>  look1  look2
#  #> 0.3126 0.5804
#  #>
#  #> $`Stopping probability under alternative`
#  #>  look1  look2
#  #> 0.3258 0.6742
#  #>
#  #> $`Probability of futility under alternative`
#  #>  look1  look2
#  #> 0.0035 0.0821
#  #>
#  #> $`Average sample size used per arm under alternative`
#  #> [1] 62.652

## -----------------------------------------------------------------------------
design_ord(prob=c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166),or=3.06,
or0=1.32,alpha=0.05,beta=0.1,k=4,frac = c(1/3,2/3,1))

## -----------------------------------------------------------------------------
design_surv(m0=20,hr0=1, hr1=0.67032, ta=40,tf=20,alpha=0.05,beta=0.1,
k=4,kappa=1,eta=0,frac=c(0.5,1))

## ----eval=FALSE---------------------------------------------------------------
#  op_fwer_surv(m0=20,alpha=0.05,beta=0.1,p=4,frac=c(1/2,1),hr0=1,hr1=0.6703,
#  nsim=10000,ta=40,tf=20,kappa=1,eta=0,seed=12)
#  #> $FWER
#  #> [1] 0.05
#  #>
#  #> $`Stagewise FWER`
#  #>  look1  look2
#  #> 0.0049 0.0460
#  #>
#  #> $`Stopping probability under null`
#  #>  look1  look2
#  #> 0.2318 0.7682
#  #>
#  #> $`Probability of futility under null`
#  #>  look1  look2
#  #> 0.2269 0.7234
#  #>
#  #> $`Average number of events happened per arm under null`
#  #> [1] 129.5902
#  #>
#  #> $`Average duration of trial(months)`
#  #> [1] 54.27148

## ----eval=FALSE---------------------------------------------------------------
#  op_power_surv(m0=20,alpha=0.05,beta=0.1,p=4,frac=c(1/2,1),hr0=1,hr1=.6703,
#  nsim=10000,ta=40,tf=20,kappa=1,eta=0,seed=12)
#  #> $Power
#  #> [1] 0.913
#  #>
#  #> $`Stagewise Power`
#  #>  look1  look2
#  #> 0.3270 0.5863
#  #>
#  #> $`Stopping probability under alternative`
#  #>  look1  look2
#  #> 0.3334 0.6666
#  #>
#  #> $`Probability of futility under alternative`
#  #>  look1  look2
#  #> 0.0059 0.0800
#  #>
#  #> $`Average number of events happened per arm under alternative`
#  #> [1] 115.0483
#  #>
#  #> $`Average duration of trial(months)`
#  #> [1] 52.27229

