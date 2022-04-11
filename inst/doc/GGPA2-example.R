### R code from vignette source 'GGPA-example.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: GGPAannot-prelim
###################################################
library("GGPA2")


###################################################
### code chunk number 3: ggpaExample-prelim
###################################################
data(simulation)
dim(simulation$pmat)
head(simulation$pmat)


###################################################
### code chunk number 4: pgraph-show (eval = FALSE)
###################################################
## adjmat <- simulation$true_G
## diag(adjmat) <- 0
## ggnet2( adjmat, label=TRUE, size=15 )


###################################################
### code chunk number 5: fig-pgraph
###################################################
adjmat <- simulation$true_G
diag(adjmat) <- 0
plot( adjmat, size=15 )


###################################################
### code chunk number 6: GGPAannot-show (eval = FALSE)
###################################################
## set.seed(12345)
## fit <- GGPA2( simulation$pmat )


###################################################
### code chunk number 7: GGPAannot-run
###################################################
set.seed(12345)
fit <- GGPA2( simulation$pmat, nBurnin=200, nMain=200 )


###################################################
### code chunk number 8: GGPAannot-show
###################################################
fit


###################################################
### code chunk number 9: GGPAannot-estimates
###################################################
str(estimates(fit))


###################################################
### code chunk number 10: GPA-assoc-ann
###################################################
assoc.marg <- assoc( fit, FDR=0.10, fdrControl="global" )
dim(assoc.marg)
apply( assoc.marg, 2, table )


###################################################
### code chunk number 11: GPA-fdr-ann
###################################################
fdr.marg <- fdr(fit)
dim(fdr.marg)
head(fdr.marg)


###################################################
### code chunk number 12: GPA-assoc-pattern-ann
###################################################
assoc.joint <- assoc( fit, FDR=0.10, fdrControl="global", i=1, j=2 )
length(assoc.joint)
head(assoc.joint)
table(assoc.joint)


###################################################
### code chunk number 13: GPA-pgraph-est-show (eval = FALSE)
###################################################
## plot(fit)


###################################################
### code chunk number 14: fig-pgraph
###################################################
plot(fit)


###################################################
### code chunk number 15: GGPAannot-pgraph-show (eval = FALSE)
###################################################
## fit <- GGPAannot( pmat, pgraph )


