library(intermediate)
data(Tmem68, package = "Tmem68")
# Find and remove Tmem68 from mediators because it is target.
m <- match("Tmem68", Tmem68$annotation$symbol)
Tmem68$annotation[m,]
target = Tmem68$target
mediator = Tmem68$mediator[,-m]
driver = Tmem68$qtl.geno
#annotation = Tmem68$annotation[-m,]
covar = Tmem68$covar

fitDefault(driver, target,, covar[,-1])$LR
fitDefault(driver, target,, covar)$LR
fitQtl2(driver, target,, covar[,-1])$LR
fitQtl2(driver, target,, covar)$LR

# Below for fitQtl2 is what we want
driver8 <- cbind(A = 1 - apply(driver, 1, sum), driver)
fitDefault(driver8, target,, covar[,-1])$LR
fitDefault(driver8, target,, covar)$LR
fitQtl2(driver8, target,, covar[,-1])$LR
fitQtl2(driver8, target,, covar)$LR

