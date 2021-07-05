library(intermediate)
data(Tmem68)
# Find and remove Tmem68 from mediators because it is target.
m <- match("Tmem68", Tmem68$annotation$symbol)
Tmem68$annotation[m,]
target = Tmem68$target
mediator = Tmem68$mediator[,-m]
driver = Tmem68$driver
#annotation = Tmem68$annotation[-m,]
covar = Tmem68$covar
colnames(covar)[1:2] <- c("Intercept", "SexM")

fitDefault(driver, target, covar[,-1, drop = FALSE])$LR
fitDefault(driver, target, covar)$LR

# Below for fitQtl2 is what we want
driver8 <- cbind(A = 1 - apply(driver, 1, sum), driver)
fitDefault(driver8, target, covar[,-1, drop = FALSE])$LR
fitDefault(driver8, target, covar)$LR

