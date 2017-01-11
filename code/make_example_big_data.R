library(ashr)
set.seed(100)
n=100000
z = rnorm(n,0,2)
#now sort z so that they are in order
z = z[order(abs(z))]

res <- ash(z,1,mixcompdist="normal",outputlevel=4)

save.image("../output/ash.bigz.RData")