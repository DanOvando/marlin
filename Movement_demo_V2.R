library(tidyverse)
library(Matrix)
n_g = 100

# Domain characteristics
lat_g = seq(55, 65, length=n_g)
Temp_g = seq(15, 5, length=n_g)

# Parameters
diffusion = 0
preference_g =  (-2 * (Temp_g - 10)^2)

preference_g <-  preference_g - min(preference_g)

preference_g[42] <- 0

# plot
matplot( cbind(lat_g,Temp_g,preference_g), type="l", lty="solid", lwd=2 )

# Movement operator
A_gg = ifelse( round(abs(outer(lat_g,lat_g,"-")),2) == round(mean(diff(lat_g)),2), 1, 0 )

# Diffusion
diffusion_gg = A_gg * diffusion
diag(diffusion_gg) = -1 * colSums(diffusion_gg)

a = Matrix::expm(diffusion_gg) %>% as.matrix()

plot(a[n_g/2,])

# Taxis
# taxis_gg = A_gg * pmax(0,outer(preference_g, preference_g, "-"))

taxis_gg = A_gg * exp(outer(preference_g, preference_g, "-"))


# taxis_gg = A_gg * outer(preference_g, preference_g, "-")


diag(taxis_gg) = -1 * colSums(taxis_gg)
# Total
mrate_gg = diffusion_gg + taxis_gg
# Annualized
mfraction_gg = Matrix::expm(mrate_gg)
hist(mfraction_gg %>% as.matrix())
image(as.matrix(mfraction_gg))
rowSums(mfraction_gg)
colSums(mfraction_gg)


n <- matrix(ncol = 1, nrow = n_g,rep(0,n_g))

n[1,1] <-  100

for (y in 1:1000){

 n <-  mfraction_gg %*% n

}

# Stationary distribution

stationary_g = eigen(mfraction_gg)$vectors[,1]

stationary_g = stationary_g / sum(stationary_g)

matplot( x=lat_g, y=cbind(as.numeric(n / sum(n)),stationary_g), type="l", col=c("black","blue") )



# n(t+1) = Mrate_gg * n(t)
