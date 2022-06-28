
#library(Matrix)
n_g = 10

# Domain characteristics
lat_g = seq(55, 65, length=n_g)
Temp_g = seq(15, 25, length=n_g)

# Parameters
diffusion = 0.8^2
preference_g = -0.1 * (Temp_g - 10)^2

# plot
matplot( cbind(lat_g,Temp_g,preference_g), type="l", lty="solid", lwd=2 )

# Movement operator
A_gg = ifelse( round(abs(outer(lat_g,lat_g,"-")),2) == round(mean(diff(lat_g)),2), 1, 0 )
 A_gg[1,2] <- 0
 A_gg[2,1] <- 0

# Diffusion
diffusion_gg = A_gg * diffusion
diag(diffusion_gg) = -1 * colSums(diffusion_gg)
# Taxis
taxis_gg = A_gg * outer(preference_g, preference_g, "-")
diag(taxis_gg) = -1 * colSums(taxis_gg)
# Total
mrate_gg = diffusion_gg + taxis_gg
# Annualized
mfraction_gg = Matrix::expm(mrate_gg)

# Stationary distribution
stationary_g = eigen(mfraction_gg)$vectors[,1]
stationary_g = stationary_g / sum(stationary_g)

#
matplot( x=lat_g, y=cbind(preference_g-min(preference_g),stationary_g), type="l", col=c("black","blue") )

# n(t+1) = Mrate_gg * n(t)
