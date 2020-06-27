set.seed(42)
library(tidyverse)

patches <- 25

ages <- 20

a = matrix(rnorm(patches * patches * ages,100), nrow = patches * patches, ncol = ages)

habitat <- expand_grid(x =1:patches, y = 1:patches) %>%
  mutate(habitat =  dnorm((x^2 + y^2), 2000,500))


habitat %>%
  ggplot(aes(x, y, fill = habitat)) +
  geom_tile()


distance <- expand_grid(x = 1:patches,y = 1:patches) %>%
  dist() %>%
  as.matrix()

distance <- dnorm(distance,0,12)

distance <- distance / rowSums(distance)


habitat_mat <-
  matrix(
    rep(habitat$habitat, patches),
    nrow = nrow(distance),
    ncol = ncol(distance),
    byrow = TRUE
  )

habitat_mat <- habitat_mat / rowSums(habitat_mat)

# image(habitat_mat)

dist_hab <-  distance * habitat_mat

dist_hab <- dist_hab / rowSums(dist_hab)

tt <- Sys.time()

b <- move_matrix_arma(dist_hab, a, 1000)

Sys.time() - tt


d <- array(0, dim = c(patches * patches,ages, 1000))

tt <- Sys.time()

for (i in 1:1000){


b <- move_matrix_arma(dist_hab, a, 1)

d[,,i] <-  b

}
Sys.time() - tt


