library(tidyverse)
library(FishLife)
library(mar)
library(Rcpp)
library(gganimate)
library(patchwork)
library(MASS)
library(doParallel)
options(dplyr.summarise.inform = FALSE)

# let's create some habitat

resolution <- 25


habitat <- expand_grid(x =1:resolution, y = 1:resolution) %>%
  mutate(habitat =  dnorm((x^2 + y^2), 600,100))


habitat %>%
  ggplot(aes(x, y, fill = habitat)) +
  geom_tile()


distance <- expand_grid(x = 1:resolution,y = 1:resolution) %>%
  dist() %>%
  as.matrix()

distance <- dnorm(distance,0,2)

distance <- distance / rowSums(distance)

wtf <- create_critter(adult_movement = 2)


habitat_mat <-
  matrix(
    rep(habitat$habitat, resolution),
    nrow = nrow(distance),
    ncol = ncol(distance),
    byrow = TRUE
  )

habitat_mat <- habitat_mat / rowSums(habitat_mat)

# habitat <- habitat_mat


image(distance)

image(habitat_mat)

dist_hab <-  distance * habitat_mat

dist_hab <- dist_hab / rowSums(dist_hab)

pop <- expand_grid(x = 1:resolution,y = 1:resolution) %>%
  mutate(n = nrow(.):1)

pop %>%
  ggplot() +
  geom_tile(aes(x, y, fill = n)) +
  scale_fill_viridis_c()

tt <- Sys.time()
pops <- list()
for (i in 1:1000){

  pop$n <- as.vector(pop$n %*% dist_hab)

  pop$i <- i

  # pop %>%
  #   ggplot() +
  #   geom_tile(aes(x, y, fill = n)) +
  #   scale_fill_viridis_c()
#
  pops[[i]] <- pop

}

Sys.time() - tt

pops <- bind_rows(pops, .id = "i") %>%
  mutate(i = as.numeric(i))

# pops %>%
#   ggplot() +
#   geom_tile(aes(x, y, fill = n)) +
#   scale_fill_viridis_c() +
#   transition_states(i) +
#   labs(title = 'i: {frame}')


huh <- pops %>%
  filter(i == max(i)) %>%
  mutate(habitat = habitat$habitat)


b <- huh %>%
  ggplot() +
  geom_tile(aes(x, y, fill = n)) +
  scale_fill_viridis_c() +
  labs(title = "loopmagic")

h <- habitat %>%
  ggplot(aes(x, y, fill = habitat)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "habitat")


 b + h

 fish <- create_critter(max_age = 20, r0 = 1e5)

 dist_hab <-  distance * habitat_mat

 dist_hab <- t(dist_hab / rowSums(dist_hab))

 a <- matrix(1, nrow = resolution * resolution, ncol = length(fish$length_at_age))

 a[,1] <- fish$r0 / (resolution * resolution)

 test <-  1000
sum(a)
 tune_ssb0 <- sim_fish(
   length_at_age = fish$length_at_age,
   weight_at_age = fish$weight_at_age,
   maturity_at_age = fish$maturity_at_age,
   steepness = 0.7,
   m = 0.2,
   patches = resolution * resolution,
   burn_steps = 100,
   r0 = fish$r0,
   ssb0 = NA,
   movement = dist_hab,
   last_n_p_a = a,
   tune_unfished = 1
 )
sum(a)


 ssb0 <- tune_ssb0$ssb0

tmp <- a

tt <- Sys.time()
tmplist <- list()
for (i in 1:100) {
  tmplist[[i]] <- sim_fish_pop(
    length_at_age = fish$length_at_age,
    weight_at_age = fish$weight_at_age,
    maturity_at_age = fish$maturity_at_age,
    steepness = 0.7,
    m = 0.2,
    patches = resolution * resolution,
    burn_steps = 0,
    r0 = fish$r0,
    ssb0 = tune_ssb0$ssb0,
    movement = dist_hab,
    last_n_p_a = tmp,
    tune_unfished = 0
  )

  tmp <- tmplist[[i]]$n_p_a
}
Sys.time() - tt

rec <- map(tmplist, "n_p_a") %>%
  map_df(~tibble(rec = .x[,1]),.id = "i") %>%
  mutate(i = as.numeric(i)) %>%
  group_by(i) %>%
  summarise(recs = sum(rec))

ssb <- map(tmplist, "ssb_p_a") %>%
  map_df(~tibble(ssb = rowSums(.x)),.id = "i") %>%
  mutate(i = as.numeric(i)) %>%
  group_by(i) %>%
  summarise(ssb = sum(ssb))

plot(ssb$ssb)

last(ssb$ssb) / ssb0

plot(ssb$ssb, rec$recs)


huh$ya <- rowSums(tmplist[[i]]$ssb_p_a)


# huh$ya <- tmplist[[i]]$n_p_a[,1]


huh %>%
  ggplot(aes(x, y, fill = ya)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "ya")

