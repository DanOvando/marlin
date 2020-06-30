library(tidyverse)
library(FishLife)
library(mar)
library(Rcpp)
library(gganimate)
library(patchwork)
library(MASS)

fish <- create_critter(max_age = 40, r0 = 1e5)

# let's create some habitat

patches <- 25


habitat <- expand_grid(x =1:patches, y = 1:patches) %>%
  mutate(habitat =  dnorm((x^2 + y^2), 200,100))


habitat %>%
  ggplot(aes(x, y, fill = habitat)) +
  geom_tile()


distance <- expand_grid(x = 1:patches,y = 1:patches) %>%
  dist() %>%
  as.matrix()

distance <- dnorm(distance,0,5)

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

pop <- expand_grid(x = 1:patches,y = 1:patches) %>%
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


 dist_hab <-  distance #* habitat_mat

 dist_hab <- t(dist_hab / rowSums(dist_hab))

 a <- matrix(1, nrow = patches * patches, ncol = length(fish$length_at_age))

 a[,1] <- fish$r0 / (patches * patches)

 d = crossprod(dist_hab,adfs)

 test <-  1000
sum(a)
 tune_ssb0 <- sim_fish_pop(
   length_at_age = fish$length_at_age,
   weight_at_age = fish$weight_at_age,
   maturity_at_age = fish$maturity_at_age,
   steepness = 0.7,
   m = 0.2,
   patches = patches * patches,
   sim_steps = 1,
   burn_steps = 100,
   r0 = fish$r0,
   ssb0 =test,
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
    patches = patches * patches,
    sim_steps = 1,
    burn_steps = 0,
    r0 = fish$r0,
    ssb0 = tune_ssb0$ssb0,
    movement = dist_hab,
    tune_unfished = 0,
    last_n_p_a = tmp
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

tmplist$n_p_a[1,] %>% plot()

huh$ya <- rowSums(tmplist$ssb_p_a)




# huh$ya <- tmplist$n_p_a[,1]


huh %>%
  ggplot(aes(x, y, fill = ya)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "ya")

