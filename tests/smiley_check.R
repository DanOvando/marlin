library(marlin)

library(tidyverse)

library(plot.matrix)

library(here)
#> here() starts at /private/var/folders/rs/t1c6lgz930g0k0cq98jsxdd00000gn/T/RtmphPRUF1/reprex-53fe7c907e0a-fishy-coral



nr <- 17
nc <- 12

bottom_value <- 1
top_value    <- 3 * bottom_value

# linear gradient values for each row (row 1 = top_value, row nr = bottom_value)
row_vals <- seq(from = top_value, to = bottom_value, length.out = nr)

# build matrix with same value across each row
m <- matrix(rep(row_vals, times = nc), nrow = nr, ncol = nc, byrow = FALSE)
colnames(m) <- paste0("V", seq_len(nc))

# helper: set cells to NA (row/col are 1-based)
set_na <- function(mat, rows, cols) {
  mat[cbind(rows, cols)] <- NA
  mat
}

eye_r0 <- 6
left_c0  <- 5
right_c0 <- 9

# Left eye: circle (3x3 ring)
ring_offsets <- expand.grid(dr = -1:1, dc = -1:1)
ring_offsets <- ring_offsets[!(ring_offsets$dr == 0 & ring_offsets$dc == 0), ]

m <- set_na(
  m,
  rows = eye_r0 + ring_offsets$dr,
  cols = left_c0 + ring_offsets$dc
)

# Right eye: X (3x3 diagonals + center)
x_offsets <- data.frame(
  dr = c(-1, -1,  0,  1, 1),
  dc = c(-1,  1,  0, -1, 1)
)

m <- set_na(
  m,
  rows = eye_r0 + x_offsets$dr,
  cols = right_c0 + x_offsets$dc
)

m <- set_na(m, rows = rep(12, 7), cols = 4:10)
m <- set_na(m, rows = rep(13, 5), cols = 5:9)
m <- set_na(m, rows = rep(14, 3), cols = 6:8)

habitat <- m

rm(m)

resolution <- rev(dim(habitat))

patch_area <- 2

years <- 20

seasons <- 4

time_step <- 1 / seasons

steps <- years * seasons


check <- habitat |>
  as.data.frame() |>
  mutate(y = n():1) |>
  tidyr::pivot_longer(-y, names_to = "x", values_to = "habitat") |>
  dplyr::mutate(x = match(x, unique(x))) |>
  dplyr::arrange(x,y)

check |>
  ggplot(aes(x,y,fill = habitat)) +
  geom_tile()


plot(habitat)

fauna <-
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      habitat = habitat,
      recruit_habitat = habitat,
      adult_home_range = 5,
      recruit_home_range = 15,
      density_dependence = "global_habitat",
      seasons = seasons,
      depletion = .25,
      resolution = resolution,
      patch_area = patch_area,
      steepness = 0.6,
      ssb0 = 4000,
      m = 0.1
    )
  )


check$b0 <- fauna$bigeye$r0s

check |>
  ggplot(aes(x,y,fill = b0)) +
  geom_tile()


n_in_space <- fauna$bigeye$grid

n_in_space$n <- rowSums(fauna$bigeye$unfished$n_p_a)

n_in_space |>
  ggplot(aes(x,y,fill = n)) +
  geom_tile()
