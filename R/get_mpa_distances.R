get_mpa_distances <- function(patch_grid){
  
  open_patches <- patch_grid |> 
    filter(!mpa)
  
  mpa_patches <- patch_grid |> 
    filter(mpa)
  
  patch_distances <- patch_grid |> 
    select(x,y) |> 
    dist(diag = TRUE) |> 
    as.matrix()
  
  mpa_distances <- rowSums(patch_distances[open_patches$patch,mpa_patches$patch])
  
  nearest_mpa <- apply(patch_distances[open_patches$patch,mpa_patches$patch],1,min)
  
  fished_distances <- rowSums(patch_distances[mpa_patches$patch,open_patches$patch])
  
  nearest_fished <- apply(patch_distances[mpa_patches$patch,open_patches$patch],1,min)
  
  mpa_patches$total_distance_to_fished <- fished_distances
  
  mpa_patches$nearest_fished <- nearest_fished
  
  open_patches$total_distance_to_mpa <- mpa_distances
  
  open_patches$nearest_mpa <- nearest_mpa
  
  open_patches |> 
    ggplot(aes(x,y,fill = total_distance_to_mpa)) + 
    geom_tile() + 
    geom_tile(data = mpa_patches, aes(x,y), fill = "red")
  
  open_patches |> 
    ggplot(aes(x,y,fill = nearest_mpa)) + 
    geom_tile() + 
    geom_tile(data = mpa_patches, aes(x,y), fill = "red")
  
  mpa_patches |> 
    ggplot(aes(x,y,fill = total_distance_to_fished)) + 
    geom_tile() + 
    geom_tile(data = open_patches, aes(x,y), fill = "red")
  
  mpa_patches |> 
    ggplot(aes(x,y,fill = nearest_fished)) + 
    geom_tile() + 
    geom_tile(data = open_patches, aes(x,y), fill = "red")
  
  
  
  return(open_patches)
  
  
}