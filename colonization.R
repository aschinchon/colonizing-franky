library(imager)
library(tidyverse)

# Location of the photograph
file <- "frankenstein.jpg"

# Load, convert to grayscale, filter image (to convert it to bw) and sample
load.image(file) %>%
  grayscale() %>%  
  threshold("45%") %>%
  as.cimg() %>%
  as.data.frame() %>%
  filter(value == 0) %>%
  select(x, y) -> franky  

# Parameters: play with them!
l     <- 25 # longitude of lines 
d_min <-  2 # minimun distance to search attractors
d_max <- 22 # maximun distance to search attractors
d_rem <- 12 # distance to remove attractors


# Sample size
n <- 3000

# Random sample of points from protograph
segments <- data.frame()

# We will compute 20 layers 
for(j in 1:20){
  
  # sample of points from photograph
  franky %>%
    sample_n(n) -> attractors
  
  # Initialization
  attractors %>%
    sample_n(1) %>%
    summarise(x = mean(x), y = mean(y)) %>% 
    mutate(parent = NA) %>%
    add_column(id = 1, .before = "x") -> nodes
  
  # let's colonize until less than 15000 attractors remain
  while(nrow(attractors) > 1500){
    
    # each attraction point is associated with the tree node that is closest to it
    merge(nodes, attractors, by = NULL, suffixes = c("_node","_attractor")) %>%
      mutate(d = sqrt((x_attractor - x_node)^2 + (y_attractor - y_node)^2)) -> distances_all
    
    # closest node to each attractor
    distances_all %>%
      group_by(x_attractor, y_attractor) %>%
      arrange(d) %>%
      slice(1) -> closests
    
    # keep those within of the radius of influence and normalize
    closests %>%
      ungroup() %>%
      filter(d > d_min, d < d_max) %>%
      mutate(x_dir = (x_attractor - x_node)/d,
             y_dir = (y_attractor - y_node)/d) %>%
      group_by(id, x_node, y_node) %>%
      summarise(x_dirs = sum(x_dir),
                y_dirs = sum(y_dir)) %>%
      mutate(x_norm = x_dirs/sqrt(x_dirs^2 + y_dirs^2),
             y_norm = y_dirs/sqrt(x_dirs^2 + y_dirs^2)) %>%
      ungroup() %>%
      transmute(x = x_node + sqrt(l)*cos(atan2(y_norm, x_norm)),
                y = y_node + sqrt(l)*sin(atan2(y_norm, x_norm)),
                parent = id) -> nodes_new
    
    # Add new nodes to the set of nodes
    add_column(nodes_new,
               id = seq(from = max(nodes$id)+1, by = 1, length.out = nrow(nodes_new)),
               .before = "x") %>% rbind(nodes) -> nodes
    
    # Remove closest attractors
    distances_all %>%
      filter(d < d_rem) %>%
      distinct(x = x_attractor, y = y_attractor) -> attractors_to_remove
    attractors %>%
      anti_join(attractors_to_remove, by = c("x", "y")) -> attractors
    
  }
  
  # Create segments
  nodes %>%
    inner_join(nodes, by = c("id" = "parent"), suffix = c("", "end")) %>%
    select(x, y, xend, yend) -> segments_tmp
  
  add_column(segments_tmp,
             drawing = j,
             .before = "x") -> segments_tmp
  
  segments %>% rbind(segments_tmp) -> segments  
}

# Draw result
ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, group = drawing),
               alpha = 0.4,
               data = segments) +
  coord_equal() +
  scale_y_reverse() +
  theme_void()
