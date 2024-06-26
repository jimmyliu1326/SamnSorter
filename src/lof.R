# Helper script for calculating LOF for a single data point within an embedding
library(Rcpp)
library(RcppArmadillo)

euclidean_distance <- function(a, b) {
  sqrt(sum((a - b)^2))
}

k_distance <- function(data, point, k) {
  distances <- apply(data, 1, function(x) euclidean_distance(x, point))
  sorted_distances <- sort(distances)
  k_dist <- sorted_distances[k + 1]
  k_neighborhood <- which(distances <= k_dist)
  return(list(k_distance = k_dist, k_neighborhood = k_neighborhood))
}

reachability_distance <- function(data, point, k, k_neighborhood) {
  reach_distances <- sapply(k_neighborhood, function(i) {
    max(k_distance(data, data[i, ], k)$k_distance, euclidean_distance(data[i, ], point))
  })
  return(reach_distances)
}

local_reachability_density <- function(reach_distances) {
  return(1 / (mean(reach_distances)))
}

lof_score <- function(lrd, lrds) {
  return(mean(lrds) / lrd)
}

lof_point <- function(data_matrix, point, k = 5) {
  
  # calculate the k-distance and k-distance neighborhood
  kd <- k_distance(data_matrix, point, k)
  k_neighborhood <- kd$k_neighborhood
  k_dist <- kd$k_distance
  
  # calculate the reachability distance
  reach_distances <- reachability_distance(data_matrix, point, k, k_neighborhood)
  
  # calculate the local reachability density (LRD)
  lrd_point <- local_reachability_density(reach_distances)
  
  # calculate the LRD for the k-neighborhood
  lrds_neighbors <- sapply(k_neighborhood, function(i) {
    kd_neighbor <- k_distance(data_matrix, data_matrix[i, ], k)
    reach_distances_neighbor <- reachability_distance(data_matrix, data_matrix[i, ], k, kd_neighbor$k_neighborhood)
    local_reachability_density(reach_distances_neighbor)
  })
  
  # calculate the LOF score
  lof_value <- lof_score(lrd_point, lrds_neighbors)
  return(lof_value)
}