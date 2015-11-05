# bottomup code

bottom_up_segment <- function(sequence, max_error) {
  seg_list <- list();

  for (i in 1:floor(length(sequence)/2)) {
    seg_list[[i]] <- list();
    seg_list[[i]]$indexes <- c(2*i-1, 2*i);
    seg_list[[i]]$values <- sequence[(2*i-1):(2*i)];
  }

  merge_costs <- rep(0, length(seg_list) - 1);

  for (i in 1:length(merge_costs)) {
    merged_seg <- merge_segs(seg_list[[i]], seg_list[[i+1]]);
    merge_costs[i] <- sum(lm(merged_seg$values ~ c(merged_seg$indexes[1]:merged_seg$indexes[2]))$residuals^2);
  }
  
  while (min(merge_costs) < max_error) {
    i <- which(merge_costs == min(merge_costs))[1];
    seg_list[[i]] <- merge_segs(seg_list[[i]], seg_list[[i+1]]);
    seg_list[[i+1]] <- NULL;
    
    first_part <- NULL;
    last_part <- NULL;
    
    if (i < length(seg_list)) {
      next_seg <- merge_segs(seg_list[[i]], seg_list[[i+1]]);
      merge_costs[i+1] <- sum(lm(next_seg$values ~ c(next_seg$indexes[1]:next_seg$indexes[2]))$residuals^2);
      last_part <- merge_costs[(i+1):length(merge_costs)];
    }
    
    if (i > 1) {
      prev_seg <- merge_segs(seg_list[[i-1]], seg_list[[i]]);
      merge_costs[i-1] <- sum(lm(prev_seg$values ~ c(prev_seg$indexes[1]:prev_seg$indexes[2]))$residuals^2);
      first_part <- merge_costs[1:i-1];
    }
    
    merge_costs <- c(first_part, last_part);    
  }
  
  res <- list();
  res$seg_list <- seg_list;
  res$merge_costs <- merge_costs;
  
  res
}

bottomupsegment_generic <- function(sequence, create_segment, compute_error, max_error) {  
  seg_list <- create_segment(sequence);
  
  merge_costs <- rep(0, length(seg_list) - 1);
  
  for (i in 1:length(merge_costs)) {
    merged_seg <- merge_segs(seg_list[[i]], seg_list[[i+1]]);
    merge_costs[i] <- compute_error(merged_seg$values, c(merged_seg$indexes[1], merged_seg$indexes[2])); 
  }
  
  while (min(merge_costs) < max_error) {
    i <- which(merge_costs == min(merge_costs))[1];
    seg_list[[i]] <- merge_segs(seg_list[[i]], seg_list[[i+1]]);
    seg_list[[i+1]] <- NULL;
    
    first_part <- NULL;
    last_part <- NULL;
    
    if (i < length(seg_list)) {
      next_seg <- merge_segs(seg_list[[i]], seg_list[[i+1]]);
      merge_costs[i+1] <- compute_error(next_seg$values, c(next_seg$indexes[1], next_seg$indexes[2]));
      last_part <- merge_costs[(i+1):length(merge_costs)];
    }
    
    if (i > 1) {
      prev_seg <- merge_segs(seg_list[[i-1]], seg_list[[i]]);
      merge_costs[i-1] <- compute_error(prev_seg$values, c(prev_seg$indexes[1], prev_seg$indexes[2]));
      first_part <- merge_costs[1:i-1];
    }
    
    merge_costs <- c(first_part, last_part);    
  }
  
  res <- list();
  res$seg_list <- seg_list;
  res$merge_costs <- merge_costs;
  
  res
}

create_segment_line <- function(sequence){
  seg_list <- list();
  
  for (i in 1:floor(length(sequence)/2)) {
    seg_list[[i]] <- list();
    seg_list[[i]]$indexes <- c(2*i-1, 2*i);
    seg_list[[i]]$values <- sequence[(2*i-1):(2*i)];
  }
  
  seg_list
}

compute_error_line <- function(sequence, indexes) {
  sum(lm(sequence ~ c(indexes[1]:indexes[2]))$residuals^2)
}

create_segment_poly2 <- function(sequence){
  seg_list <- list();
  
  for (i in 1:floor(length(sequence)/3)) {
    seg_list[[i]] <- list();
    seg_list[[i]]$indexes <- c(3*i-2, 3*i);
    seg_list[[i]]$values <- sequence[(3*i-2):(3*i)];
  }
  
  seg_list
}

compute_error_poly2 <- function(sequence, indexes) {
  sum(lm(sequence ~ c(indexes[1]:indexes[2])^2 + c(indexes[1]:indexes[2]))$residuals^2)
}

create_segment_poly3 <- function(sequence){
  seg_list <- list();
  
  for (i in 1:floor(length(sequence)/4)) {
    seg_list[[i]] <- list();
    seg_list[[i]]$indexes <- c(4*i-3, 4*i);
    seg_list[[i]]$values <- sequence[(4*i-3):(4*i)];
  }
  
  seg_list
}

compute_error_poly3 <- function(sequence, indexes) {
  sum(lm(sequence ~ c(indexes[1]:indexes[2])^3 + c(indexes[1]:indexes[2])^2 + c(indexes[1]:indexes[2]))$residuals^2)
}

merge_segs <- function(seg1, seg2) {
  new_seg <- list();
  new_seg$indexes <- c(seg1$indexes[1], seg2$indexes[2]);
  new_seg$values <- c(seg1$values, seg2$values);
    
  new_seg
}