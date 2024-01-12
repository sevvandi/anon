# THIS FILE DOES THE FOLLOWING
# OPTIMIZATION PROBLEM = 1st ONE (0 \leq Su \leq f_u(d)
# IMPLEMENTATION DETAILS: dense matrix in lpSolve and
#                         standard solving in lpSolve
# WEIGHT SCHEMES:given in weights_opt
# PREDICTION/FORECASTING STEP: h in T+h
# PREDICTED NUMBER OF NODES: conf_level1
#                            (only consider increasing in size graphs)
# CONSTRAINTS: One constraint for degree of each node
#              One additional constraint for total number of edges
# OPTIMIZATION: Make the optimization tighter.
#               Remove columns with zero weights in the big matrix
#               because the zero weights mean no edge.
#               This will make a smaller optimization problem.
#               Remove rows with zero rhs.
# FORMULATION: Add input parameter for formulation 1 and formulation 2
#


predict_graph <- function(graphlist,
                          formulation = 2,
                          conf_level1 = NULL,
                          conf_level2 = 90,
                          dense_opt = 2,
                          weights_opt = 2,
                          weights_param = 0.001,
                          h = 1){
  # graphlist is the list of graphs
  # conf_level1 ==> this is to get different number of new nodes.
  #                 Right now, we're getting the mean predicted graph
  # conf_level2 ==> this is the upper bound of the vertex degree predictions
  #                 Right now, we're working only with the upper bound degree.
  # dense_opt   ==> this dictates of dense option needs to be used.
  #                 dense_opt = 1 for standard constraint matrix
  #                 dense_opt = 2 for dense matrix
  # weights_opt ==> this indicates the weights option
  #                 weights_opt = 1 is uniform weights
  #                 weights_opt = 2 is binary weights
  #                        1 for all except no-edges in past
  #                 weights_opt = 3 is binary weights (more selective version)
  #                        1 for all except no-edges and some 1s for new
  #                 weights_opt = 4 proportional weights

  grout_upper <- grout_lower <- graph_mean <- graph_lower <- graph_upper <- NULL

  # Step 1 - forecast the number of nodes
  new_nodes_list <-  predict_num_nodes(graphlist, conf_level1, h)
  new_nodes <- new_nodes_list$new_nodes
  if(!is.null(conf_level1)){
    # lower confidence level
    new_nodes_lower <- new_nodes_list$lower_conf
    # upper confidence level
    new_nodes_upper <- new_nodes_list$upper_conf
  }

  # Step 2 - Forecast the degree of the old nodes
  probj <- predict_old_nodes_degree(graphlist, conf_level2, h)

  # Step 3 - Using the above predict the mean graph
  grout <- predict_graph_internal(graphlist,
                                  formulation,
                                  conf_level1,
                                  conf_level2,
                                  dense_opt,
                                  weights_opt,
                                  weights_param,
                                  h,
                                  probj,
                                  new_nodes)

  # Step 4 - If conf_level1 is not NULL then predict the upper and lower graphs
  if(!is.null(conf_level1)){
    # lower confidence level
    grout_lower <- predict_graph_internal(graphlist,
                                          formulation,
                                          conf_level1,
                                          conf_level2,
                                          dense_opt,
                                          weights_opt,
                                          weights_param,
                                          h,
                                          probj,
                                          new_nodes_lower)
    # upper confidence level
    grout_upper <- predict_graph_internal(graphlist,
                                          conf_level1,
                                          formulation,
                                          conf_level2,
                                          dense_opt,
                                          weights_opt,
                                          weights_param,
                                          h,
                                          probj,
                                          new_nodes_upper)
  }

  list(
    graph_mean = grout,
    graph_lower = grout_lower,
    graph_upper = grout_upper
  )


}


predict_graph_internal <- function(graphlist,
                                   formulation,
                                   conf_level1 = NULL,
                                   conf_level2 = 90,
                                   dense_opt = 2,
                                   weights_opt = 2,
                                   weights_param = 0.25,
                                   h = 1,
                                   probj,
                                   new_nodes){
  # graphlist is the list of graphs
  # conf_level1 ==> this is to get different number of new nodes.
  #                 Right now, we're getting the mean predicted graph
  # conf_level2 ==> this is the upper bound of the vertex degree predictions
  #                 Right now, we're working only with the upper bound degree.
  # dense_opt   ==> this dictates of dense option needs to be used.
  #                 dense_opt = 1 for standard constraint matrix
  #                 dense_opt = 2 for dense matrix
  # weights_opt ==> this indicates the weights option
  #                 weights_opt = 1 is uniform weights
  #                 weights_opt = 2 is binary weights
  #                        1 for all except no-edges in past
  #                 weights_opt = 3 is binary weights (more selective version)
  #                        1 for all except no-edges and some 1s for new
  #                 weights_opt = 4 proportional weights


  # Step 1 - forecast the number of nodes - done in predict_graph function


  # Step 2 - Forecast the degree of the old nodes - done in predict_graph function
  # But the other hilo values need to be taken from probj
  f_rhs_hi <- probj$degree_hi
  f_rhs_lo <- probj$degree_lo
  freq <- probj$freq
  vertex <- probj$vertex
  total_hi <- probj$total_edges_upper
  total_low <- probj$total_edges_lower

  # Step 3 - Forecast the degree of new nodes
  # For new nodes forecast the average degree of new nodes
  rm_verts <- new_deg <- NULL
  rhsobj <- compute_rhs(graphlist, formulation, new_nodes, f_rhs_hi, f_rhs_lo = NULL, vertex, freq, total_hi, total_low)
  f_rhs <- rhsobj$f_rhs
  new_deg <- rhsobj$new_deg
  rm_verts <- rhsobj$rm_verts
  signs <- rhsobj$signs


  # Step 4 - Construct the union graph, find the constraint matrix and the weights
  biggr <- construct_union_graph(graphlist, new_nodes, new_deg, rm_verts, weights_opt = weights_opt, weights_param = weights_param)
  if(dense_opt == 1){
    f_con <- constraint_matrix(biggr, doubleup = FALSE)
  }else if(dense_opt == 2){
    f_con <- constraint_matrix_dense(biggr, doubleup = FALSE)
  }
  wts <- get_weights(biggr, f_con, weights_opt, dense_opt, doubleup = FALSE)
  # Step 4.1 - Remove zero weights and associated entries for optimization
  tighter_obj <- remove_zero_rows_and_columns(wts, f_rhs, f_con, f_obj, signs, dense_opt)
  f_con <- tighter_obj$f_con
  f_rhs <- tighter_obj$f_rhs
  signs <- tighter_obj$signs
  f_obj <- tighter_obj$f_obj


  if(dense_opt == 1){
    # Step 5 - Remove zero rows and Run lpSolve
    # Remove zero rows
    rhs_nonzero <- remove_zero_rhs(f_con, f_rhs, signs)
    f_con <- rhs_nonzero$f_con
    f_rhs <- rhs_nonzero$f_rhs
    signs <- rhs_nonzero$signs
    # # Remove >= rows. Because it didn't work for one example. NOT HAPPY ABOUT THIS!
    # rhs_geq <- remove_less_than(f_con, f_rhs, signs)
    # f_con <- rhs_geq$f_con
    # f_rhs <- rhs_geq$f_rhs
    # signs <- rhs_geq$signs
  }

  if(dense_opt == 1){
    lpobj <-  run_lpsolve(f_con, f_obj, f_rhs, signs)
  }else if(dense_opt == 2){
    lpobj <-  run_lpsolve_dense(f_con, f_obj, f_rhs, signs)
  }


  # Step 6: Get the graph relating to the values
  # Using the solution get the bit string - add the zeros that were there originally
  bitstring <- rep(0, length(wts))
  bitstring[which(wts!=0)] <- lpobj$solution

  if(dense_opt == 1){
    adj <- adjacency_from_solution(bitstring, gorder(biggr))
    grout <- graph_from_adjacency_matrix(adj, mode = "undirected")
  }else if (dense_opt == 2){
    edges <- edgelist_from_solution(bitstring, gorder(biggr))
    grout <- igraph::graph_from_edgelist(edges, directed = FALSE)
    add_verts <- vcount(biggr) - vcount(grout)
    if(add_verts > 0){
      grout <- igraph::add_vertices(grout, add_verts)
    }


  }


  grout

}

remove_zero_rows_and_columns <- function(wts, f_rhs, f_con, f_obj, signs, dense_opt){

  # Columns - keep only non-zero ones
  nonzero_col_inds <- which(wts != 0)
  f_obj <- wts[nonzero_col_inds]
  col_renamed_df <- data.frame(original = nonzero_col_inds, renamed = 1:length(nonzero_col_inds))

  # Rows - keep only non-zero ones
  nonzero_row_inds <- which(f_rhs !=0)
  f_rhs <- f_rhs[nonzero_row_inds]
  signs <- signs[nonzero_row_inds]
  row_renamed_df <- data.frame(original = nonzero_row_inds, renamed = 1:length(nonzero_row_inds))

  # Remove and rename the zero columns
  if(dense_opt == 1){
    f_con <- f_con[ ,nonzero_col_inds]
  }else if(dense_opt == 2){
    f_con <- f_con %>%
      filter(V3 !=0)
    f_con <- f_con %>%
      filter(V2 %in% col_renamed_df$original) %>%
      rename(original = V2) %>%
      left_join(col_renamed_df) %>%
      rename(V2 = renamed ) %>%
      select(-original) %>%
      select(V1, V2, V3)
  }


  # Remove and rename the zero rows
  if(dense_opt == 1){
    f_con <- f_con[nonzero_row_inds, ]
  }else if(dense_opt == 2){
    f_con <- f_con %>%
      filter(V1 %in% row_renamed_df$original) %>%
      rename(original = V1) %>%
      left_join(row_renamed_df) %>%
      rename(V1 = renamed) %>%
      select(-original) %>%
      select(V1, V2, V3)
  }

  list(
    f_con = f_con,
    f_rhs = f_rhs,
    signs = signs,
    f_obj = f_obj)
}


get_weights <- function(biggr, f_con, weights_opt, dense_opt, doubleup = FALSE){
  # The length of the objective function
  if(dense_opt == 1){
    len <- dim(f_con)[2]
  }else if(dense_opt == 2){
    len <- max(f_con[ ,2])
  }

  # weights schemes
  f_obj <- rep(0, len)
  if(weights_opt == 1){
    f_obj <- rep(1, len)
  }else if((weights_opt == 2)|(weights_opt == 3)){
    if(dense_opt == 1){
      # standard constraint matrix
      colsums <- apply(f_con, 2, sum)
      f_obj[which(colsums != 0)] <- 1
    }else if(dense_opt == 2){
      # dense constraint matrix
      colnames(f_con) <- c("constraint", "column", "value")
      cols <- f_con %>%
        filter(value > 0) %>%
        distinct(column) %>%
        pull(column)
      f_obj[cols] <- 1
    }
  }else if(weights_opt == 4){
    # Weights df
    weights <- cbind.data.frame(igraph::as_edgelist(biggr), E(biggr)$weight)
    colnames(weights) <- c("From", "To", "weight")
    weights <- weights %>%
      arrange(From, To)
    # Edges df as in the constraint matrix
    n <- vcount(biggr)
    from_seq <- c()
    to_seq <- c()
    for( i in 1:(n-1)){
      t_seq <- (i+1):n
      fr_seq <- rep(i, length(t_seq))
      from_seq <- c(from_seq, fr_seq)
      to_seq <- c(to_seq, t_seq)
    }
    df_edges <- data.frame("From" = from_seq, "To" = to_seq)
    f_obj <- left_join(df_edges, weights) %>%
      arrange(From, To) %>%
      tidyr::replace_na(list(weight = 0)) %>%
      pull(weight)
  }
  f_obj
}



remove_zero_rhs <- function(f_con, f_rhs, signs){
  inds <- which(f_rhs == 0)
  if(length(inds) > 0){
    f_rhs <- f_rhs[-inds]
    f_con <- f_con[-inds, ]
    signs <- signs[-inds]
  }
  list(
    f_con = f_con,
    f_rhs = f_rhs,
    signs = signs)

}

remove_less_than <- function(f_con, f_rhs, signs){
  inds <- which(signs == '>=')
  if(length(inds) > 0){
    f_rhs <- f_rhs[-inds]
    f_con <- f_con[-inds, ]
    signs <- signs[-inds]
  }
  list(
    f_con = f_con,
    f_rhs = f_rhs,
    signs = signs)

}

compute_rhs <- function(graphlist,
                        formulation,
                        new_nodes,
                        f_rhs_hi,
                        f_rhs_lo,
                        vertex,
                        freq,
                        total_hi,
                        total_lo){

  new_deg <- f_rhs <- rm_verts <- NULL

  if(new_nodes > 0){
    new_deg <- predict_new_nodes_degree(graphlist)
    if(is.null(f_rhs_lo)){
      f_rhs <- c(f_rhs_hi, rep(new_deg, new_nodes))
      signs <- c(rep("<=", length(f_rhs_hi)),  rep("<=", new_nodes) )
    }else{
      f_rhs <- c(f_rhs_hi, rep(new_deg, new_nodes), f_rhs_lo, rep(0, new_nodes))
      signs <- c(rep("<=", length(f_rhs_hi)),  rep("<=", new_nodes), rep(">=", length(f_rhs_lo)), rep(">=", new_nodes) )
    }
  }else if(new_nodes == 0){
    f_rhs <- c(f_rhs_hi, f_rhs_lo)
    signs <- c(rep("<=", length(f_rhs_hi)), rep(">=", length(f_rhs_lo)) )
  }else if(new_nodes < 0){
    rm_num <- -1*new_nodes
    signs <- c(rep("<=", length(f_rhs_hi)), rep(">=", length(f_rhs_lo)) )
    df_temp <- data.frame(vertex = vertex, degree = c(f_rhs_hi, f_rhs_lo), signs = signs, freq = freq)
    df_temp_rm <- df_temp %>%
      filter(degree == 0) %>%
      arrange(freq) %>%
      slice(1:rm_num)
    df_temp2 <- dplyr::setdiff(df_temp, df_temp_rm)
    f_rhs <- df_temp2 %>%
      pull(degree)
    rm_verts <- df_temp_rm %>%
      pull(vertex)
    signs <- df_temp2 %>%
      pull(signs)
  }

  # If only formulation = 2
  # Add the upper bound for total edges
  if(formulation == 2){
    f_rhs <- c(f_rhs, total_hi)
    signs <- c(signs, "<=")
  }

  list(
    new_deg = new_deg,
    f_rhs = f_rhs,
    rm_verts = rm_verts,
    signs = signs
  )
}

predict_num_nodes <- function(graphlist, conf_level = 90, h = 1){
  # Forecast the number of nodes
  NN <- length(graphlist)
  num_nodes <- rep(0, NN)
  lower_conf <- upper_conf <- NULL
  for(kk in 1:NN){
    num_nodes[kk] <- gorder(graphlist[[kk]])
  }

  df <- data.frame(index = 1:NN, nodes = num_nodes)

  if(NN < 15){
    lmmod <- lm(nodes ~index, data = df)
    next_nodes <- predict(lmmod, newdata = data.frame(index = (NN+h)))
  }else{
    next_nodes <- df %>% tsibble::as_tsibble(index = index) %>%
      fabletools::model(ets = fable::ETS(nodes)) %>%
      forecast(h = h) %>%
      pull(.mean)
    fit <- df %>% tsibble::as_tsibble(index = index) %>%
      fabletools::model(ets = fable::ETS(nodes)) %>%
      forecast(h = h)
    if(!is.null(conf_level)){
      dfhilo <- fit %>%
        hilo(level = conf_level)
      lower_conf <-  floor(unlist(dfhilo[ ,5])[1]) - gorder(graphlist[[NN]])
      upper_conf <-  ceiling(unlist(dfhilo[ ,5])[2]) - gorder(graphlist[[NN]])
    }
  }


  new_nodes <-  ceiling(next_nodes) - gorder(graphlist[[NN]])
  list(
    new_nodes = new_nodes[h],
    lower_conf = lower_conf,
    upper_conf = upper_conf
  )

}


construct_union_graph <- function(graphlist,
                                  new_nodes,
                                  new_deg,
                                  rm_verts = NULL,
                                  weights_opt = 2,
                                  weights_param = 0.001){
  # Construct the union graph
  # new_nodes = number of new nodes
  # new_deg =  average degree of new nodes
  # rm_verts = if new nodes is negative, then the set of nodes to remove
  # weights_opt = 2: Binary weights - new nodes connected to all old nodes
  # weights_opt = 3: Binary weights - new nodes connected to most connected old nodes
  # weights_opt = 4: Proportional weights

  NN <- length(graphlist)

  if(weights_opt == 4){
    # Graph needs weights
    for(ii in 1:NN){
      grr <- graphlist[[ii]]
      E(grr)$weight <- 1
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
        E(biggr)$weight_1[is.na(E(biggr)$weight_1)] <- 0
        E(biggr)$weight_2[is.na(E(biggr)$weight_2)] <- 0
        E(biggr)$weight <- E(biggr)$weight_1 +  E(biggr)$weight_2
      }
    }
    E(biggr)$weight <- E(biggr)$weight/max(E(biggr)$weight)
  }else{
    # For t = 1 to n
    for(ii in 1:NN){
      grr <- graphlist[[ii]]
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
      }
    }
  }

  if(new_nodes > 0){
    # Add new nodes
    biggr <- add_vertices(biggr, new_nodes)
    new_vertices <- V(biggr)[( gorder(biggr)-new_nodes+1):gorder(biggr)]
    old_vertices <- V(biggr)[1:(gorder(biggr)-new_nodes)]

    new_deg2 <- min(2*new_deg, vcount(biggr))
    verts <- order(degree(biggr), decreasing = TRUE)[1:new_deg2]

    if((weights_opt == 1)|(weights_opt == 2)){
      # Binary weights - new nodes connected to all old nodes
      # Add new edges from new nodes to all old nodes
      possible_edges <- c(rbind(rep(old_vertices, new_nodes), rep(new_vertices, each = length(old_vertices)) ))
    }else if(weights_opt == 3){
      # Binary weights - new nodes connected to most connected old nodes
      # Add new edges from new nodes to mostly connected vertices
      possible_edges <- c(rbind(rep(verts, new_nodes), rep(new_vertices, each = length(verts)) ))
    }else if(weights_opt == 4){
      # Proportional weights - new nodes connected to all old nodes
      # But the weights will be much smaller
      possible_edges <- c(rbind(rep(old_vertices, new_nodes), rep(new_vertices, each = length(old_vertices)) ))
      new_weights <- quantile(E(biggr)$weight, probs = weights_param)
    }

    if(weights_opt == 4){
      biggr <- biggr %>%
        add_edges(possible_edges,weight = new_weights)
    }else{
      biggr <- biggr %>%
        add_edges(possible_edges)
    }

  }else if(new_nodes < 0){
    # remove rm_verts from biggr
    biggr <- biggr %>%
      delete_vertices(rm_verts)
  }
  biggr
}


predict_old_nodes_degree <- function(graphlist, conf_level2, h){
  # Forecast the degree of the old nodes
  NN <- length(graphlist)
  for(jj in 1:NN){
    gr <- graphlist[[jj]]
    if(is.null(names(V(gr)))){
      v_names <- unclass(V(gr))
    }else{
      v_names <- names(V(gr))
    }
    df <- data.frame(time = jj, vertex = v_names, degree = degree(gr))
    df_sum <- data.frame(time = jj, edges = ecount(gr))
    if(jj == 1){
      dfall <- df
      dfall_sum <- df_sum
    }else{
      dfall <- rbind.data.frame(dfall, df)
      dfall_sum <- rbind.data.frame(dfall_sum, df_sum)
    }
  }

  dfmerge <- tibble::as_tibble(dfall) %>%
    tsibble::as_tsibble(index = time , key = vertex) %>%
    dplyr::arrange(time, vertex) %>%
    tsibble::fill_gaps(degree = 0)

  dffreq <- dfall %>%
    group_by(vertex) %>%
    summarize(freq = n())


  # Fit ARIMA models
  fit <- dfmerge %>%
    fabletools::model(arima = fable::ARIMA(degree),
                      naive = fable::NAIVE(degree)) %>%
    forecast(h = h) %>%
    full_join(dffreq)

  # Fit ARIMA for total edges
  fit_total <-  tibble::as_tibble(dfall_sum) %>%
    tsibble::as_tsibble(index = time) %>%
    fabletools::model(arima = fable::ARIMA(edges)) %>%
    forecast(h = h)

  # Get hilo for vertices separately
  dfhilo <- fit %>%
    hilo(level = conf_level2)
  uniq_times <- unique(dfhilo$time)

  colnames(dfhilo)[7] <- 'hilow'
  dfhilo <- dfhilo %>%
    mutate(lower = floor(hilow$lower), upper = ceiling(hilow$upper))
  dfhilo_vtx <- dfhilo %>%
    filter(time == uniq_times[h]) %>%
    group_by(vertex) %>%
    summarize(lower2 = max(min(lower, na.rm = TRUE), 0),
              upper2 = max(upper, na.rm = TRUE),
              mean2 = round(mean(abs(.mean), na.rm = TRUE)),
              freq = mean(freq))

  # Get hilo for total edges
  dfhilo_total <- fit_total %>%
    hilo(level = conf_level2)


  # Forecast of the old nodes
  # upper limit
  f_rhs_up <- dfhilo_vtx %>%
    pull(upper2)
  # lower limit
  f_rhs_lo <- dfhilo_vtx %>%
    pull(lower2)
  # mean
  f_rhs_mean <- dfhilo_vtx %>%
    pull(mean2)
  # frequency
  freq <- dfhilo_vtx %>%
    pull(freq)
  # vertex label
  vertex <- dfhilo_vtx %>%
    pull(vertex)

  inf_inds <- which(is.infinite(f_rhs_up))
  if(length(inf_inds) > 0){
    f_rhs_up[inf_inds] <- 3*f_rhs_mean[inf_inds]
  }

  #  Total edges predictions
  total_edges_mean <- dfhilo_total %>%
    pull(.mean) %>%
    ceiling()
  colnames(dfhilo_total)[5] <- 'hilow'
  dfhilo_total <- dfhilo_total %>%
    mutate(lower = floor(hilow$lower), upper = ceiling(hilow$upper))
  total_edges_lower <- dfhilo_total %>%
    pull(lower)
  total_edges_upper <- dfhilo_total %>%
    pull(upper)


  list(
    vertex = vertex,
    degree_hi = f_rhs_up,
    degree_lo = f_rhs_lo,
    degree_mean = f_rhs_mean,
    freq = freq,
    total_edges_mean = total_edges_mean[h],
    total_edges_lower = total_edges_lower[h],
    total_edges_upper = total_edges_upper[h]

  )

}


predict_new_nodes_degree <- function(graphlist){
  # Forecast the degree of new nodes
  # For new nodes forecast the average degree of new nodes
  NN <- length(graphlist)
  average_new_degree <-node_count <- edge_count <- rep(0, NN)
  for(jj in 1:NN){
    node_count[jj] <- vcount(graphlist[[jj]])
    edge_count[jj] <- ecount(graphlist[[jj]])
    new_node_count <- diff(node_count)
    if(jj == 1){
      average_new_degree[jj] <- 0
    }else{
      average_new_degree[jj] <- mean(degree(graphlist[[jj]])[(node_count[jj] - new_node_count[(jj-1)] + 1):node_count[jj]])
    }
  }
  df2 <- data.frame(index = 1:NN, avg_deg = average_new_degree)

  new_deg <- df2 %>% tsibble::as_tsibble(index = index) %>%
    fabletools::model(arima = fable::ARIMA(avg_deg)) %>%
    forecast(h = 1) %>%
    mutate(forecast = ceiling(.mean)) %>%
    pull(forecast)
  new_deg
}


run_lpsolve <- function(f_con, f_obj, f_rhs, signs){
  f_dir <- signs
  lpobj <- lpSolve::lp("max", f_obj, f_con, f_dir, f_rhs, all.bin = TRUE)
  lpobj
}


run_lpsolve_dense <- function(f_con, f_obj, f_rhs, signs){
  f_con <- as.matrix(f_con)
  f_dir <- signs
  lpobj <- lpSolve::lp(direction = "max",
                       objective.in = f_obj,
                       const.dir = f_dir,
                       const.rhs = f_rhs,
                       dense.const = f_con,
                       all.bin = TRUE)
  lpobj
}

generate_graph <- function(gr, del_edge, new_nodes, edge_increase = 0.1){
  # gr - graph to start with
  # del_edge - between 0 and 1. The proportion of edges to delete
  # new_nodes - if less than 1 then it is a proportion, else it is the number of nodes to add
  # edge_increase - the proportion of edges to add

  if(new_nodes < 1){
    new_nodes <- ceiling(vcount(gr)*(new_nodes))
  }
  num_edges <- ecount(gr)
  edges <- E(gr)

  # Removing edges
  num_edges1 <- floor(del_edge*num_edges)
  edges_remove <- sample(E(gr), num_edges1)
  gr2 <- igraph::delete_edges(gr, edges_remove)

  # Adding vertices
  gr2 <- igraph::add_vertices(gr2, new_nodes)

  # Add edges
  num_edges2 <- ceiling(num_edges*(1 + edge_increase)) -  ecount(gr2)
  gr2 <- gr2 + igraph::edge(sample(V(gr2), num_edges2*2, replace = T))
  gr2 <- igraph::simplify(gr2)
  return(gr2)
}

growth_metrics <- function(graphlist){
  n <- length(graphlist)
  deg_cosine_sim <- rep(0, (n-1))
  edges <- vertices <- rep(0, n)
  for( i in 1:n){
    vertices[i] <- vcount(graphlist[[i]])
    edges[[i]] <- ecount(graphlist[[i]])
  }
  for(i in 1:(n-1)){
    deg_i <- igraph::degree(graphlist[[i]])
    deg_ip1 <- igraph::degree(graphlist[[i+1]])[1:length(deg_i)]
    deg_cosine_sim[i] <- (sum(deg_i*deg_ip1))/(norm(as.matrix(deg_i), type = "F")*norm(as.matrix(deg_ip1), type = "F"))

  }
  vertex_growth <- vertices[2:n]/vertices[1:(n-1)] - 1
  edges_growth <- edges[2:n]/edges[1:(n-1)] - 1

  list(
    v_growth = vertex_growth,
    e_growth = edges_growth,
    v_growth_mean = mean(vertex_growth[11:(n-1)]),
    v_growth_sd = sd(vertex_growth[11:(n-1)]),
    e_growth_mean = mean(edges_growth[11:(n-1)]),
    e_growth_sd = sd(edges_growth[11:(n-1)]),
    deg_cosine_sim = deg_cosine_sim
  )
}


eval_metrics <- function(actualgr, predgr){

  new_vert <- NULL
  metric1 <- (vcount(predgr) - vcount(actualgr))/vcount(actualgr)

  if(metric1 > 0){
    # vcount(predgr) > vcount(actualgr)
    new_vert <- vcount(predgr) - vcount(actualgr)
    actualgr <-igraph::add_vertices(actualgr,nv =new_vert)
  }else if(metric1 < 0){
    # vcount(predgr) < vcount(actualgr)
    new_vert <-  vcount(actualgr) - vcount(predgr)
    predgr <-igraph::add_vertices(predgr,nv =new_vert)
  }
  adj_act <-igraph::as_adj(actualgr)
  adj_pred <- igraph::as_adj(predgr)
  metric2 <- (Matrix::norm((adj_pred - adj_act), type = "F"))/(Matrix::norm(adj_act, type = "F"))

  metric3 <- diff_metrics(as.vector(adj_act), as.vector(adj_pred))

  deg_act <- igraph::degree(actualgr)
  deg_prd <- igraph::degree(predgr)
  deg_cosine_sim <- (sum(deg_act*deg_prd))/(norm(as.matrix(deg_act), type = "F")*norm(as.matrix(deg_prd), type = "F"))
  # lsa::cosine(deg_act, deg_prd)

  lap_act <- igraph::laplacian_matrix(actualgr)
  lap_prd <- igraph::laplacian_matrix(predgr)

  eig_act <- eigen(lap_act)$values
  eig_prd <- eigen(lap_prd)$values
  eig_cosine_sim <- (sum(eig_act*eig_prd))/(norm(as.matrix(eig_act), type = "F")*norm(as.matrix(eig_prd), type = "F"))


  list(
    node_error = metric1,
    adjacency_error = metric2,
    degree_cosine = deg_cosine_sim,
    lap_eigen_cosine = eig_cosine_sim,
    precision = metric3$precision,
    recall = metric3$recall,
    fmeasure = metric3$fmeasure,
    sensitivity = metric3$sensitivity,
    specificity = metric3$specificity,
    gmean = metric3$gmean,
    accuracy = metric3$accuracy

  )
}

diff_metrics <- function(act, pred) {
  # positives to be denoted by 1 and negatives with 0
  stopifnot(length(act) == length(pred))
  n <- length(act)
  tp <- sum((act == 1) & (pred == 1))
  tn <- sum((act == 0) & (pred == 0))
  fp <- sum((act == 0) & (pred == 1))
  fn <- sum((act == 1) & (pred == 0))
  prec <- (tp + tn) / n
  sn <- tp / (tp + fn)
  sp <- tn / (tn + fp)
  precision <- if_else(
    (tp + fp) == 0,
    0,
    tp / (tp + fp)
  )
  recall <- tp / (tp + fn)
  fmeasure <- if_else(
    (precision == 0) & (recall == 0),
    0,
    2 * precision * recall / (precision + recall)
  )
  tibble(
    N = n,
    true_pos = tp,
    true_neg = tn,
    false_pos = fp,
    false_neg = fn,
    accuracy = prec,
    sensitivity = sn,
    specificity = sp,
    gmean = sqrt(sn * sp),
    precision = precision,
    recall = recall,
    fmeasure = fmeasure
  )
}


constraint_matrix <- function(gr, doubleup = FALSE){
  # gr is an igraph
  adj <- as_adjacency_matrix(gr)
  num_cols <- dim(adj)[1]*(dim(adj)[1] - 1)/2
  NN <- num_rows <- dim(adj)[1]
  constr_mat <- matrix(0, nrow = num_rows, ncol = num_cols)
  st_col <- 1
  for(jj in 1:(num_rows-1)){
    vals <- adj[jj,(jj+1):NN]
    en_col <- st_col + length(vals) - 1
    # Update relevant row in the constraint matrix
    constr_mat[jj, st_col:en_col] <- vals
    # Fill the diagonal in the submatrix below
    if(dim(diag(vals))[1] == 0){
      # Last entry - no matrix, only 1 entry
      constr_mat[(jj+1):NN,st_col:en_col] <- vals
    }else{
      constr_mat[(jj+1):NN,st_col:en_col] <-  diag(vals)
    }

    # Update st_col
    st_col <- en_col + 1

  }

  # All edges constraint
  all_edges <- rep(1, dim(constr_mat)[2])
  constr_mat <- rbind.data.frame(constr_mat, all_edges)

  constr_mat
}


constraint_matrix_dense <- function(gr, doubleup = FALSE){
  # gr is an igraph
  edges <- igraph::as_edgelist(gr)
  colnames(edges) <- c("From", "To")
  n <- vcount(gr)
  n_seq <- c(0,seq(n-1, 1))
  n_seq2 <- cumsum(n_seq)
  for( i in 1:n){
    vertex <- i
    dat_edges <- edges %>%
      as.data.frame() %>%
      dplyr::filter(From == i)
    if(NROW(dat_edges) > 0){
      df <- matrix(0, nrow = 2*NROW(dat_edges), ncol = 3)
      to_edges <- dat_edges[ ,2]
      cols <- n_seq2[i] + to_edges-i
      len <- NROW(dat_edges)
      # The first part of the constraints
      # en <- st + len - 1
      df[1:len, 1] <- i
      df[1:len, 2] <- cols
      df[1:len, 3] <- 1
      df[(len+1):(2*len), 1] <- to_edges
      df[(len+1):(2*len), 2] <- cols
      df[(len+1):(2*len), 3] <- 1
    }else{
      df <- matrix(c(i, n_seq2[i], 0), nrow = 1)
    }
    if(i == 1){
      dense_mat <- df
    }else{
      dense_mat <- rbind.data.frame(dense_mat, df)
    }
  }

  # Total edges constraint
  num_rows <- max(dense_mat[ ,2])
  edges_const <- data.frame(constraint = (n+1), coef = 1:num_rows, val = 1)
  colnames(edges_const) <- colnames(dense_mat)
  dense_mat <- rbind.data.frame(dense_mat, edges_const)

  dense_mat <- dense_mat %>%
    as.data.frame() %>%
    arrange(V1, V2)

  dense_mat

}


constraint_matrix_with_newnodes <- function(gr, newnodes){
  # gr is an igraph
  # newnodes is the number of new nodes
  adj <- as_adjacency_matrix(gr)
  num_cols <- dim(adj)[1]*(dim(adj)[1] - 1)/2
  NN <- num_rows <- dim(adj)[1]
  constr_mat <- matrix(0, nrow = num_rows, ncol = num_cols)
  new_positions <- rep(0, num_cols)
  st_col <- 1
  for(jj in 1:(num_rows-1)){
    vals <- adj[jj,(jj+1):NN]
    en_col <- st_col + length(vals) - 1
    # Update relevant row in the constraint matrix
    constr_mat[jj, st_col:en_col] <- vals
    # Fill the diagonal in the submatrix below
    if(dim(diag(vals))[1] == 0){
      # Last entry - no matrix, only 1 entry
      constr_mat[(jj+1):NN,st_col:en_col] <- vals
    }else{
      constr_mat[(jj+1):NN,st_col:en_col] <-  diag(vals)
    }

    # New positions vector
    new_positions[(en_col-newnodes + 1):en_col] <- 1

    # Update st_col
    st_col <- en_col + 1

  }
  list(
    constr_mat = constr_mat,
    new_positions = new_positions)

}



adjacency_from_solution <- function(x, n){
  # x is the solution vector
  # n is the number of vertices
  adj <- matrix(0, nrow = n, ncol = n)
  adj[lower.tri(adj)] <- x
  adj2 <- t(adj)
  adj2[lower.tri(adj2)] <- x
  adj2
}

edgelist_from_solution <- function(x, n){
  # x is the solution vector
  # n is the number of vertices
  v_seq <- 1:n
  for(i in 1:(n-1)){
    v1 <- rep(i,(n-i))
    v2 <- (i+1):n
    df <- cbind.data.frame(v1, v2)
    if(i == 1){
      dfall <- df
    }else{
      dfall <- rbind.data.frame(dfall, df)
    }
  }
  edges <-as.matrix(dfall[which(x == 1), ])
  edges

}


eval_metrics2 <- function(actualgr, predgr){

  new_vert <- NULL
  metric1 <- (vcount(predgr) - vcount(actualgr))/vcount(actualgr)

  if(metric1 > 0){
    # vcount(predgr) > vcount(actualgr)
    new_vert <- vcount(predgr) - vcount(actualgr)
    actualgr <-igraph::add_vertices(actualgr,nv =new_vert)
  }else if(metric1 < 0){
    # vcount(predgr) < vcount(actualgr)
    new_vert <-  vcount(actualgr) - vcount(predgr)
    predgr <-igraph::add_vertices(predgr,nv =new_vert)
  }
  adj_act <-igraph::as_adj(actualgr)
  adj_pred <- igraph::as_adj(predgr)

  # arrange rows by row sums and columns by column sums
  adj_act2 <- arrange_matrix(adj_act)
  adj_pred2 <- arrange_matrix(adj_pred)
  metric2 <- (Matrix::norm((adj_pred2 - adj_act2), type = "F"))/(Matrix::norm(adj_act2, type = "F"))

  metric3 <- diff_metrics(as.vector(adj_act), as.vector(adj_pred))

  deg_act <- igraph::degree(actualgr)
  deg_prd <- igraph::degree(predgr)
  deg_cosine_sim <- (sum(deg_act*deg_prd))/(norm(as.matrix(deg_act), type = "F")*norm(as.matrix(deg_prd), type = "F"))
  # lsa::cosine(deg_act, deg_prd)

  lap_act <- igraph::laplacian_matrix(actualgr)
  lap_prd <- igraph::laplacian_matrix(predgr)

  eig_act <- eigen(lap_act)$values
  eig_prd <- eigen(lap_prd)$values
  eig_cosine_sim <- (sum(eig_act*eig_prd))/(norm(as.matrix(eig_act), type = "F")*norm(as.matrix(eig_prd), type = "F"))


  list(
    node_error = metric1,
    adjacency_error = metric2,
    degree_cosine = deg_cosine_sim,
    lap_eigen_cosine = eig_cosine_sim,
    precision = metric3$precision,
    recall = metric3$recall,
    fmeasure = metric3$fmeasure,
    sensitivity = metric3$sensitivity,
    specificity = metric3$specificity,
    gmean = metric3$gmean,
    accuracy = metric3$accuracy

  )
}


arrange_matrix <- function(mat){
  rowsums <- apply(mat, 1, sum)
  ord <- order(rowsums, decreasing = TRUE)
  mat2 <- mat[ord, ]
  colsums <- apply(mat, 2, sum)
  ord <- order(colsums, decreasing = TRUE)
  mat3 <- mat2[ ,ord]
  mat3
}

