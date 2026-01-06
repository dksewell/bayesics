#' Create survival curves
#' 
#' Use the semi-parametric piecewise exponential survival model 
#' to fit a survival curve
#' 
#' 
#' @references 
#' Qing Y, Thall PF, Yuan Y. A Bayesian piecewise exponential phase II design for monitoring a time-to-event endpoint. Pharm Stat. 2023 Jan;22(1):34-44. doi: 10.1002/pst.2256. Epub 2022
#' 
#' @export


survfit_b = function(formula,
                     data,
                     prior_shape,
                     prior_rate,
                     max_n_time_bins,
                     n_time_bins){
  
  
  # Extract
  if(missing(data)){
    m = 
      model.frame(formula)
  }else{
    m = 
      model.frame(formula,data)
  }
  ## Get response variable
  Surv_object = 
    model.response(m)
  time = Surv_object[,1]
  status = Surv_object[,2]
  max_time = max(time)
  
  
  # If missing prior hyperparameters, use data to get reasonable guess based on t_i\sim exp(\lambda)
  if(missing(prior_shape) | missing(prior_rate)){
    # Based on posterior updating with unit information gain
    prior_shape = 
      0.001 + mean(status)
    prior_rate = 
      0.001 + mean(time)
  }
  
  
  ## Determine if multiple groups are to be analyzed
  single_group_analysis = (length(attr(terms(m),"term.labels")) == 0)
  
  
  # Create helper for computing posterior parameters
  # combined single group and multi-group into one to fix note  survfit_b: multiple local function definitions for ‘get_post_parms’
  # with different formal arguments
  get_post_parms = function(breakpoints, gr = NULL){
    
    if(single_group_analysis){
      bpt_assignment = 
        cut(time,
            breakpoints) |> 
        as.integer()
      n_j = 
        bpt_assignment[which(dplyr::near(status,1))] |> 
        factor(levels = 1:(length(breakpoints) - 1)) |> 
        table() |> 
        as.vector()
      bpt_diffs = diff(breakpoints)
      
      e_ij = matrix(0.0,length(time),length(breakpoints) - 1)
      for(j in 1:ncol(e_ij)){
        e_ij[which(bpt_assignment > j),j] = 
          bpt_diffs[j]
        
        happened_in_this_bin = 
          which(bpt_assignment == j)
        e_ij[happened_in_this_bin,j] = 
          time[happened_in_this_bin] -
          breakpoints[j]
      }
      
      cbind(a_j_tilde =
              prior_shape + n_j,
            b_j_tilde = 
              prior_rate + colSums(e_ij)
      ) 
      
    } else {
      
      bpt_assignment = 
        cut(time[group_assignments[[gr]]],
            breakpoints) |> 
        as.integer()
      n_j = 
        bpt_assignment[which(dplyr::near(status[group_assignments[[gr]]],1))] |> 
        factor(levels = 1:(length(breakpoints) - 1)) |> 
        table() |> 
        as.vector()
      bpt_diffs = diff(breakpoints)
      
      e_ij = matrix(0.0,length(group_assignments[[gr]]),length(breakpoints) - 1)
      for(j in 1:ncol(e_ij)){
        e_ij[which(bpt_assignment > j),j] = 
          bpt_diffs[j]
        
        happened_in_this_bin = 
          which(bpt_assignment == j)
        e_ij[happened_in_this_bin,j] = 
          time[group_assignments[[gr]]][happened_in_this_bin] -
          breakpoints[j]
      }
      
      cbind(a_j_tilde =
              prior_shape + n_j,
            b_j_tilde = 
              prior_rate + colSums(e_ij)
      ) 
      
    }
  }
  
  # Begin single group analysis
  if(single_group_analysis){
    
    
    # Create helper for computing posterior parameters
    # get_post_parms = function(breakpoints){
    #   bpt_assignment = 
    #     cut(time,
    #         breakpoints) |> 
    #     as.integer()
    #   n_j = 
    #     bpt_assignment[which(dplyr::near(status,1))] |> 
    #     factor(levels = 1:(length(breakpoints) - 1)) |> 
    #     table() |> 
    #     as.vector()
    #   bpt_diffs = diff(breakpoints)
    #   
    #   e_ij = matrix(0.0,length(time),length(breakpoints) - 1)
    #   for(j in 1:ncol(e_ij)){
    #     e_ij[which(bpt_assignment > j),j] = 
    #       bpt_diffs[j]
    #     
    #     happened_in_this_bin = 
    #       which(bpt_assignment == j)
    #     e_ij[happened_in_this_bin,j] = 
    #       time[happened_in_this_bin] -
    #       breakpoints[j]
    #   }
    #   
    #   cbind(a_j_tilde =
    #           prior_shape + n_j,
    #         b_j_tilde = 
    #           prior_rate + colSums(e_ij)
    #   ) 
    # }
    
    # Try out multiple breakpoints, select the optimal via marginal likelihood
    if(missing(n_time_bins)){
      trimmed_time = c(time[which(time < max_time)],max_time)
      
      if(missing(max_n_time_bins))
        max_n_time_bins = max(2,floor((NROW(time)/5)))
      
      ml_values = c(NA,numeric(max_n_time_bins - 1))
      for(J in 2:max_n_time_bins){
        bpts = quantile(trimmed_time,seq(0,1,l = 1 + J))
        bpts[1] = 0.0
        a_b_tilde = 
          get_post_parms(bpts)
        
        ml_values[J] = 
          sum(lgamma(a_b_tilde[,1])) -
          length(time) * lgamma(prior_shape) +
          length(time) * prior_shape * log(prior_rate) - 
          sum(a_b_tilde[,1] * log(a_b_tilde[,2]))
      }
      
      J_opt = which.max(ml_values)
    }else{
      J_opt = n_time_bins
    }
    
    # Get best (or pre-specified) breakpoints
    ret = list()
    bpts = quantile(trimmed_time,seq(0,1,l = 1 + J_opt))
    bpts[1] = 0.0
    ret$posterior_parameters = 
      get_post_parms(bpts)
    
    # Get intervals
    ret$intervals = 
      rep(bpts, 
          c(1,rep(2,length(bpts) - 2),1)) |> 
      matrix(ncol = 2,
             byrow = TRUE)
    
    # Get marginal likelihood
    ret$marginal_likelihood = 
      sum(lgamma(ret$posterior_parameters[,1])) -
      length(time) * lgamma(prior_shape) +
      length(time) * prior_shape * log(prior_rate) - 
      sum(ret$posterior_parameters[,1] * log(ret$posterior_parameters[,2]))
    
    # Return data also
    ret$data = m
    ret$single_group_analysis = single_group_analysis
    
    return(structure(ret,
                     class = "survfit_b"))
    
  }else{#End: single group analysis
    
    # Begin multiple group analysis
    
    # Ensure RHS are factors
    tt <- terms(formula)
    v_names = attr(tt, "term.labels")
    X <- m[, v_names, drop = FALSE]
    for(j in seq_along(X)){
      if(!is.factor(X[[j]])){
        X[[j]] <- factor(X[[j]])
      }
    }
    
    # Combine factors into a single factor
    factor_levels = 
      do.call(expand.grid,
              lapply(1:ncol(X),function(j) levels(X[[j]])))
    colnames(factor_levels) = colnames(X)
    factor_levels = 
      factor_levels |> 
      dplyr::mutate(
        group_names = 
          sapply(1:ncol(factor_levels),
                 function(j){
                   paste(colnames(X)[j],
                         factor_levels[[j]],
                         sep = ": ")
                 }) |> 
          apply(1,paste,collapse = "; ") |> 
          factor()
      )
    X = 
      dplyr::left_join(X,
                       factor_levels,
                       by = colnames(X))
    
    group_assignments = 
      split(seq_len(nrow(X)), X$group_names)
    G = length(group_assignments)
    
    
    # Create helper for computing posterior parameters
    # get_post_parms = function(breakpoints,gr){
    #   bpt_assignment = 
    #     cut(time[group_assignments[[gr]]],
    #         breakpoints) |> 
    #     as.integer()
    #   n_j = 
    #     bpt_assignment[which(dplyr::near(status[group_assignments[[gr]]],1))] |> 
    #     factor(levels = 1:(length(breakpoints) - 1)) |> 
    #     table() |> 
    #     as.vector()
    #   bpt_diffs = diff(breakpoints)
    #   
    #   e_ij = matrix(0.0,length(group_assignments[[gr]]),length(breakpoints) - 1)
    #   for(j in 1:ncol(e_ij)){
    #     e_ij[which(bpt_assignment > j),j] = 
    #       bpt_diffs[j]
    #     
    #     happened_in_this_bin = 
    #       which(bpt_assignment == j)
    #     e_ij[happened_in_this_bin,j] = 
    #       time[group_assignments[[gr]]][happened_in_this_bin] -
    #       breakpoints[j]
    #   }
    #   
    #   cbind(a_j_tilde =
    #           prior_shape + n_j,
    #         b_j_tilde = 
    #           prior_rate + colSums(e_ij)
    #   ) 
    # }
    
    # Try out multiple breakpoints, select the optimal via marginal likelihood
    J_opt = 
      a_b_tilde = 
      trimmed_time = 
      list()
    if(missing(n_time_bins))
      n_time_bins = NULL
    if(missing(max_n_time_bins))
      max_n_time_bins = NULL
    for(g in 1:G){
      if(is.null(n_time_bins)){
        trimmed_time[[g]] = 
          c(time[group_assignments[[g]]][which(time[group_assignments[[g]]] < max_time)],
            max_time)
        
        if(is.null(max_n_time_bins))
          max_n_time_bins = max(2,floor((NROW(time[group_assignments[[g]]])/5)))
        
        ml_values = c(NA,numeric(max_n_time_bins - 1) - Inf)
        for(J in 2:max_n_time_bins){
          bpts = quantile(trimmed_time[[g]],seq(0,1,l = 1 + J))
          if(anyDuplicated(bpts))
            next
          bpts[1] = 0.0
          a_b_tilde[[g]] = 
            get_post_parms(bpts,g)
          
          ml_values[J] = 
            sum(lgamma(a_b_tilde[[g]][,1])) -
            length(time) * lgamma(prior_shape) +
            length(time) * prior_shape * log(prior_rate) - 
            sum(a_b_tilde[[g]][,1] * log(a_b_tilde[[g]][,2]))
        }
        
        J_opt[[g]] = which.max(ml_values)
        
        max_n_time_bins = n_time_bins = NULL
      }else{
        J_opt[[g]] = n_time_bins
      }
    }
    
    
    # Get best (or pre-specified) breakpoints
    ret = list()
    for(g in 1:G){
      ret[[g]] = list()
    }
    bpts = 
      lapply(1:G,
             function(g){
               ret = 
                 quantile(trimmed_time[[g]],seq(0,1,l = 1 + J_opt[[g]]))
               ret[1] = 0.0
               return(ret)
             }
      )
    for(g in 1:G){
      ret[[g]]$posterior_parameters = 
        get_post_parms(bpts[[g]],g)
    }
    
    # Get intervals
    for(g in 1:G){
      ret[[g]]$intervals = 
        rep(bpts[[g]], 
            c(1,rep(2,length(bpts[[g]]) - 2),1)) |> 
        matrix(ncol = 2,
               byrow = TRUE)
    }
    
    names(ret) = 
      levels(factor_levels$group_names)
    
    # Get marginal likelihood
    ret$marginal_likelihood = 0.0
    for(g in 1:G){
      ret$marginal_likelihood = 
        ret$marginal_likelihood + 
        sum(lgamma(ret[[g]]$posterior_parameters[,1])) -
        length(group_assignments[[g]]) * lgamma(prior_shape) +
        length(group_assignments[[g]]) * prior_shape * log(prior_rate) - 
        sum(ret[[g]]$posterior_parameters[,1] * log(ret[[g]]$posterior_parameters[,2]))
    }
    
    
    # Return data also
    ret$data = m
    ret$single_group_analysis = single_group_analysis
    ret$group_names = levels(factor_levels$group_names)
    
    return(structure(ret,
                     class = "survfit_b"))
    
  }#End: multiple group analysis
  
}
