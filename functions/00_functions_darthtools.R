check_sum_of_transition_array <- function (a_P, n_rows = NULL, n_states = NULL, n_cycles, err_stop = TRUE, 
          verbose = TRUE) 
{
  if (!is.null(n_rows) & !is.null(n_states)) {
    stop("Pick either n_rows or n_states, not both.")
  }
  if (is.null(n_rows) & is.null(n_states)) {
    stop("Need to specify either n_rows or n_states, but not both.")
  }
  if (!is.null(n_rows)) {
    n_states <- n_rows
  }
  a_P <- as.array(a_P)
  d <- length(dim(a_P))
  if (d == 2) {
    valid <- sum(rowSums(a_P))
    if (valid != n_states) {
      if (err_stop) {
        stop("This is not a valid transition Matrix")
      }
      if (verbose) {
        warning("This is not a valid transition Matrix")
      }
    }
  }
  else {
    valid <- (apply(a_P, d, function(x) sum(rowSums(x))) == 
                n_states)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_cycles)))) {
      if (err_stop) {
        stop("This is not a valid transition Matrix")
      }
      if (verbose) {
        warning("This is not a valid transition Matrix")
      }
    }
  }
}



check_transition_probability <- function (a_P, err_stop = FALSE, verbose = FALSE) 
{
  a_P <- as.array(a_P)
  n_dim <- length(dim(a_P))
  if (n_dim < 3) {
    a_P <- array(a_P, dim = list(nrow(a_P), ncol(a_P), 1), 
                 dimnames = list(rownames(a_P), colnames(a_P), "Time independent"))
  }
  m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                 dim(a_P))
  if (dim(m_indices_notvalid)[1] != 0) {
    v_rows_notval <- rownames(a_P)[m_indices_notvalid[, 1]]
    v_cols_notval <- colnames(a_P)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 
                                                             3]]
    df_notvalid <- data.frame(`Transition probabilities not valid:` = matrix(paste0(paste(v_rows_notval, 
                                                                                          v_cols_notval, sep = "->"), "; at cycle ", v_cycles_notval), 
                                                                             ncol = 1), check.names = FALSE)
    if (err_stop) {
      stop("Not valid transition probabilities\n", paste(capture.output(df_notvalid), 
                                                         collapse = "\n"))
    }
    if (verbose) {
      warning("Not valid transition probabilities\n", paste(capture.output(df_notvalid), 
                                                            collapse = "\n"))
    }
  }
}


prob_to_odds <- function (p) 
{
  odds <- p/(1 - p)
  return(odds)
}


prob_to_prob <- function (p, t = 1) {
  if ((sum(p > 1) > 0) | (sum(p < 0) > 0)) {
    stop("probability not between 0 and 1")
  }
  p_new <- 1 - (1 - p)^(t)
  return(p_new)
}


rate_to_prob <- function (r, t = 1) {
  if ((sum(r < 0) > 0)) {
    stop("rate not greater than or equal to 0")
  }
  p <- 1 - exp(-r * t)
  return(p)
}

prob_to_rate <- function (p, t = 1){
  if ((sum(p > 1) > 0) | (sum(p < 0) > 0)) {
    stop("probability not between 0 and 1")
  }
  r = -(1/t) * log(1 - p)
  return(r)
}


get_os <- function () 
{
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") 
      os <- "MacOSX"
  }
  else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)) 
      os <- "osx"
    if (grepl("linux-gnu", R.version$os)) 
      os <- "linux"
  }
  tolower(os)
}
