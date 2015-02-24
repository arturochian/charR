#' Simlate Charcoal Records
#'
#' User a user-specified charcoal series to simulate new charcoal series that
#' are consistent with the user's series.
sim_char <- function(timevar, char, pois.offset = NULL,
                    front.cutoff = 2000, end.cutoff = 12500, threshold = 1,
                    polynomial = 10, phi = 0.085, sigma = NULL, size = NULL,
                    peak_freq = NA, ex = NA, vx = NA, pois.distn = TRUE,
                    seed = NULL, n.sims = 1) {
  # timevar is the estimated age of each observation
  # char is the charcoal accumulation rate of each observation
  # front.cutoff is the first age we'll accept observations from
  # end.cutoff is the last age we'll accept observations from
  # threshold is the cutoff point for background noise
  # polynomial is the order of the polynomial to use when fitting
  # phi is the estimate of the autocorrelation; defaults to value from 10th order polynomial
  # sigma is the standard deviation of the background noise
  ## if not supplied, defaults to being simulated from the fitted model
  # peak_freq is the rate at which we observe peaks
  ## if not supplied, defaults to the proportion of observations above the threshold
  # ex is the average CHAR in peak years
  ## if not supplied, defaults to average CHAR above threshold
  # vx is the variance of CHAR in peak years
  ## if not supplied, defaults to variance of CHAR above threshold
  # seed is an optional argument.  If specified, sets the seed so the simulation can be reproduced
  # n.sims is the number of simulations to run; defaults to 1

  # setting the seed
  if(!is.null(seed)) set.seed(seed)
  # setting up data.frame
  char.dat <- data.frame(timevar = timevar, CHAR = char)
  # removing data outside cutoffs
  char.include <- subset(char.dat, timevar > front.cutoff & timevar < end.cutoff)

  # cutting off peaks
  char.noise <- subset(char.include, CHAR < threshold)
  if(!is.null(pois.offset)) char.offset <- pois.offset[timevar > front.cutoff & timevar < end.cutoff & char < threshold]

  # setting up model matrix of background noise
  time_poly4 <- data.frame(cbind(char.noise$CHAR, 1, poly(char.noise$timevar, polynomial)))
  names(time_poly4)[1] <- "CHAR"
  n <- nrow(time_poly4)

  if(!pois.distn) {
    # fitting the polynomial with AR1 correlation
    poly_fit <- gls(CHAR ~ . - 1, data = time_poly4, correlation = corAR1())

    # subtracting off the number of coefficients and 2 (for sigma and phi)
    resid_df <- n - length(coef(poly_fit)) - 2

    ###### Simulating from Fitted Model #################

    ###  drawing variance-covariance matrix for simulated data

    # simulating sigma^2
    if(is.null(sigma)) sim_sigma2 <- poly_fit$sigma^2 * resid_df / rchisq(1, resid_df)
    else sim_sigma2 <- sigma^2

    # setting up correlation matrix for simulated data
    corr_mat <- matrix(NA, nrow = n, ncol = n)
    # comment later
    chron_lag <- Lag(char.noise$timevar)
    chron_diff <- char.noise$timevar - chron_lag
    chron_bar <- mean(chron_diff, na.rm = T)
    for(i in 1:n) {
      for(j in 1:n) {
        corr_mat[i, j] <- phi ^
          (abs(sum(chron_diff[1:i], na.rm = T) -
                 sum(chron_diff[1:j], na.rm = T)) / chron_bar)
      }
    }
    # multiplying to find covariance matrix
    sim_varcov <- sim_sigma2 * corr_mat

    ### drawing mean for simulated data

    # simulating coefficients for time trend from fitted model; using draw of sigma^2
    sim_beta2 <- mvrnorm(n.sims,
                         coef(poly_fit),
                         sim_sigma2 * summary(poly_fit)$corBeta)

    # multiplying to find mean for simulated data
    if(n.sims > 1) {
      mean4 <- apply(sim_beta2, 1, function(x) as.matrix(time_poly4[, -1]) %*% matrix(x))
    }
    else {
      mean4 <- as.matrix(time_poly4[, -1]) %*% matrix(sim_beta2)
    }


    # drawing CHAR series from MVN(mu = mean4, var = sim_varcov)
    sim_char4 <- apply(mean4, 2, function(x) {
      mvrnorm(1,
              x,
              sim_sigma2 * corr_mat) } )
    # setting negatives to zero
    sim_char4[sim_char4 < 0] <- 0
  } else {
    # fitting the model
    response <- round(char.offset * time_poly4$CHAR)
    poly_fit <- glm(response ~ . - 1 + offset(log(char.offset)),
                    data = time_poly4[, -1], family = quasipoisson())
    # simulating from model
    if(is.null(size)) size <- fitted(poly_fit) / (summary(poly_fit)$dispersion - 1)
    sim_beta2 <- arm::sim(poly_fit, n.sims = n.sims)@coef
    mean4 <- exp(model.matrix(poly_fit) %*% t(sim_beta2))
    sim_char4 <- apply(mean4, 2, function(mycol) rnbinom(length(mycol),
                                                         size = size,
                                                         mu = mycol*char.offset))
    sim_char4 <- apply(sim_char4, 2, function(mycol) mycol / char.offset)
  }

  ######  adding peaks ###################

  # determining peak frequency, mean and variance of peaks if not supplied
  if(is.na(peak_freq)) peak_freq <- 1 - nrow(char.noise) / nrow(char.include)
  if(is.na(ex)) ex <- mean(char.include$CHAR[char.include$CHAR > threshold])
  ex <- ex - threshold
  if(is.na(vx)) vx <- var(char.include$CHAR[char.include$CHAR > threshold])

  # drawing for which years will be peaks
  peak_year <- replicate(n.sims, rbinom(nrow(time_poly4), 1, peak_freq))

  # drawing values of CHAR in peak years
  alpha <- ex^2 / vx
  beta <- alpha / ex
  peak_char <- replicate(n.sims, rgamma(nrow(peak_year), alpha, beta) + threshold)


  # adding peaks to simulated data
  sim_char5 <- sim_char4
  for(i in 1:ncol(sim_char5)) {
    sim_char5[peak_year[, i] == 1, i] <- peak_char[peak_year[, i] == 1, i]
  }

  if(n.sims > 1) {
    sim.char <- data.frame(
      timevar = char.noise$timevar,
      sim_char5
    )
    char.melt <- reshape2::melt(sim.char, id.vars = "timevar", value.name = "CHAR", variable.name = "Simulation")

    sim.peaks <- data.frame(
      timevar = char.noise$timevar,
      peak_year
    )
    peak.melt <- reshape2::melt(sim.peaks, id.vars = "timevar", value.name = "peaks", variable.name = "Simulation")

    sim.data <- base::merge(char.melt, peak.melt)
    sim.data$peaks <- factor(sim.data$peaks)
    levels(sim.data$peaks) <- c("no.peak", "peak")
  } else {
    # setting up data.frame to return
    sim.data <- data.frame(
      timevar = char.noise$timevar,
      CHAR = sim_char5,
      peaks = factor(peak_year))
    levels(sim.data$peaks) <- c("no.peak", "peak")
  }

  return(list(data = sim.data, model = poly_fit,
              params = list(front.cutoff = front.cutoff, end.cutoff = end.cutoff,
                            threshold = threshold, polynomial = polynomial, phi = phi,
                            sigma = ifelse(exists("sigma2"), is.sqrt(sigma2), NA),
                            size = summary(size),
                            ex = ex, vx = vx, pois.distn = pois.distn, seed = seed, n.sims = n.sims)))
}
