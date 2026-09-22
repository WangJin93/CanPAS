# synthetic cohort resembling CanPAS merged data
mk_cohort <- function(n = 200, prefix = "S", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- data.frame(ID = paste0(prefix, seq_len(n)),
                  OS_time = round(rexp(n, 0.15) + 0.1, 3),
                  OS_status = rbinom(n, 1, 0.5))
  d$GAPDH <- rnorm(n, 12, 1)
  d$TNS1  <- rnorm(n, 8, 1.2)
  d$PTEN  <- rnorm(n, 10, 1)
  d$age   <- round(runif(n, 30, 85))
  d$sex   <- sample(c("male", "female"), n, TRUE)
  d$stage <- factor(sample(c("I", "II", "III", "IV"), n, TRUE,
                           prob = c(0.4, 0.3, 0.2, 0.1)))
  d
}
