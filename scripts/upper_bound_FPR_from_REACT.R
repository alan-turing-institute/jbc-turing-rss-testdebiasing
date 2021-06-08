d <- read.csv("data/region.csv")
d$phat <- d$nr / d$Nr
prev_pri_list <- list(unif_0_10 = seq(0, .1, length.out = 100), 
                      unif_0_5 = seq(0, .05, length.out = 100), 
                      unif_0_3 = seq(0, .03, length.out = 100), 
                      nonpar_eb = quantile(d$phat[!is.na(d$phat)], seq(0.01, .99, by = .01)))

par(mfrow = c(2, 2))
for(prinamc in names(prev_pri_list)) {
  prev_pri <- prev_pri_list[[prinamc]]
  pseq <- 10^seq(-5, -3, len = 100)
  prev <- 0
  pall <- 0
  for (prev in prev_pri) {
    pdens <- matrix(NA, nrow(d), length(pseq), dimnames = list(1:nrow(d), pseq))
    for(pc in pseq) {
      pdens[, as.character(pc)] <- dbinom(x = d$nr, size = d$Nr, prob = prev + pc, log = F)
    }  
    pall <- pall + pdens
  }
  tot_log_lik <- colSums(log(pall))
  plot(pseq, tot_log_lik, log = "x", ty = "l", xlab = "FPR (alpha)", ylab = "Log marginal likelihood")
  mtext(side = 3, text = paste0("Prevalence prior = ", prinamc), line = 2.5)
  labat <- pseq[which.max(tot_log_lik)]
  abline(v = labat)
  axis(side = 3, at = labat, label = signif(pseq[which.max(tot_log_lik)], 2))
}

