library(parallel)
library(foreach)
library(doParallel)

source("scripts/SIR_utils.R")
control_debias <- prevdebiasr::get_control_parameters()

npt.inf <- 50
eps <- 1e-10
ncore.use <- 15

bin.d <- control_debias$bin.d
mseqc <- control_debias$I_seq
mseq.all <- 0:(2 * max(mseqc))
binmap <- data.frame(m = mseq.all, bin = findInterval(mseq.all, bin.d$e.l))
binmap$e.l <- bin.d[binmap$bin, "e.l"]
binmap$e.u <- bin.d[binmap$bin, "e.u"]
binmap$mu <- bin.d[binmap$bin, "mu"]
if("cl" %in% ls())
  stopCluster(cl)
cl <- makeCluster(ncore.use)
registerDoParallel(cl)
trmat <- matrix(0, length(control_debias$I_seq), length(control_debias$I_seq), 
                dimnames = list(control_debias$I_seq, control_debias$I_seq))
temp_dir <- tempdir()

trl <- foreach(Rc = control_SIR$R_grid, 
               .errorhandling = "stop", 
               .verbose = T) %dopar% {
                 print(Rc)
                 file_out <- file.path(temp_dir, paste0(Rc, ".RDS"))
                 replace_files <- F
                 if(!file.exists(file_out) | replace_files) {
                   for(m1 in mseqc){
                     print(m1)
                     inds.prob.gz <- 1 + bin.d[match(m1, bin.d$mu), "e.l"]:bin.d[match(m1, bin.d$mu), "e.u"]
                     m1probtab <- data.frame(m = mseq.all[inds.prob.gz], prob = 1)
                     if(nrow(m1probtab) > npt.inf)
                       m1probtab <- m1probtab[floor(seq(1, nrow(m1probtab), len = npt.inf)), ]
                     m1probtab$prob <- m1probtab$prob / sum(m1probtab$prob)
                     prob_recover_in_one_time_unit <- stats::pexp(q = 1, rate = control_SIR$epi_gamma)
                     hmm.beta <- prob_recover_in_one_time_unit * Rc
                     newinf.seq <- 0:qpois(1 - eps, hmm.beta * m1)
                     newrec.seq <- 0:qbinom(1 - eps, max(m1probtab$m), prob_recover_in_one_time_unit)
                     newinf.mass <- dpois(newinf.seq, hmm.beta * m1)
                     newinf.mass.keep <- newinf.mass[newinf.mass > eps]
                     newinf.seq.keep <- newinf.seq[newinf.mass > eps]
                     newinf.mass.keep <- newinf.mass.keep / sum(newinf.mass.keep)
                     newrec.mass <- dbinom(newrec.seq, m1, prob_recover_in_one_time_unit)
                     newrec.mass.keep <- newrec.mass[newrec.mass > eps]
                     newrec.seq.keep <- newrec.seq[newrec.mass > eps]
                     newrec.mass.keep <- newrec.mass.keep / sum(newrec.mass.keep)
                     new.diff.mat <- outer(newinf.seq.keep, newrec.seq.keep, '-')
                     new.diff.range <- range(new.diff.mat)
                     new.diff.seq <- new.diff.range[1]:new.diff.range[2]
                     newnum.prob <- outer(newinf.mass.keep, newrec.mass.keep, '*')
                     new.diff.dat <- data.frame(new = c(new.diff.mat)[order(c(new.diff.mat))], prob = c(newnum.prob)[order(c(new.diff.mat))])
                     new.diff.dat$cumprob <- cumsum(new.diff.dat$prob)
                     newprobs.cum <- rev(new.diff.dat$cumprob)[match(new.diff.seq, rev(new.diff.dat$new))]
                     newprobs <- diff(c(0, newprobs.cum))
                     new.diff.dat.cum <- data.frame(new = new.diff.seq, prob = newprobs)
                     add.prob.dat <- data.frame(m.new = (min(m1probtab$m) + min(new.diff.dat.cum$new)):(max(m1probtab$m) + max(new.diff.dat.cum$new)), prob = 0)
                     for(j in 1:nrow(m1probtab)){
                       newnumun <- m1probtab$m[j] + new.diff.dat.cum$new
                       newprobs <- new.diff.dat.cum$prob
                       add.prob.dat[match(newnumun, add.prob.dat$m), "prob"] <- add.prob.dat[match(newnumun, add.prob.dat$m), "prob"] + m1probtab$prob[j] * newprobs
                     }
                     if(any(add.prob.dat$m.new < 0)){
                       add.prob.dat$prob[add.prob.dat$m.new == 0] <- sum(add.prob.dat$prob[add.prob.dat$m.new <= 0])
                       add.prob.dat <- add.prob.dat[add.prob.dat$m.new >= 0, ]
                     }
                     if(any(add.prob.dat$m.new > max(mseq.all))){
                       add.prob.dat$prob[add.prob.dat$m.new == max(mseq.all)] <- sum(add.prob.dat$prob[add.prob.dat$m.new >= max(mseq.all)])
                       add.prob.dat <- add.prob.dat[add.prob.dat$m.new <= max(mseq.all), ]
                     }
                     add.prob.dat$bin <- binmap[1 + add.prob.dat$m.new, "bin"]
                     add.prob.dat$cumprob <- cumsum(add.prob.dat$prob)
                     binprobs <- rev(rev(add.prob.dat$cumprob)[match(rev(1:control_debias$nbin), rev(add.prob.dat$bin))])
                     binprobs[is.na(binprobs)] <- 0
                     bin.d.prob.add <- bin.d
                     bin.d.prob.add$prob <- diff(c(0, binprobs))
                     bin.d.prob.add$prob[bin.d.prob.add$prob < 0] <- 0
                     table(binmap[add.prob.dat$m.new, "bin"]) * add.prob.dat$prob[1]
                     trmat[as.character(m1), ] <- bin.d.prob.add$prob
                     if(sum(bin.d.prob.add$prob) > 1 + 1e-10)
                       stop()
                   }
                   saveRDS(trmat, file_out)
                 } else {
                   trmat <- readRDS(file_out)
                 }
                 return(trmat)
               }
stopCluster(cl)

names(trl) <- control_SIR$R_grid
out_file <- paste0("transmats/poisson_SIR_epi_gamma_", control_SIR$epi_gamma, ".RDS")
dir.create("transmats")
saveRDS(trl, out_file)
