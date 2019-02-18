require(MASS)

qst.sim <- function(coal.mat, nn=1){
  cov.mat <- coal.mat
  for(ii in 1:ncol(coal.mat))
    for(jj in 1:ncol(coal.mat))
      cov.mat[ii, jj] <- mean(coal.mat[ii,]) + mean(coal.mat[ii,]) - mean(coal.mat) - coal.mat[ii,jj]
  sq.devs <- mvrnorm(nn, rep(0, ncol(coal.mat)), Sigma=cov.mat)^2
  var.between <- apply(sq.devs, 1, mean)
  return(var.between / (var.between + mean(diag(coal.mat))))
}

cm.to.fst <- function(coal.mat){
  tt <- mean(coal.mat)
  t0 <- mean(diag(coal.mat))
  return((tt - t0)/tt)
}

lk.sim <- function(coal.mat, nd, nn=1){
  fst <- cm.to.fst(coal.mat)
  return(fst * rchisq(n=nn, df=nd-1) / (nd - 1))
}

ring.matrix <- function(N, m, d){
  result <- matrix(data=0, nrow=d, ncol=d)
  for (ii in 1:d) {
    for (jj in 1:d) {
      b <- abs(ii - jj)
      result[ii,jj] <- 2*N*d + (d - b) * b/(2 * m)
    }
  }
  return(result)
}

nospace.matrix <- function(N, m, d){
  result <- matrix(data=0, nrow=d, ncol=d)
  for (ii in 1:d) {
    for (jj in 1:d) {
      if(ii == jj){
        result[ii, jj] <- 2*N*d
      } else{
        result[ii, jj] <- 2*N*d + (d-1)/(2*m)
      }
    }
  }
  return(result)
}

get.m <- function(N, FST, d){
  return((1-FST)*(d^2-1)/(24*FST*N*d))
}

get.m.nospace <- function(N, FST, d){
  return((1-FST)*(d-1)^2/(4*FST*N*d^2))
}

## Simulate popluations with the same FST but differing structure
FST.sim <- 0.1
d.set <- c(2, 4, 6, 10, 20, 50)
m.set <- get.m(1000, FST.sim, d.set)
test.mats <- list()
for(ii in 1:length(d.set)) test.mats[[ii]] <- ring.matrix(1000, m.set[ii], d.set[ii])
qst.norms <- list()
for(ii in 1:length(d.set)) qst.norms[[ii]] <- qst.sim(test.mats[[ii]], nn=1e6)
qst.lks <- list()
for(ii in 1:length(d.set)) qst.lks[[ii]] <- lk.sim(test.mats[[ii]], nd=d.set[ii], nn=2e6)

norm.col <- "#5e3c99"
lk.col <- "#e66101"

ff <- 1/2
pdf("qst_deme_circle.pdf", 12*ff, 15*ff)
par(mfrow=c(3, 2), mai=c(.5, .45, .4, .15), mgp=c(2.3,1,0))
for(ii in 1:length(d.set)){
  yrange <- c(0, max(density(qst.lks[[ii]])$y, 
                     density(qst.norms[[ii]])$y))
  plot(density(qst.lks[[ii]]), type="l", col=lk.col, lty=2, xlim=c(0,0.3), 
       main=paste(d.set[ii], " demes"), lwd=2, xlab="Qst", cex.lab=1.3, ylim=yrange)
  points(density(qst.norms[[ii]]), type="l", lwd=2, col=norm.col)
  points(density(qst.lks[[ii]]), type="l", lwd=2, col=lk.col, lty=2)
  abline(v=mean(qst.norms[[ii]]), lwd=2, col=norm.col)
  abline(v=mean(qst.lks[[ii]]), lwd=2, col=lk.col, lty=2)
  if(ii == 1) legend("topright", legend=c("normal model null", "LK null"), 
                     lwd=1.5, lty=c(1,2), col=c(norm.col, lk.col))
}
dev.off()

norm.quants <- sapply(qst.norms, function(X) quantile(X, .95))
lk.quants <- sapply(qst.lks, function(X) quantile(X, .95))

pdf("qst_deme_percentile.pdf", width=5, height=5)
par(mgp=c(2.3,1,0))
plot(d.set, lk.quants, type="b", log="x", col="#e66101", 
     xlab="number of demes", ylab=expression('95'^'th'*' percentile'), cex.lab=1.3, lwd=2)
points(d.set, norm.quants, type="b", col="#5e3c99", lwd=2)
legend("topright", legend=c("normal model null", "LK null"), col=c("#5e3c99", "#e66101"), 
       lwd=c(2,2))
dev.off()

## Simulate popluations with the same FST but differing structure in population with no spatial component
m.set.nospace <- get.m.nospace(1000, FST.sim, d.set)
test.mats.nospace <- list()
for(ii in 1:length(d.set)) test.mats.nospace[[ii]] <- nospace.matrix(1000, m.set.nospace[ii], d.set[ii])
qst.norms.nospace <- list()
for(ii in 1:length(d.set)) qst.norms.nospace[[ii]] <- qst.sim(test.mats.nospace[[ii]], nn=1e6)
qst.lks.nospace <- list()
for(ii in 1:length(d.set)) qst.lks.nospace[[ii]] <- lk.sim(test.mats.nospace[[ii]], nd=d.set[ii], nn=2e6)

sp.col <- "#0571b0"
nosp.col <- "#92c5de"
lk.col <- "#ca0020"

pdf("qst_deme_compare_foo.pdf", 12*ff, 15*ff)
par(mfrow=c(3, 2), mai=c(.5, .45, .4, .15), mgp=c(2.3,1,0))
for(ii in 1:length(d.set)){
  yrange <- c(0, max(density(qst.lks.nospace[[ii]])$y, 
                     density(qst.norms.nospace[[ii]])$y,
                     density(qst.norms[[ii]])$y))
  plot(density(qst.lks.nospace[[ii]]), type="l", col=lk.col, lty=2, xlim=c(0,0.3), 
       main=paste(d.set[ii], " demes"), lwd=2, xlab=expression("Q"[ST]), cex.lab=1.3, ylim=yrange)
  points(density(qst.norms.nospace[[ii]]), type="l", lwd=2, col=nosp.col)
  points(density(qst.norms[[ii]]), type="l", lwd=2, col=sp.col)
  points(density(qst.lks.nospace[[ii]]), type="l", lwd=2, col=lk.col, lty=2)
  abline(v=mean(qst.norms.nospace[[ii]]), lwd=2, col=nosp.col)
  abline(v=mean(qst.norms[[ii]]), lwd=2, col=sp.col)
  abline(v=mean(qst.lks.nospace[[ii]]), lwd=2, col=lk.col, lty=2)
  if(ii == 1) legend("topright", legend=c("Island model", "Stepping-stone model", "LK null"), 
                     lwd=1.5, lty=c(1, 1, 2), col=c(nosp.col, sp.col, lk.col))
}
dev.off()

norm.quants.nospace <- sapply(qst.norms.nospace, function(X) quantile(X, .95))
lk.quants.nospace <- sapply(qst.lks.nospace, function(X) quantile(X, .95))

pdf("qst_deme_percentile_nospace.pdf", width=5, height=5)
par(mgp=c(2.3,1,0))
plot(d.set, lk.quants.nospace, type="b", log="x", col=lk.col, 
     xlab="number of demes", ylab=expression('95'^'th'*' percentile of Q'[ST]), 
     cex.lab=1.3, lwd=2, lty=2)
points(d.set, norm.quants.nospace, type="b", col=nosp.col, lwd=2)
points(d.set, norm.quants, type="b", col=sp.col, lwd=2)
legend("topright", legend=c("no spatial structure","spatial structure", "LK null"), 
       col=c(nosp.col, sp.col, lk.col), 
       lwd=c(2,2,2), lty=c(1,1,2))
dev.off()
