
source("https://raw.githubusercontent.com/rnorouzian/m/master/m.r")

source("https://raw.githubusercontent.com/izeh/m/master/m.r")

dat <- data.frame(C1 = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
                  C2 = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

lst <- lapply(0:8, function(x) {dat[0:x, ] <- 0; dat })
ESL <- sapply(lst, sum)
EFL <- sum(lengths(dat)) - ESL

inter <- sapply(lst, irr)


plot(1:9, inter[1,], type = "b", ylim = c(min(inter[1,]),  0.7), ylab = "Kappa", panel.f = abline(h = 0, col = 8), xaxt = "n", yaxt = "n",
     xlab = "ESL:EFL Frequency", font.lab = 2, pch = 19, las = 1, mgp = c(1.5, .55, 0), padj = .4, cex.axis = .7, cex.lab = .7)
axis(1, at = 1:9, labels = paste0("(", ESL, ":", EFL, ")"), cex.axis = .7, mgp = c(1, .3, 0))
axis(2, at = c(-.1, 0:6*.1), cex.axis = .7, mgp = c(1.5, .55, 0), las = 1, padj = .4)


plot(1:9, inter[1,], type = "b", ylim = c(min(inter[1,]),  0.7), ylab = "IRR", panel.f = abline(h = 0, col = 8), xaxt = "n", yaxt = "n",
     xlab = "ESL:EFL Frequency", font.lab = 2, pch = 19, las = 1, mgp = c(1.5, .55, 0), padj = .4, cex.axis = .7, cex.lab = .7)
axis(1, at = 1:9, labels = paste0("(", ESL, ":", EFL, ")"), cex.axis = .7, mgp = c(1, .3, 0))

axis(2, at = c(-.1, 0:6*.1), cex.axis = .7, mgp = c(1.5, .55, 0), las = 1, padj = .4)
lines(1:9, inter[2,], type = "b", col = 2, lty = 2, pch = 19)
legend("top", c("Kappa", "S index"), pch = 19, col = 1:2, lty = c(1,2), bty = "n", 
       text.font = 2, cex = .8, horiz = T)





jpeg(file = "S1.jpeg", width = 5.5, height = 3.3, res = 600, units = "in")
par(mgp=c(1.5, .55, 0), tcl = -0.4, mar = c(3.3,3.6,1.1,1.1))
plot(1:9, inter[1,], type = "b", ylim = c(min(inter[1,]),  0.7), ylab = "Kappa", panel.f = abline(h = 0, col = 8), xaxt = "n", yaxt = "n",
     xlab = "ESL:EFL Frequency", font.lab = 2, pch = 19, las = 1, mgp = c(1.5, .55, 0), padj = .4, cex.axis = .7, cex.lab = .7)
axis(1, at = 1:9, labels = paste0("(", ESL, ":", EFL, ")"), cex.axis = .7, mgp = c(1, .3, 0))
axis(2, at = c(-.1, 0:6*.1), cex.axis = .7, mgp = c(1.5, .55, 0), las = 1, padj = .4)

dev.off()





jpeg(file = "S2.jpeg", width = 5.5, height = 3.3, res = 600, units = "in")
par(mgp=c(1.5, .55, 0), tcl = -0.4, mar = c(3.3,3.6,1.1,1.1))
plot(1:9, inter[1,], type = "b", ylim = c(min(inter[1,]),  0.7), ylab = "IRR", panel.f = abline(h = 0, col = 8), xaxt = "n", yaxt = "n",
     xlab = "ESL:EFL Frequency", font.lab = 2, pch = 19, las = 1, mgp = c(1.5, .55, 0), padj = .4, cex.axis = .7, cex.lab = .7)
axis(1, at = 1:9, labels = paste0("(", ESL, ":", EFL, ")"), cex.axis = .7, mgp = c(1, .3, 0))

axis(2, at = c(-.1, 0:6*.1), cex.axis = .7, mgp = c(1.5, .55, 0), las = 1, padj = .4)
lines(1:9, inter[2,], type = "b", col = 2, lty = 2, pch = 19)
legend("top", c("Kappa", "S index"), pch = 19, col = 1:2, lty = c(1,2), bty = "n", 
       text.font = 2, cex = .8, horiz = T)
dev.off()





irr <- int <- function(X, nsim = 1e3, useNA = "ifany", level = .95, digits = 6, raw = TRUE, cats = NULL) 
{
  
  if(!inherits(X, c("data.frame", "matrix", "table"))) stop("Ratings must be 'data.frame', 'matrix', and if not raw, a 'table'.", call. = FALSE)
 
   cats <- if(is.null(cats)) levels(factor(unlist(X))) else cats
  
  if(raw) X <- table(row(X), factor(unlist(X), levels = cats), useNA = useNA) 
  
  X2 <- X * (X - 1)
  sumcol <- colSums(X)
  sumrow <- rowSums(X)
  nc <- ncol(X)
  nr <- nrow(X)
  tot <- sum(X)
  pij <- X2/(sumrow * (sumrow - 1))
  pi <- rowSums(pij)
  p <- mean(pi)
  pj <- sumcol/tot
  pj2 <- pj^2
  pe <- sum(pj2)
  KAPPA <- (p - pe)/(1 - pe)
  s <- (nc * p - 1)/(nc - 1)
  pi.v.boot <- replicate(nsim, pi.boot <- sample(pi, size = nr, replace = TRUE))
  p.boot <- apply(pi.v.boot, 2, mean)
  s.boot <- sapply(seq_len(nsim), function(i) (nc * p.boot[i] - 1)/(nc - 1))
  
  p <- (1 - level) / 2
  s.boot.ci <- quantile(s.boot, probs = c(p, 1-p), na.rm = TRUE)
  
  return(round(c(KAPPA = KAPPA, 
                 S.index = s, 
                 lower = s.boot.ci[[1]], 
                 upper = s.boot.ci[[2]], 
                 conf.level = level), digits))
}                             
