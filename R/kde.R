## kernel density estimates

library(KernSmooth)

bw4bkde <-
function (x, kernel = "normal", canonical = FALSE, bandwidth,
    gridsize = 401L, range.x, truncate = TRUE) {
    if (!missing(bandwidth) && bandwidth <= 0)
        stop("'bandwidth' must be strictly positive")
    kernel <- match.arg(kernel, c("normal", "box", "epanech",
        "biweight", "triweight"))
    n <- length(x)
    del0 <- switch(kernel, normal = (1/(4 * pi))^(1/10), box = (9/2)^(1/5),
        epanech = 15^(1/5), biweight = 35^(1/5), triweight = (9450/143)^(1/5))
    h <- if (missing(bandwidth))
        del0 * (243/(35 * n))^(1/5) * sqrt(var(x))
    else if (canonical)
        del0 * bandwidth
    else bandwidth
    h
}

n <- 1000
x <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 1, 1))
y <- 0.8*x+rnorm(n, 0, 3)
## 1D
bw <- c(x=bw4bkde(x), y=bw4bkde(y))
dx <- bkde(x, bandwidth=bw["x"])
dy <- bkde(y, bandwidth=bw["y"])
d2 <- bkde2D(cbind(x, y), bandwidth=bw)
qx <- quantile(x, c(0, 0.25, 0.5, 0.75, 1))
qy <- quantile(y, c(0, 0.25, 0.5, 0.75, 1))
#Cvol <- (d2$x1[2]-d2$x1[1]) * (d2$x2[2]-d2$x2[1])
#Vsum <- sum(d2$fhat * Cvol)
o <- order(d2$fhat)
i <- order(o)
f <- d2$fhat
f <- f/sum(f)
cs <- cumsum(f[o])[i]
dim(cs) <- dim(d2$fhat)

ColA <- "#1b9e77"
ColB <- "#7570b3"
ColC <- "#000000"
Col1 <- paste0(ColA, "80")
Col2 <- paste0(ColA, "40")
Col3 <- paste0(ColB, "80")
Col4 <- paste0(ColB, "40")
Col5 <- paste0(ColC, "20")
Pal <- colorRampPalette(c("#FFFFFF", Col2))
S <- 4
mat <- matrix(1, S, S)
mat[1,] <- 2
mat[,S] <- 3
mat[1,S] <- 4
layout(mat)
op <- par(mar=c(6, 5, 0.5, 0.5)+0.1)
plot(x, y, axes=FALSE, type="n", xlim=range(dx$x), ylim=range(dy$x))
image(d2$x1, d2$x2, d2$fhat, col=Pal(100), add=TRUE)
points(x, y, col=Col5, pch=19)
contour(d2$x1, d2$x2, cs, add=TRUE, col=ColB, levels=c(0.25, 0.5, 0.75))
box(col=Col2)
axis(1)
axis(2)
par(mar=c(0, 5, 2, 0.5)+0.1)
yl <- range(dx$y)
#yl[1] <- yl[1]-0.1*diff(yl)
plot(dx$x, dx$y, type="n", axes=FALSE, ann=FALSE, ylim=yl)
polygon(dx$x, dx$y, col=Col2, border=Col1)
lines(qx[c(2, 4)], c(0, 0), col=ColB, lwd=7, lend=2, xpd=TRUE)
lines(qx[c(1, 5)], c(0, 0), col=ColB, lwd=2, lend=2, xpd=TRUE)
points(qx[3], 0, col=ColB, pch=3, xpd=TRUE, cex=2)
par(mar=c(6, 0.5, 0.5, 2)+0.1)
xl <- range(dx$y)
#xl[1] <- xl[1]-0.1*diff(xl)
plot(dy$y, dy$x, type="n", axes=FALSE, ann=FALSE, xlim=xl)
polygon(dy$y, dy$x, col=Col2, border=Col1)
lines(c(0, 0), qy[c(2, 4)], col=ColB, lwd=7, lend=2, xpd=TRUE)
lines(c(0, 0), qy[c(1, 5)], col=ColB, lwd=2, lend=2, xpd=TRUE)
points(0, qy[3], col=ColB, pch=3, xpd=TRUE, cex=2)
par(op)