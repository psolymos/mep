## kernel density estimates

library(KernSmooth)
library(ash)

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
op <- par(mar=c(6, 5, 0, 0)+0.1)
plot(x, y, axes=FALSE, type="n", xlim=range(dx$x), ylim=range(dy$x))
image(d2$x1, d2$x2, d2$fhat, col=Pal(100), add=TRUE)
points(x, y, col=Col5, pch=19)
contour(d2$x1, d2$x2, cs, add=TRUE, col=ColB, levels=c(0.25, 0.5, 0.75))
#box(col=Col2)
axis(1)
axis(2)
par(mar=c(0, 5, 2, 0)+0.1)
yl <- range(dx$y)
#yl[1] <- yl[1]-0.1*diff(yl)
plot(dx$x, dx$y, type="n", axes=FALSE, ann=FALSE, ylim=yl)
polygon(dx$x, dx$y, col=Col2, border=Col1)
lines(qx[c(2, 4)], c(0, 0), col=ColB, lwd=7, lend=2, xpd=TRUE)
lines(qx[c(1, 5)], c(0, 0), col=ColB, lwd=2, lend=2, xpd=TRUE)
points(qx[3], 0, col=ColB, pch=3, xpd=TRUE, cex=2)
par(mar=c(6, 0, 0, 2)+0.1)
xl <- range(dx$y)
#xl[1] <- xl[1]-0.1*diff(xl)
plot(dy$y, dy$x, type="n", axes=FALSE, ann=FALSE, xlim=xl)
polygon(dy$y, dy$x, col=Col2, border=Col1)
lines(c(0, 0), qy[c(2, 4)], col=ColB, lwd=7, lend=2, xpd=TRUE)
lines(c(0, 0), qy[c(1, 5)], col=ColB, lwd=2, lend=2, xpd=TRUE)
points(0, qy[3], col=ColB, pch=3, xpd=TRUE, cex=2)
par(op)

rugden <-
function (x, ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"),
    col.ticks = "white", quiet = getOption("warn") < 0, ...)
{
    Pal <- colorRampPalette(c("#FFFFFF", col))(100)
    dx <- bkde(x)
    Col <- cut(dx$y, 100, include.lowest=TRUE, labels=FALSE)
    #plot(dx$x, dx$y, col=Pal[Col])
    xi <- dx$x
    st <- 0.5*(xi[2]-xi[1])

    x <- as.vector(x)
    ok = is.finite(x)
    x <- x[ok]
    if (!quiet) {
        u <- par("usr")
        if (side%%2 == 1) {
            if (par("xlog"))
                u <- 10^u[1:2]
            else u <- u[1:2]
        }
        else {
            if (par("ylog"))
                u <- 10^u[3:4]
            else u <- u[3:4]
        }
        if (any(x < u[1] | x > u[2]))
            warning("Some values will be clipped")
    }
    u <- par("usr")
    par("pin")
    if (ticksize < 0.5)
        tic <- min(diff(u[3:4]), diff(u[1:2])) * ticksize
    else tic <- ifelse(side%%2 == 1, diff(u[3:4]), diff(u[1:2])) *
        ticksize
    if (ticksize < 0)
        opar <- par(xpd = TRUE)
    switch(as.character(side),
        `1` = polygon(u[c(1, 2, 2, 1,
            1)], u[3] + c(0, 0, tic, tic, 0), col = col, border = NA, ...),
        `2` = polygon(u[1] + c(0, 0, tic, tic, 0), u[c(3,
            4, 4, 3, 3)], col = col, border = NA, ...),
        `3` = polygon(u[c(1,
            2, 2, 1, 1)], u[4] + c(0, 0, -tic, -tic, 0), col = col,
            border = NA, ...),
        `4` = polygon(u[2] + c(0, 0, -tic,
            -tic, 0), u[c(3, 4, 4, 3, 3)], col = col, border = NA, ...))
    if (ticksize < 0)
        par(opar)
    invisible(x)
}

## add rug related to: 1D density
## add rug related to: boxplot style

## join the dots
stairs <- function(x, y=NULL, ...) {
    if (is.null(y)) {
        if (is.list(x)) {
            y <- x[[2L]]
            x <- x[[1L]]
        } else {
            if (!is.null(dim(x))) {
                y <- x[,2L]
                x <- x[,1L]
            } else {
                stop("y must be specified")
            }
        }
    }
    n <- length(x)
    o <- order(x)
    x <- x[o]
    y <- y[o]
    d <- 0.5*diff(x)
    list(
        x=rep(c(x[1]-d[1], x[-1]-d, x[n]+d[n-1]), c(1, rep(2, n-1), 1)),
        y=rep(y, each=2))
}
set.seed(1)
x <- cumsum(runif(10, 0.5, 1))
y <- runif(10, 0, 1)
plot(x,y)
lines(x, y)
lines(spline(x, y, n = 10 * length(x)), col=2)
lines(stairs(x, y), col=4)

## boxplot

x <- rnorm(1000)
qbox <- function(x, type="|I", at=0, w=0.5, tick=0.5, horiz=TRUE, add=FALSE, ...) {
    q <- quantile(x, seq(0, 1, 0.25))
    types <- c(
        " -", "--", "|-",
        " +", "-+", "|+",
        " =", "-=", "|=",
        " I", "-I", "|I")
    type <- match.arg(type, types)
    if (nchar(type) < 2)
        stop("type must be fully specified")
    a <- substr(type, 1, 1)
    b <- substr(type, 2, 2)
    if (!add)
        plot(NA, xlim=range(q), ylim=at+c(1, -1)*w, axes=FALSE, ann=FALSE)#, ...)
    if (a == "|") {
        lines(q[c(1,1)], at+c(tick, -tick)*w)
        lines(q[c(5,5)], at+c(tick, -tick)*w)
    }
    if (a %in% c("-", "|")) {
        lines(q[c(1,2)], c(at, at))
        lines(q[c(4,5)], c(at, at))
    }
    if (b %in% c("-", "+")) {
        lines(q[c(2,4)], c(at, at))
    }
    if (b  %in% c("=", "I")) {
        polygon(q[c(2,4,4,2)], at+c(1, 1, -1, -1)*w)
    }
    if (b %in% c("+", "I")) {
        lines(q[c(3,3)], at+c(w, -w))
    }
    invisible(NULL)
}

## d=density, h=hist, s=strip
dplot <- function(x, type="d", at=0, w=0.5, alpha=0, col=1, border=NA, add=FALSE) {
    type <- match.arg(type, c("d", "h", "s"))
    if (!add)
        plot(NA, xlim=range(x), ylim=at+c(1, -1)*w, axes=FALSE, ann=FALSE)#, ...)
    if (type == "d") {
        z <- KernSmooth::bkde(x)
        r <- z$y / max(z$y)
    }
    if (type %in% c("s", "h")) {
        h <- hist(x, plot=FALSE)
        z <- stairs(list(x=h$mids, y=h$density))
        r <- h$density / max(h$density)
    }
    if (type == "s") {
        z$y <- rep(1, length(z$y))
    } else {
        z$y <- z$y/max(z$y)
    }
    y1 <- at+w*(1-alpha)*z$y
    y2 <- at-w*(1-alpha)*z$y
    p <- list(x=c(z$x, rev(z$x)), y=c(y1, rev(y2)))
    ## need to distinguish full color from gradient, for type=s
    if (is.function(col)) {
        Cols <- col(100)
    } else {
        if (length(col) > 1L) {
            Cols <- colorRampPalette(col)(100)
        } else {
            Cols <- colorRampPalette(c("#FFFFFF", col))(100)
        }
    }
    i <- Cols[cut(r, 100, include.lowest=TRUE, labels=FALSE)]
    if (type == "d") {
        for (j in 2:length(r)) {
            k1 <- c(j-1, j)
            k2 <- c(j, j-1)
            polygon(z$x[c(k1, k2)], c(y1[k1], y2[k2]), border=NA, col=i[j])
        }
    } else {
        for (j in 1:length(r)) {
            k1 <- c(j*2-1, j*2)
            k2 <- c(j*2, j*2-1)
            polygon(z$x[c(k1, k2)], c(y1[k1], y2[k2]), border=NA, col=i[j])
        }
    }
    ## this does not work for density -- check that, can't use stairs

    polygon(p, col=NA, border=border)
    invisible(NULL)
}
