## kernel density estimates

library(KernSmooth)
#library(ash)

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

qbox <-
function(x, type="|I", at=0, w=0.5, tick=0.5,
horiz=TRUE, add=FALSE,
lwd=1, lwd_med=2, col=par("fg"), col_med=col, col_box=NA, ...)
{
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
    if (!add) {
        if (horiz)
            plot(q, rep(at, 5), type="n", axes=FALSE, ann=FALSE, ...)
        if (!horiz)
            plot(rep(at, 5), q, type="n", axes=FALSE, ann=FALSE, ...)
    }
    if (a == "|") {
        if (horiz) {
            lines(q[c(1,1)], at+c(tick, -tick)*w, lend=1, lwd=lwd, col=col)
            lines(q[c(5,5)], at+c(tick, -tick)*w, lend=1, lwd=lwd, col=col)
        } else {
            lines(at+c(tick, -tick)*w, q[c(1,1)], lend=1, lwd=lwd, col=col)
            lines(at+c(tick, -tick)*w, q[c(5,5)], lend=1, lwd=lwd, col=col)
        }
    }
    if (a %in% c("-", "|")) {
        if (horiz) {
            lines(q[c(1,2)], c(at, at), lend=1, lwd=lwd, col=col)
            lines(q[c(4,5)], c(at, at), lend=1, lwd=lwd, col=col)
        } else {
            lines(c(at, at), q[c(1,2)], lend=1, lwd=lwd, col=col)
            lines(c(at, at), q[c(4,5)], lend=1, lwd=lwd, col=col)
        }
    }
    if (b %in% c("-", "+")) {
        if (horiz) {
            lines(q[c(2,4)], c(at, at), lend=1, lwd=lwd, col=col)
        } else {
            lines(c(at, at), q[c(2,4)], lend=1, lwd=lwd, col=col)
        }
    }
    if (b  %in% c("=", "I")) {
        if (horiz) {
            polygon(q[c(2,4,4,2)], at+c(1, 1, -1, -1)*w,
                lwd=lwd, col=col_box, border=col)
        } else {
            polygon(at+c(1, 1, -1, -1)*w, q[c(2,4,4,2)],
                lwd=lwd, col=col_box, border=col)
        }
    }
    if (b %in% c("+", "I")) {
        if (horiz) {
            lines(q[c(3,3)], at+c(w, -w), lend=1, lwd=lwd_med, col=col_med)
        } else {
            lines(at+c(w, -w), q[c(3,3)], lend=1, lwd=lwd_med, col=col_med)
        }
    }
    out <- list(
        quantiles=q,
        type=type,
        at=at,
        w=q,
        tick=tick,
        horiz=horiz)
    invisible(out)
}

Type <- c(
    " -", "--", "|-",
    " +", "-+", "|+",
    " =", "-=", "|=",
    " I", "-I", "|I")
set.seed(1)
x <- rnorm(1000)
qbox(x, type=Type[1], at=1, w=0.25, horiz=FALSE,
    xlim=c(0, length(Type)+1), ylim=range(x))
title(main="qbox types", ylab="x")
axis(2)
axis(1, at=1:length(Type), labels=paste0("'", Type, "'"),
    family="mono", lwd=0, cex.axis=0.8)
for (i in 2:length(Type))
    qbox(x, type=Type[i], at=i, w=0.25, horiz=FALSE, add=TRUE)

## d=density, h=hist, s=strip
dplot <-
function(x, type="d", at=0, w=0.5, shift=0.5, col=par("fg"),
border=NA, add=FALSE, horiz=TRUE, ...)
{
    type <- match.arg(type, c("d", "h", "s"))
    ## take the range of z$x instead of x
    if (type == "d") {
        z <- KernSmooth::bkde(x)
        r <- z$y / max(z$y)
    }
    if (type %in% c("s", "h")) {
        h <- hist(x, plot=FALSE)
        z <- stairs(list(x=h$mids, y=h$density))
        r <- h$density / max(h$density)
    }
    if (!add) {
        if (horiz)
            plot(z$x, rep(at, length(z$x)), type="n", axes=FALSE, ann=FALSE, ...)
        if (!horiz)
            plot(rep(at, length(z$x)), z$x, type="n", axes=FALSE, ann=FALSE, ...)
    }
    if (type == "s") {
        z$y <- rep(1, length(z$y))
    } else {
        z$y <- z$y/max(z$y)
    }
    if (shift < 0 || shift > 1)
        stop("shift must be in [0, 1]")
    y1 <- at+2*w*(1-shift)*z$y
    y2 <- at-2*w*shift*z$y
    p <- list(x=c(z$x, rev(z$x)), y=c(y1, rev(y2)))
    if (is.function(col)) {
        Cols <- col(100)
    } else {
        Cols <- if (type == "s")
            colorRampPalette(c("#FFFFFF", col[1L]))(100) else col[1L]
    }
    if (length(Cols) < 2) {
        if (horiz) {
            polygon(p$x, p$y, col=col, border=border)
        } else {
            polygon(p$y, p$x, col=col, border=border)
        }
    } else {
        i <- Cols[cut(r, 100, include.lowest=TRUE, labels=FALSE)]
        if (type == "d") {
            for (j in 2:length(r)) {
                k1 <- c(j-1, j)
                k2 <- c(j, j-1)
                if (horiz) {
                    polygon(z$x[c(k1, k2)], c(y1[k1], y2[k2]),
                        border=i[j], col=i[j])
                } else {
                    polygon(c(y1[k1], y2[k2]), z$x[c(k1, k2)],
                        border=i[j], col=i[j])
                }
            }
        } else {
            for (j in 1:length(r)) {
                k1 <- c(j*2-1, j*2)
                k2 <- c(j*2, j*2-1)
                if (horiz) {
                    polygon(z$x[c(k1, k2)], c(y1[k1], y2[k2]),
                        border=i[j], col=i[j])
                } else {
                    polygon(c(y1[k1], y2[k2]), z$x[c(k1, k2)],
                        border=i[j], col=i[j])
                }
            }
        }
        if (horiz) {
            polygon(p$x, p$y, col=NA, border=border)
        } else {
            polygon(p$y, p$x, col=NA, border=border)
        }
    }
    out <- list(
        polygon=p,
        type=type,
        at=at,
        w=q,
        shift=shift,
        horiz=horiz)
    invisible(out)
}


rug2 <-
function (x, ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"),
    col.ticks = "white", quiet = getOption("warn") < 0, ...)
{
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
    #par("pin")
    if (ticksize < 0.5) {
        tic <- min(diff(u[3:4]), diff(u[1:2])) * ticksize
    } else {
        tic <- ifelse(side%%2 == 1, diff(u[3:4]), diff(u[1:2])) * ticksize
    }
    tic <- ifelse(side%%2 == 1, diff(u[3:4]), diff(u[1:2])) * ticksize
#    if (ticksize < 0)
        opar <- par(xpd = TRUE)
    side <- match.arg(as.character(side), c("1", "2", "3", "4"))
    if (side == "1") {
#        dplot(x, type="s", horiz=TRUE, add=TRUE, at=u[3]+tic/2, w=tic/2)
        dplot(x, type="d", horiz=TRUE, add=TRUE, at=u[3], w=tic, shift=0)
    }
    if (side == "2") {

    }
    if (side == "3") {

    }
    if (side == "4") {

    }
#    switch(as.character(side), `1` = polygon(u[c(1, 2, 2, 1,
#        1)], u[3] + c(0, 0, tic, tic, 0), col = col, border = NA,
#        ...), `2` = polygon(u[1] + c(0, 0, tic, tic, 0), u[c(3,
#        4, 4, 3, 3)], col = col, border = NA, ...), `3` = polygon(u[c(1,
#        2, 2, 1, 1)], u[4] + c(0, 0, -tic, -tic, 0), col = col,
#        border = NA, ...), `4` = polygon(u[2] + c(0, 0, -tic,
#        -tic, 0), u[c(3, 4, 4, 3, 3)], col = col, border = NA,
#        ...))
#    switch(as.character(side), `1` = sapply(x, function(z) lines(c(z,
#        z), u[3] + c(0, tic), col = col.ticks, lwd = lwd)), `2` = sapply(x,
#        function(z) lines(u[1] + c(0, tic), c(z, z), col = col.ticks,
#            lwd = lwd)), `3` = sapply(x, function(z) lines(c(z,
#        z), u[4] + c(0, -tic), col = col.ticks, lwd = lwd)),
#        `4` = sapply(x, function(z) lines(u[2] + c(0, -tic),
#            c(z, z), col = col.ticks, lwd = lwd)))

#    if (ticksize < 0)
        par(opar)
    invisible(x)
}

hist(x)
rug2(x)

## bivariate plot with h/d at margin & bivariate quantiles

quantiles4bkde <- function(x, y) {
    bw <- c(x=bw4bkde(x), y=bw4bkde(y))
    d2 <- bkde2D(cbind(x, y), bandwidth=bw)
    o <- order(d2$fhat)
    i <- order(o)
    f <- d2$fhat
    f <- f/sum(f)
    cs <- cumsum(f[o])[i]
    dim(cs) <- dim(d2$fhat)
    d2$cs <- cs
    d2
}


biplot <- function(x, y, col=1, ...) {
    Pal <- colorRampPalette(c("#FFFFFF", col))(100)
    d <- quantiles4bkde(x, y)

    op <- par(xpd = TRUE, mar=c(5,4,4,4)+0.1)
    plot(x, y, type="n", xlim=range(d$x1), ylim=range(d$x2), axes=FALSE, ...)
    u <- par("usr")
    image(d$x1, d$x2, d$fhat, col=Pal[1:66], add=TRUE)
    #points(x, y, col=paste0(Pal[50], "80"), pch=".")
    contour(d$x1, d$x2, d$cs, add=TRUE, col=Pal[100], levels=c(0.25, 0.5, 0.75))
    box(col=Pal[33])
    axis(1, col=Pal[100])
    axis(2, col=Pal[100])
    dplot(x, horiz=TRUE, add=TRUE, type="d", col=Pal[33], border=NA,
        at=u[4], shift=0, w=0.05*diff(u[3:4]))
    qbox(x, type=" I", horiz=TRUE, col=Pal[66],
        col_med="#FFFFFF", col_box=Pal[66], lwd=2, lwd_med=3,
        at=u[4], w=0.01*diff(u[3:4]), add=TRUE)
    dplot(y, horiz=FALSE, add=TRUE, type="d", col=Pal[33], border=NA,
        at=u[2], shift=0, w=0.05*diff(u[1:2]))
    qbox(y, type=" I", horiz=FALSE, col=Pal[66],
        col_med="#FFFFFF", col_box=Pal[66], lwd=2, lwd_med=3,
        at=u[2], w=0.01*diff(u[1:2]), add=TRUE)
    par(op)

    invisible(NULL)
}

n <- 1000
x <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 1, 1))
y <- 0.8*x+rnorm(n, 0, 3)
biplot(x, y, col=1)


## multivariate boxplot/violin

x1 <- c(rnorm(n/2, -2, 0.5), rnorm(n/2, 1, 1))
x2 <- 0.8*x+rnorm(n, 0, 3)
x3 <- x1 + 3 + 0.5*x2 -2*x1

x <- cbind(A=x1,C=x2,B=x3)

if (!is.list(x)) {
    xx <- lapply(1:NCOL(x), function(i) x[,i])
    names(xx) <- colnames(x)
}
for (i in seq_len(m))
    xx[[i]] <- xx[[i]][!is.na(xx[[i]])]
m <- length(xx)
if (is.null(names(xx)))
    names(xx) <- paste0("V", seq_len(m))

plot(c(1-0.5, m+0.5), range(unlist(xx)), ann=FALSE, axes=FALSE, type="n")
axis(1, seq_len(m), names(xx))
axis(2)
for (i in seq_len(m)) {
    dplot(xx[[i]], horiz=FALSE, at=i, w=0.45, type="d", add=TRUE, col=Pal[50])
    qbox(xx[[i]], type=" I", horiz=FALSE, col=Pal[50],
        col_med=Pal[100], col_box=Pal[1], lwd=2, lwd_med=3,
        at=i, w=0.05, add=TRUE)
}
