## TODO:
## - conditional plots for fx2
## - heatmap with +/- corr
## - pairs plot

## Exploration

find_bw <-
function (x, kernel = "normal", canonical = FALSE, bandwidth,
    gridsize = 401L, range.x, truncate = TRUE)
{
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
density1 <-
function(x, bw=NULL, ...)
{
    bw <- if (is.null(bw))
        find_bw(x, ...) else find_bw(x, bandwidth=bw, ...)
    out <- KernSmooth::bkde(x, bandwidth=bw, ...)
    names(out)[names(out) == "y"] <- "f"
    out$bw <- bw
    out
}

cramer_v <-
function(x)
{
    n <- sum(x)
    E <- outer(rowSums(x), colSums(x), "*") / n
    ChiSq <- sum((x - E)^2 / E)
    sqrt((ChiSq / n) / min(nrow(x) - 1L, ncol(x) - 1L))
}

sum_fx <- function(x) {
    f <- x$f
    if (is.null(dim(f))) {
        d <- diff(as.numeric(x$x))
        d <- c(d[1L], d)
        out <- sum(f * d)
    } else {
        d1 <- diff(as.numeric(x$x1))
        d1 <- c(d1[1L], d1)
        d2 <- diff(as.numeric(x$x2))
        d2 <- c(d2[1L], d2)
        out <- sum(t(f * d1) * d2)
    }
    out
}

get_type <- function(x) {
    if (is.logical(x))
        return("int")
    if (is.character(x))
        return("nom")
    type <- if (!is.factor(x))
        "cnt" else if (is.ordered(x))
            "ord" else "nom"
    if (type == "cnt")
        if (isTRUE(all.equal(as.vector(x), as.integer(round(x + 0.001)))))
            type <- "int"
    type
}

fx1 <-
function(x, ...)
{
    if (!is.null(dim(x)))
        stop("x must be a vector")
    out <- list(
        x=NULL,
        f=NULL,
        n=NULL,
        bw=NULL,
        name=deparse(substitute(x)),
        range=NULL,
        type=NULL)
    x <- x[!is.na(x)]
    out$n <- length(x)
    if (is.character(x))
        x <- as.factor(x)
    if (is.logical(x))
        x <- as.integer(x)
    out$type <- get_type(x)
    if (out$type == "cnt") {
        z <- density1(x, ...)
        out$x <- z$x
        out$f <- z$f
        out$bw <- z$bw
        out$range <- range(x)
    }
    if (out$type == "int") {
        z <- table(x) / out$n
        out$x <- sort(unique(x))
        out$f <- as.numeric(z[match(as.character(out$x), names(z))])
        out$range <- range(x)
    }
    if (out$type %in% c("ord", "nom")) {
        z <- table(x) / out$n
        out$x <- factor(levels(x), levels(x), ordered=out$type == "ord")
        out$f <- as.numeric(z[match(levels(x), names(z))])
        out$range <- range(as.numeric(out$x))
    }
    class(out) <- "fx1"
    out
}

axis_range <- function(range, type, w=0.5, dim=1) {
    lim <- range
    if (dim == 1) {
        if (type == "int") {
            lim[1] <- lim[1] - 0.05*diff(range)
            lim[2] <- lim[2] + 0.05*diff(range)
        }
        if (type %in% c("ord", "nom")) {
            lim[1] <- lim[1] - w*1.05
            lim[2] <- lim[2] + w*1.05
        }
    } else {
        if (type != "cnt")
            lim[1] <- lim[1] - 0.5
            lim[2] <- lim[2] + 0.5
    }
    lim
}

plot.fx1 <-
function(x, col=par("col"), las=1, ...)
{
    #col <- colorRampPalette(c("#FFFFFF", col))(100)[66]
    w <- if (x$type == "ord")
        0.4 else 0.5
    xlim <- axis_range(x$range, x$type, w=w, dim=1)
    plot(as.numeric(x$x), x$f, xaxs = "i", yaxs = "i",
        type="n", ann=FALSE, axes=FALSE, xlim=xlim, ...)
    if (x$type == "cnt") {
#        i <- x$x >= x$range[1] & x$x <= x$range[2]
#        xx <- x$x[i]
#        xx <- c(xx[1], xx, xx[length(xx)])
#        xf <- c(0, x$f[i], 0)
#        polygon(xx, xf, col=col, border=col)
        polygon(x$x, x$f, col=col, border=col)
        a1 <- axis(1, tick=FALSE, las=las)
        rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
    }
    if (x$type == "int") {
        segments(x0=x$x, y0=rep(0, length(x$f)), y1=x$f, lend=1,
            col=col, lwd=3)
        if (length(x$x) <= 5) {
            a1 <- axis(1, x$x, x$x, tick=FALSE, las=las)
            rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
        } else {
            a1 <- axis(1, tick=FALSE, las=las) # plot(table()) has values only at x$x
            rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
        }
    }
    if (x$type %in% c("ord", "nom")) {
        mid <- as.numeric(x$x)
        for (i in seq_len(length(x$x))) {
            if (x$f[i] > 0)
                polygon(mid[i]+c(-w, -w, w, w),
                    c(0, x$f[i], x$f[i], 0), border=col, col=col)
        }
        a1 <- axis(1, mid, levels(x$x), tick=FALSE, las=las)
        rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
        if (x$type == "ord")
            axis(1, mid[-1]-0.5, rep("<", nlevels(x$x)-1), tick=FALSE, las=0)
    }
        a2 <- axis(2, tick=FALSE, las=las)
        rug(a2, -0.05, side=2, lwd=1, quiet=TRUE)
        title(ylab=paste0("f(", x$name, ")"), xlab=x$name)
    invisible(x)
}

if (FALSE) {

## continuous
x1 <- rnorm(100)
xx <- fx1(x1)
sum_fx(xx)
## character = nominal
x2 <- sample(LETTERS[1:4], 100, replace=TRUE, prob=4:1)
xx <- fx1(x2)
sum_fx(xx)
## factor = nominal
xx <- fx1(as.factor(x2))
sum_fx(xx)
## ordered = ordinal
xx <- fx1(as.ordered(x2))
sum_fx(xx)
## count data = integer
x3 <- rpois(100, 2) - 1
xx <- fx1(x3)
sum_fx(xx)

col <- 2
op <- par(mfrow=c(2,2))
plot(fx1(x1), col=col)
plot(fx1(x2), col=col)
plot(fx1(as.ordered(x2)), col=col)
plot(fx1(x3), col=col)
par(op)

x1 <- rnorm(100)
x2 <- as.factor(sample(LETTERS[1:4], 100, replace=TRUE, prob=4:1))
x1 <- as.factor(sample(LETTERS[5:10], 100, replace=TRUE))
str(density2(x1, x1))
str(density2(x1, x2))
str(density2(x2, x1))
str(density2(x2, x2))

}

fx2 <-
function(x1, x2, bw1=NULL, bw2=NULL, m1=51L, m2=51L, ...)
{
    if (!is.null(dim(x1)))
        stop("x1 must be a vector")
    if (!is.null(dim(x2)))
        stop("x2 must be a vector")
    out <- list(
        x1=NULL,
        x2=NULL,
        f=NULL,
        n=NULL,
        bw1=NULL,
        bw2=NULL,
        cor=NULL,
        range1=NULL,
        range2=NULL,
        name1=deparse(substitute(x1)),
        name2=deparse(substitute(x2)),
        type1=NULL,
        type2=NULL)
    keep <- !is.na(x1) & !is.na(x2)
    x1 <- x1[keep]
    x2 <- x2[keep]
    out$n <- sum(keep)
    if (is.character(x1))
        x1 <- as.factor(x1)
    if (is.character(x2))
        x2 <- as.factor(x2)
    if (is.logical(x1))
        x1 <- as.integer(x1)
    if (is.logical(x2))
        x2 <- as.integer(x2)
    out$type1 <- get_type(x1)
    out$type2 <- get_type(x2)
    out$range1 <- range(as.numeric(x1))
    out$range2 <- range(as.numeric(x2))

    if (out$type1 != "cnt")
        bw1 <- 0.25
    if (out$type2 != "cnt")
        bw2 <- 0.25
    bw1 <- if (is.null(bw1))
        find_bw(x1) else find_bw(x1, bandwidth=bw1)
    bw2 <- if (is.null(bw2))
        find_bw(x2) else find_bw(x2, bandwidth=bw2)
    z <- KernSmooth::bkde2D(
        cbind(as.numeric(x1), as.numeric(x2)),
        bandwidth=c(bw1, bw2),
        gridsize=c(m1, m2), ...)
    out$x1 <- z$x1
    out$x2 <- z$x2
    out$f <- z$fhat

    if (out$type1 == "cnt") {
        i1 <- seq_len(length(out$x1))
        d1 <- diff(z$x1[1:2])
    } else {
        d1 <- 1
    }
    if (out$type1 == "int") {
        x1n <- seq(min(x1), max(x1), by=1)
        i1 <- sapply(x1n, function(z) which.min(abs(z - out$x1)))
    }
    if (out$type1 %in% c("ord", "nom")) {
        x1n <- factor(levels(x1), levels(x1), ordered=out$type1 == "ord")
        x1c <- as.numeric(x1n)
        i1 <- sapply(x1c, function(z) which.min(abs(z - out$x1)))
    }
    if (out$type2 == "cnt") {
        i2 <- seq_len(length(out$x2))
        d2 <- diff(z$x2[1:2])
    } else {
        d2 <- 1
    }
    if (out$type2 == "int") {
        x2n <- seq(min(x2), max(x2), by=1)
        i2 <- sapply(x2n, function(z) which.min(abs(z - out$x2)))
    }
    if (out$type2 %in% c("ord", "nom")) {
        x2n <- factor(levels(x2), levels(x2), ordered=out$type2 == "ord")
        x2c <- as.numeric(x2n)
        i2 <- sapply(x2c, function(z) which.min(abs(z - out$x2)))
    }
    if (out$type1 != "cnt")
        out$x1 <- x1n
    if (out$type2 != "cnt")
        out$x2 <- x2n
    ## subset & standardize
    out$f <- out$f[i1, i2, drop=FALSE]
    out$f <- out$f / (d1 * d2 * sum(out$f))
    out$bw1 <- bw1
    out$bw2 <- bw2

    if (out$type1 != "nom" && out$type2 != "nom")
        out$cor <- cor(as.numeric(x1), as.numeric(x2), method="spearman")
    if (out$type1 != "nom" && out$type2 == "nom") {
        LM <- lm(as.numeric(x1) ~ x2)
        out$cor <- sqrt(summary(LM)$r.squared)
    }
    if (out$type1 == "nom" && out$type2 != "nom") {
        LM <- lm(as.numeric(x2) ~ x1)
        out$cor <- sqrt(summary(LM)$r.squared)
    }
    if (out$type1 == "nom" && out$type2 == "nom")
        out$cor <- cramer_v(out$f * out$n)
    out$condition <- 0L
    class(out) <- "fx2"
    out
}

flip_fx2 <- function(x)
{
    tmp <- x
    x$x1 <- tmp$x2
    x$x2 <- tmp$x1
    x$f <- t(x$f)
    x$bw1 <- tmp$bw2
    x$bw2 <- tmp$bw1
    x$range1 <- tmp$range2
    x$range2 <- tmp$range1
    x$type1 <- tmp$type2
    x$type2 <- tmp$type1
    x$name1 <- tmp$name2
    x$name2 <- tmp$name1
    if (x$condition > 0)
        x$condition <- switch(as.character(x$condition),
            "1" = 2L,
            "2" = 1L)
    x
}

plot.fx2 <-
function(x, flip=FALSE, col=par("col"), las=1, ...)
{
    if (flip)
        x <- flip_fx2(x)
    if (!is.function(col) && length(col) < 2L) {
        col <- colorRampPalette(c(par("bg"), col))(100)
    } else {
        if (is.function(col))
            col <- col(100)
    }
    xlim <- axis_range(x$range1, x$type1, dim=2)
    ylim <- axis_range(x$range2, x$type2, dim=2)

    image(as.numeric(x$x1), as.numeric(x$x2), x$f, col=col,
        ann=FALSE, axes=FALSE, xlim=xlim, ylim=ylim, ...)
    if (x$condition > 0) {
        if (x$condition == 1L) {
            xlab <- paste("|", x$name1)
            ylab <- x$name2
        } else {
            xlab <- x$name1
            ylab <- paste("|", x$name2)
        }
    } else {
        xlab <- x$name1
        ylab <- x$name2
    }
    title(xlab=xlab, ylab=ylab)
    if (x$type1 %in% c("cnt", "int")) {
        if (x$type1 == "int" && length(x$x1) <= 5) {
            a1 <- axis(1, x$x1, x$x1, tick=FALSE, las=las)
            rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
        } else {
            a1 <- axis(1, tick=FALSE, las=las) # plot(table()) has values only at x$x
            rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
        }
    } else {
        a1 <- axis(1, as.numeric(x$x1), levels(x$x1), tick=FALSE, las=las)
        rug(a1, -0.05, side=1, lwd=1, quiet=TRUE)
        if (x$type1 == "ord")
            axis(1, as.numeric(x$x1)[-1]-0.5,
                rep("<", nlevels(x$x1)-1), tick=FALSE, las=0)
    }
    if (x$type2 %in% c("cnt", "int")) {
        if (x$type2 == "int" && length(x$x2) <= 5) {
            a2 <- axis(2, x$x2, x$x2, tick=FALSE, las=las)
            rug(a2, -0.05, side=1, lwd=1, quiet=TRUE)
        } else {
            a2 <- axis(2, tick=FALSE, las=las)
            rug(a2, -0.05, side=2, lwd=1, quiet=TRUE)
        }
    } else {
        a2 <- axis(2, as.numeric(x$x2), levels(x$x2), tick=FALSE, las=las)
        rug(a2, -0.05, side=2, lwd=1, quiet=TRUE)
        if (x$type2 == "ord")
            axis(2, as.numeric(x$x2)[-1]-0.5,
                rep("<", nlevels(x$x2)-1), tick=FALSE, las=0)
    }
    invisible(x)
}

if (FALSE) {
z <- fx2(rnorm(100),rnorm(100))
sum_fx(z)
z <- fx2(rnorm(100),rpois(100, 2) - 1)
sum_fx(z)
x2 <- sample(LETTERS[1:4], 100, replace=TRUE, prob=4:1)
x1 <- sample(LETTERS[1:4], 100, replace=TRUE, prob=4:1)
z <- fx2(x1,x2)
sum_fx(z)
z <- fx2(rnorm(100),x2)
sum_fx(z)
z <- fx2(rnorm(100),as.ordered(x2))
sum_fx(z)
image(as.numeric(z[[1]]), as.numeric(z[[2]]), z[[3]])
}

cnd_fx2 <-
function(x, condition, breaks=4, ...)
{
    if (x$condition > 0)
        stop("already conditioned")
    if (missing(condition))
        stop("specifiy condition")
    if (length(condition) != 1L)
        stop("specify one condition")
    if (!is.numeric(condition))
        condition <- match(as.character(condition),
            unlist(x[c("name1", "name2")]))
    if (!(condition %in% 1:2))
        stop("condition must be in 1:2")
    if (condition == 1) # only condition on x2
        x <- flip_fx2(x)

    if (x$type2 == "cnt") {
        ct <- cut(x$x2, breaks)
        tb <- as.numeric(table(ct))
        mm <- model.matrix(~ct-1)
        fun <- function(v) drop(t(mm) %*% v) / tb
        #x$x2 <- unname(fun(x$x2))
        x$x2 <- factor(levels(ct), levels(ct), ordered=TRUE)
        x$f <- unname(t(apply(x$f, 1, fun)))
        x$range2 <- range(as.numeric(x$x2))
        x$type2 <- "ord" # cell size is now irrelevant
    }
    for (i in seq_along(x$x2))
        x$f[,i] <- x$f[,i] / sum(x$f[,i])
    x$condition <- 2L

    if (condition == 1)
        flip_fx2(x) else x
}

## slice 'n' dice ~ snd
snd <-
function(x, ...)
{
    fx1 <- lapply(colnames(x), function(i) fx1(x[,i]))
    names(fx1) <- colnames(x)
    for (i in seq_len(length(fx1)))
        out[[i]]$name <- colnames(x)[i]
    K <- ncol(x)
    D <- diag(1, K, K)
    ID <- data.frame(
        row = row(D)[lower.tri(D)],
        col = col(D)[lower.tri(D)])
    fx2 <- lapply(rownames(ID), function(i) fx2(x[,ID[i,1]], x[,ID[i,2]]))
    for (i in seq_len(length(fx2))) {
        fx2[[i]]$name1 <- colnames(x)[ID[i,1]]
        fx2[[i]]$name2 <- colnames(x)[ID[i,2]]
    }
    D[lower.tri(D)] <- sapply(fx2, "[[", "cor")
    D <- t(D)
    D[lower.tri(D)] <- sapply(fx2, "[[", "cor")
    dimnames(D) <- list(colnames(x), colnames(x))
    out <- list(
        fx1=fx1,
        fx2=fx2,
        id = data.frame(
            row=factor(colnames(x)[ID$row], colnames(x)),
            col=factor(colnames(x)[ID$col], colnames(x))),
        cor=D)
    class(out) <- "snd"
    out
}

## --- OK so far ---

n <- 100
dat <- data.frame(
    cnt1 = rnorm(n),
    cnt2 = rgamma(n, 2),
    int1 = rpois(n, 2) - 1,
    int2 = rbinom(n, 1, 1/3),
    nom1 = factor(sample(LETTERS[1:3], n, replace=TRUE, prob=3:1), LETTERS[3:1]),
    nom2 = as.factor(sample(letters[1:10], n, replace=TRUE, prob=sqrt(1:10))))
dat$ord1 <- as.ordered(sample(dat$nom1))
dat$ord2 <- as.ordered(sample(dat$nom2))
summary(dat)

x <- fx2(dat$cnt1, dat$cnt2)
str(x)
str(cnd_fx2(x, 1))
str(cnd_fx2(x, 2))
plot(x)
plot(cnd_fx2(x, 2))
plot(cnd_fx2(x, 1))


pairs.default(dat)
pairs.default(dat,
    panel=function(x,y,...) plot(fx2(x,y,...)),
    diag.panel=function(x, ...) plot(fx1(x,...)))



plot_slice <-
function(x, name=NULL, ask, ...)
{
    x <- x$fx1
    vars <- names(x)
    if (is.null(name))
        name <- vars
    name <- if (is.character(name)) {
        name[match(vars, name)]
    } else {
        vars[name]
    }
    name <- name[!is.na(name)]
    vars <- vars[name]
    np <- length(vars)
    if (np < 1)
        stop("must define at least one variable")
    if (missing(ask))
        ask <- prod(par("mfcol")) < np && dev.interactive()
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for (i in seq_len(np)) {
        plot(x[[i]], ...)
    }
    invisible(x)
}

x <- slice(dat)
#plot(x)

plot_dice <-
function(x, name1=NULL, name2=NULL, ask, ...)
{
        x <- x$fx2
    if (is.null(name1) && is.null(name2)) {
        if (missing(ask))
            ask <- prod(par("mfcol")) < nrow(x$id) && dev.interactive()
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        for (i in seq_len(length(x$fx))) {
            plot(x$fx[[i]], ...)
        }
    } else {
        vars <- colnames(x$cor)
        ID <- lapply(x$id, as.numeric)
        if (length(name1) != 1L)
            stop("name1 must be of length 1")
        if (length(name2) != 1L)
            stop("name2 must be of length 1")
        if (!is.numeric(name1))
            name1 <- which(vars == as.character(name1))
        if (!is.numeric(name2))
            name2 <- which(vars == as.character(name2))
        if (name1 == name2)
            stop("use slice, not dice")
        if (name1 > name2) {
            plot(x$fx[[which(ID$row == name1 & ID$col == name2)]], ...)
        } else {
            plot(flip_fx2(x$fx[[which(ID$row == name2 & ID$col == name1)]]), ...)
        }
    }
    invisible(x)
}

x <- dice(dat)
plot(x, 1, 2)
plot(x, "int1", "nom2")
plot(x, "nom2", "int1")

h <- hclust(dist(x$cor))
image(x$cor[h$order, h$order][,K:1])


i=1;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=2;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=3;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=4;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=5;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=6;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=7;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=8;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=9;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=10;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=11;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=12;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=13;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=14;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=15;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=16;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=17;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=18;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=19;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=20;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=21;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=22;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=23;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))
i=24;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]])) # !
i=25;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]])) # !
i=26;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]])) # !
i=27;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]])) # !
i=28;str(fx2(dat[,ID[i,1]], dat[,ID[i,2]]))

## need to separate plotting region, axes, and title
## so that the plot() can use those
## and pairs can utilize as well
pairs.snd <-
function (x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc,
    lower.panel = panel, upper.panel = panel, diag.panel = NULL,
    text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3,
    cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1,
    log = "")
{
    if (doText <- missing(text.panel) || is.function(text.panel))
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x,
            y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main,
        oma, ...) {
        xpd <- NA
        if (side%%2L == 1L && xl[j])
            xpd <- FALSE
        if (side%%2L == 0L && yl[i])
            xpd <- FALSE
        if (side%%2L == 1L)
            Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]]))
                x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    }
    else if (!is.numeric(x))
        stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel))
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2L)
        stop("only one column in the argument to 'pairs'")
    if (!all(horInd >= 1L && horInd <= nc))
        stop("invalid argument 'horInd'")
    if (!all(verInd >= 1L && verInd <= nc))
        stop("invalid argument 'verInd'")
    if (doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels))
                labels <- paste("var", 1L:nc)
        }
        else if (is.null(labels))
            doText <- FALSE
    }
    oma <- if ("oma" %in% nmdots)
        dots$oma
    main <- if ("main" %in% nmdots)
        dots$main
    if (is.null(oma))
        oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
    opar <- par(mfrow = c(length(horInd), length(verInd)), mar = rep.int(gap/2,
        4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    xl <- yl <- logical(nc)
    if (is.numeric(log))
        xl[log] <- yl[log] <- TRUE
    else {
        xl[] <- grepl("x", log)
        yl[] <- grepl("y", log)
    }
    for (i in if (row1attop)
        verInd
    else rev(verInd)) for (j in horInd) {
        l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y",
            ""))
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE,
            type = "n", ..., log = l)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2L) || !has.upper || !has.lower))
                localAxis(1L + 2L * row1attop, x[, j], x[, i],
                  ...)
            if (i == nc && (j%%2L || !has.upper || !has.lower))
                localAxis(3L - 2L * row1attop, x[, j], x[, i],
                  ...)
            if (j == 1 && (!(i%%2L) || !has.upper || !has.lower))
                localAxis(2L, x[, j], x[, i], ...)
            if (j == nc && (i%%2L || !has.upper || !has.lower))
                localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag)
                  localDiagPanel(as.vector(x[, i]), ...)
                if (doText) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  xlp <- if (xl[i])
                    10^0.5
                  else 0.5
                  ylp <- if (yl[j])
                    10^label.pos
                  else label.pos
                  text.panel(xlp, ylp, labels[i], cex = cex.labels,
                    font = font.labels)
                }
            }
            else if (i < j)
                localLowerPanel(as.vector(x[, j]), as.vector(x[,
                  i]), ...)
            else localUpperPanel(as.vector(x[, j]), as.vector(x[,
                i]), ...)
            if (any(par("mfg") != mfg))
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots)
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots)
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main,
            font = font.main)
    }
    invisible(NULL)
}



## kernel density estimates

library(KernSmooth)
#library(ash)

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
border=NA, add=FALSE, horiz=TRUE, limit, ...)
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
    if (missing(limit))
        limit <- range(z$x)
    ii <- z$x >= limit[1] & z$x <= limit[2]
    z$x <- z$x[ii]
    z$y <- z$y[ii]

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
    d2$bw <- bw
    d2
}


biplot <- function(x, y, col=1, type="d", ...) {
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
    dplot(x, horiz=TRUE, add=TRUE, type=type, col=Pal[33], border=NA,
        at=u[4], shift=0, w=0.05*diff(u[3:4]), limit=u[1:2])
    qbox(x, type=" I", horiz=TRUE, col=Pal[66],
        col_med="#FFFFFF", col_box=Pal[66], lwd=2, lwd_med=3,
        at=u[4], w=0.01*diff(u[3:4]), add=TRUE)
    dplot(y, horiz=FALSE, add=TRUE, type=type, col=Pal[33], border=NA,
        at=u[2], shift=0, w=0.05*diff(u[1:2]), limit=u[3:4])
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

biplot(iris[,1], iris[,2])
biplot(iris[,3], iris[,4])

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

col=1
Pal <- colorRampPalette(c("#FFFFFF", col))(100)
plot(c(1-0.5, m+0.5), range(unlist(xx)), ann=FALSE, axes=FALSE, type="n")
axis(1, seq_len(m), names(xx), tick=FALSE)
axis(2)
for (i in seq_len(m)) {
    dplot(xx[[i]], horiz=FALSE, at=i, w=0.45, type="d", add=TRUE, col=Pal[50])
    qbox(xx[[i]], type=" I", horiz=FALSE, col=Pal[50],
        col_med=Pal[100], col_box=Pal[1], lwd=2, lwd_med=3,
        at=i, w=0.05, add=TRUE)
}
#box(col=Pal[50])
rug(seq_len(m), -0.03, lwd=1)



library(microbenchmark)

## taken from https://github.com/juba/questionr
Cramer_V <- function(x) {
  chid <- stats::chisq.test(x, correct=FALSE)$statistic
  dim <- min(nrow(x),ncol(x)) - 1
  as.numeric(sqrt(chid / (sum(x) * dim)))
}

## chisq.test based value takes 5.6x longer
x1 <- sample(LETTERS[1:4], 100, replace=TRUE, prob=4:1)
x2 <- sample(LETTERS[1:4], 100, replace=TRUE, prob=4:1)
x <- table(x1, x2)
b <- microbenchmark(
    this=cramer_v(x),
    questionr=Cramer_V(x)
)

## kde is 3x faster
x <- rnorm(10000)
b <- microbenchmark(
    fft=density(x),
    kde=KernSmooth::bkde(x)
)
library(ggplot2)
autoplot(b)
