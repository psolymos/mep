## siplot: selection index plot, without model fit

siplot <- function (object, ...)
    UseMethod("siplot")

siplot.default <-
function(y, x, which=NULL, ask, ylab, subset=NULL, ...)
{
    if (NCOL(x) < 2L)
        x <- as.data.frame(x=x)
    object <- .get_frame(formula=y ~ ., data=x, type="numeric")
    .mep_engine(object, sip=TRUE,
        which=which, ask=ask, ylab=ylab, subset=subset, ...)
}

siplot.formula <-
function(formula, data, which=NULL, ask, ylab, subset=NULL, ...)
{
    object <- .get_frame(formula=formula, data=data, type="numeric")
    .mep_engine(object, sip=TRUE,
        which=which, ask=ask, ylab=ylab, subset=subset, ...)
}

