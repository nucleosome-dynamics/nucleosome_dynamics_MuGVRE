#!/usr/bin/Rscript

myFilter <- function(x, f, ...)
    x[f(x, ...)]

partial <- function(f, ...)
{
    capture <- list(...)
    function(x) do.call(f, c(list(x), capture))
}

compose <- function(...)
{
    comp2 <- function(f, g) {
        force(f)
        force(g)
        function(x) f(g(x))
    }
    Reduce(comp2, list(...))
}

flip2args <- function(f)
    function(x, y) f(y, x)

procFun <- function(f, g)
{
    force(f)
    g(f)
}
