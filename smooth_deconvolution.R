source("methods.R")
# Generic functions
fit = function(obj, ...) {
    UseMethod("fit", obj)
}
plot = function(obj, ...) {
    UseMethod("plot", obj)
}
initialize = function(obj) {
    UseMethod("initialize", obj)
}