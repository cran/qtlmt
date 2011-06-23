.First.lib <- function(lib, pkg) {
   library.dynam("qtlmt", pkg, lib)
}
.onLoad <- function(lib, pkg) print("R/qtlmt is loaded")
.noGenerics <- TRUE
.onUnload <- function(libpath) library.dynam.unload("qtlmt", libpath)

