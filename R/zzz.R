#Copyright 2004, W. Wolski, all rights reserved.
.First.lib <- function(lib, pkg) library.dynam("pairseqsim",pkg,lib)
.Last.lib <- function(libpath) library.dynam.unload("pairseqsim", libpath)

