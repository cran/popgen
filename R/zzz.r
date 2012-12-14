
.onLoad <- function(lib, pkg){
     # do whatever needs to be done when the package is loaded
      library.dynam("popgen", pkg, lib)
}
