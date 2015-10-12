test_intersection <- function(a,b){
    #reads in two rasters and tests for overlap T or F
    # if returns TRUE then there is overlap
    !class(try(intersect(a,b),T ))=='try-error'
}


