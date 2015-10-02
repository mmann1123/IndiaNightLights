 # compute squared euclidean distance from each sample to each cluster center

clusters <- function(x, centers) {

  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  max.col(-t(tmp))  # find index of min distance
}



# Example
# create a simple data set with two clusters
#set.seed(1)
#x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#colnames(x) <- c("x", "y")
#x_new <- rbind(matrix(rnorm(10, sd = 0.3), ncol = 2),
#               matrix(rnorm(10, mean = 1, sd = 0.3), ncol = 2))
#colnames(x_new) <- c("x", "y")
#
#cl <- kmeans(x, centers=2)
#
#all.equal(cl[["cluster"]], clusters(x, cl[["centers"]]))
## [1] TRUE
#clusters(x_new, cl[["centers"]])
## [1] 2 2 2 2 2 1 1 1 1 1
#
#plot(x, col=cl$cluster, pch=3)
#points(x_new, col= clusters(x_new, cl[["centers"]]), pch=19)
#points(cl[["centers"]], pch=4, cex=2, col="blue")

