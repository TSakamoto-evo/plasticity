data <- read.table("fitness.txt")

head(data)
sep <- 1
index <- seq(1, dim(data)[1], sep)

pop_size <- 10000

plot(NA, xlim=c(0, 2e+8),
     type="l", ylim=c(0, 1), xlab="", ylab="", xaxt="n", yaxt="n")

axis(1, at=c(0, 1e+8, 2e+8))
axis(2, at=c(0, 0.5, 1.0))
points(data[index,1], data[index,2], col="red", type="l")
points(data[index,1], data[index,3], col="blue", type="l")

