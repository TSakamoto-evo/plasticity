data <- read.table("regi_fitness.txt")

head(data)

plot(data[,1], data[,2], col="red", type="l", ylim=c(0, 1), xlab="time", ylab="fitness")
points(data[,1], data[,3], col="blue", type="l")