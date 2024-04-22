data <- read.table("simulation/regi_fit.txt")
head(data)

#postscript("tmp.eps", width=4, height=3.2)
#postscript("tmp.eps", width=4, height=2.4)

plot(NA, type="l", xlim=c(0, 5e+5), ylim=c(0, 1), xlab="", ylab="", xaxt="n", yaxt="n")
#plot(NA, type="l", xlim=c(0e+5, 0.6e+5), ylim=c(0, 1.3), xlab="", ylab="", yaxt="n")

axis(1, at=c(0, 2.5e+5, 5e+5))
axis(2, at=c(0, 0.5, 1.0))

points(data[,1], data[,2], col="red", type="l")
points(data[,1], data[,3], col="blue", type="l")

abline(v=0, lty="dashed")
abline(v=5e+5, lty="dashed")
abline(v=1e+5, lty="dashed")

#abline(h=0.95)
abline(v=1.8e+5)
#dev.off()
