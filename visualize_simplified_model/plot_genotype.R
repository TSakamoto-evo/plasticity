data <- read.table("simulation/ind_regi.txt")
head(data)

cue <- 5

time <- 100000
subdata <- as.matrix(data[data$V1 == time, ])
num_g <- subdata[2]

g <- subdata[1, 3:(2+num_g)]
b <- subdata[1, (3+num_g):(2+num_g*(num_g+1))]
b <- matrix(b, num_g, num_g, byrow=T)
c <- subdata[1, (3+num_g*(num_g+1)):(2+num_g*(num_g+2))]

data2 <- read.table("basin.txt")
head(data2)

#postscript("tmp.eps", width=2.85, height=3.2)
plot(NA, xlim=c(0, 5), ylim=c(0, 5))
for(i in 0:20){
  for(j in 0:20){
    x1 <- i / 4
    x2 <- j / 4
    
    px1 <- (1-0.2)*x1 + 0.5 + 0.5 * tanh(b[1, 1]*x1 + b[1, 2]*x2)
    px2 <- (1-0.2)*x2 + 0.5 + 0.5 * tanh(b[2, 1]*x1 + b[2, 2]*x2)
    
    arrows(x1, x2, px1, px2, length=0.02, col=gray(0.5), lwd=0.5)
  }
}

points(g[1], g[2], cex=1.5, lwd=2, pch=4)

if(length(unique(data2$V3)) == 2){
  regi_x <- numeric(0)
  regi_y <- numeric(0)
  
  sep <- data2$V2[2] - data2$V2[1]
  
  for(i in unique(data2$V1)){
    subdata <- data2[data2$V1 == i, ]
    equ <- unique(subdata$V3)
    
    if(length(equ) > 1){
      equ_diff <- subdata$V3[-1] - subdata$V3[-length(subdata$V3)]
      y <- subdata$V2[-1][equ_diff != 0] - sep / 2
      
      regi_x <- c(regi_x, i)
      regi_y <- c(regi_y, y)
    }
  }
  
  l <- smooth.spline(regi_x, regi_y)
  points(predict(l, regi_x), type="l", lwd=2, lty="dashed")
  
  data3 <- read.table("equilibria.txt")
  if(data3$V2[1] < data3$V2[2]){
    points(data3$V2[1], data3$V3[1], col = 4, cex=1, lwd=2)
    points(data3$V2[2], data3$V3[2], col = 2, cex=1, lwd=2)
  }else{
    points(data3$V2[1], data3$V3[1], col = 2, cex=1, lwd=2)
    points(data3$V2[2], data3$V3[2], col = 4, cex=1, lwd=2)
  }
}else{
  data3 <- read.table("equilibria.txt")
  points(data3$V2[1], data3$V3[1], col = "purple", cex=1, lwd=2)
}

p1 <- c(g[1])
p2 <- c(g[2])

p <- g

for(i in 1:40){
  effect <- b %*% p
  p <- 0.8 * p + 0.5 + 0.5 * tanh(effect)
  p1 <- c(p1, p[1])
  p2 <- c(p2, p[2])
}

points(p1, p2, type="l", col=4, lwd=1.5)

p1 <- c(g[1])
p2 <- c(g[2])
p <- g

for(i in 1:10){
  effect <- b %*% p + c*cue
  p <- 0.8 * p + 0.5 + 0.5 * tanh(effect)
  p1 <- c(p1, p[1])
  p2 <- c(p2, p[2])
}

points(p1[length(p1)], p2[length(p2)], pch=4, col=2)

for(i in 11:40){
  effect <- b %*% p
  p <- 0.8 * p + 0.5 + 0.5 * tanh(effect)
  p1 <- c(p1, p[1])
  p2 <- c(p2, p[2])
}

points(p1, p2, type="l", col=2, lwd=1.5)
#dev.off()