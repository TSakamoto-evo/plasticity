data <- read.table("simulation/ind_regi.txt")
head(data)

cue <- 5
step <- 30

max_value_g <- 5
max_value_b <- 0.48
max_value_c <- 0.43

time <- 400000
subdata <- as.matrix(data[data$V1 == time, ])
num_g <- subdata[2]

g <- subdata[1, 3:(2+num_g)]
b <- subdata[1, (3+num_g):(2+num_g*(num_g+1))]
b <- matrix(b, num_g, num_g, byrow=T)
c <- subdata[1, (3+num_g*(num_g+1)):(2+num_g*(num_g+2))]

#postscript("tmp.eps", width = 4.3, height=2.9)
plot(NA, xlim=c(0, num_g + 7.5), ylim=c(0, num_g), axes=F, ann=F)
for(i in 1:num_g){
  polygon(c(num_g+2.1, num_g+2.1, num_g+3.1, num_g+3.1), 
          c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), 
          col=rgb(g[i]/max_value_g, 0, 0), border="white", lwd=0.8)
  
  if(c[i] >= 0){
    col <- rgb(c[i]/max_value_c, 0, 0)
  }else{
    col <- rgb(0, -c[i]/max_value_c, 0)
  }
  polygon(c(num_g+0.6, num_g+0.6, num_g+1.6, num_g+1.6), 
          c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), col=col, 
          border="white", lwd=0.8)
}

for(i in 1:num_g){
  for(j in 1:num_g){
    if(b[i, j] >= 0){
      col <- rgb(b[i, j]/max_value_b, 0, 0)
    }else{
      col <- rgb(0, -b[i, j]/max_value_b, 0)
    }
    polygon(c(j, j, j-1, j-1), c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), 
            col=col, border="white", lwd=0.8)
  }
}

p <- g
for(i in 1:step){
  effect <- b %*% p
  p <- 0.8 * p + 0.5 + 0.5 * tanh(effect)
}
for(i in 1:num_g){
  polygon(c(num_g+5.85, num_g+5.85, num_g+6.85, num_g+6.85), 
          c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), 
          col=rgb(p[i]/max_value_g, 0, 0), border="white", lwd=0.8)
}

p <- g
for(i in 1:10){
  effect <- b %*% p + c*cue
  p <- 0.8 * p + 0.5 + 0.5 * tanh(effect)
}

for(i in 11:step){
  effect <- b %*% p
  p <- 0.8 * p + 0.5 + 0.5 * tanh(effect)
}
for(i in 1:num_g){
  polygon(c(num_g+4.1, num_g+4.1, num_g+5.1, num_g+5.1), 
          c(num_g-i, num_g-i+1, num_g-i+1, num_g-i), 
          col=rgb(p[i]/max_value_g, 0, 0), border="white", lwd=0.8)
}

#dev.off()

timedata <- data[data$V1 <= time, ]
g <- timedata[, 3:(2+num_g)]
b <- timedata[, (3+num_g):(2+num_g*(num_g+1))]
c <- timedata[, (3+num_g*(num_g+1)):(2+num_g*(num_g+2))]

max(abs(g))
max(abs(b))
max(abs(c))
