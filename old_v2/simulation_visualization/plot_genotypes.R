num_g <- 8

data <- read.table("regi_genotype.txt")

head(data)

max_value_g <- 5
max_value_e <- max(abs(data[,(num_g+3):(2*num_g+2)]))
max_value_b <- max(abs(data[,(2*num_g+3):(num_g**2+2*num_g+2)]))

max_g <- 5
time <- 5e+7

rm(subdata)
subdata <- data[data$V1 == time,]
subdata <- subdata[order(subdata$V2, decreasing=T)[1], ]

plot(NA, xlim=c(0, num_g + 7.5), ylim=c(0, num_g), axes=F, ann=F)

if(num_g > 0){
  for(j in 1:num_g){
    color <- rgb(subdata[1, 2+j]/max_value_g, 0, 0)
    polygon(c(num_g+1.5, num_g+1.5, num_g+2.5, num_g+2.5), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color)
  }
}

for(j in 1:num_g){
  index <- num_g + 2 + j
  if(subdata[index] < 0){
    color <- rgb(0, -subdata[1, index]/max_value_e, 0)
  }else{
    color <- rgb(subdata[1, index]/max_value_e, 0, 0)
  }
  polygon(c(num_g+0.25, num_g+0.25, num_g+1.25, num_g+1.25), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color)
}

for(j in 0:(num_g - 1)){
  for(k in 0:(num_g - 1)){
    index <- 2 * num_g + j * num_g + k + 3
    if(subdata[1, index] < 0){
      color <- rgb(0, -subdata[1, index]/max_value_b, 0)
    }else{
      color <- rgb(subdata[1, index]/max_value_b, 0, 0)
    }
    polygon(c(k, k, k+1, k+1), c(num_g-j, num_g-j-1, num_g-j-1, num_g-j), col=color)
  }
}

subdata <- as.matrix(subdata)
g_vec <- subdata[1, 3:(num_g+2)]
c_vec <- subdata[1, (num_g+3):(2*num_g+2)]
b_mat <- subdata[,(2*num_g+3):(num_g**2+2*num_g+2)]
b_mat <- matrix(b_mat, nrow=num_g, ncol=num_g, byrow=T)

p1 <- g_vec
p2 <- g_vec

regi_p1 <- p1
regi_p2 <- p2

for(dev in 1:10){
  effect1 <- b_mat %*% p1 + max_g * c_vec
  effect2 <- b_mat %*% p2
  
  p1 <- 0.8 * p1 + (0.5 + 0.5 * tanh(effect1))
  p2 <- 0.8 * p2 + (0.5 + 0.5 * tanh(effect2))
  
  regi_p1 <- cbind(regi_p1, p1)
  regi_p2 <- cbind(regi_p2, p2)
}

for(dev in 11:25){
  effect1 <- b_mat %*% p1
  effect2 <- b_mat %*% p2
  
  p1 <- 0.8 * p1 + (0.5 + 0.5 * tanh(effect1))
  p2 <- 0.8 * p2 + (0.5 + 0.5 * tanh(effect2))
  
  regi_p1 <- cbind(regi_p1, p1)
  regi_p2 <- cbind(regi_p2, p2)
}

for(j in 1:num_g){
  color1 <- rgb(p1[j, 1] / max_value_g, 0, 0)
  color2 <- rgb(p2[j, 1] / max_value_g, 0, 0)
  
  polygon(c(num_g+2.75, num_g+2.75, num_g+3.75, num_g+3.75), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color1)
  polygon(c(num_g+4, num_g+4, num_g+5, num_g+5), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color2)
}

for(k in 0:100){
  p1 <- g_vec
  regi_p1 <- p1

  for(dev in 1:10){
    effect1 <- b_mat %*% p1 + max_g * c_vec * (100-k)/100
    p1 <- 0.8 * p1 + (0.5 + 0.5 * tanh(effect1))
    regi_p1 <- cbind(regi_p1, p1)
  }
  
  for(dev in 11:25){
    effect1 <- b_mat %*% p1
    p1 <- 0.8 * p1 + (0.5 + 0.5 * tanh(effect1))
    regi_p1 <- cbind(regi_p1, p1)
  }
  
  for(j in 1:num_g){
    color1 <- rgb(p1[j, 1] / max_value_g, 0, 0)
    polygon(c(k*2.0/101+num_g+5.5, k*2.0/101+num_g+5.5, (k+1)*2.0/101+num_g+5.5, (k+1)*2.0/101+num_g+5.5), 
            c(num_g-j, num_g-j+0.6, num_g-j+0.6, num_g-j), col=color1, border=color1, lwd=0.1)
  }
}

