num_g <- 8

data <- read.table("regi_genotype.txt")

head(data)

max_value_g <- 5

max_g <- 5
time <- 5e+7

rm(subdata)
subdata <- data[data$V1 == time,]
subdata <- subdata[order(subdata$V2, decreasing=T)[1], ]

plot(NA, xlim=c(0, num_g*2+2), ylim=c(-25.5, 0.5), axes=F, ann=F)


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

for(i in 1:26){
  for(j in 1:num_g){
    color1 <- rgb(regi_p1[j, i] / max_value_g, 0, 0)
    color2 <- rgb(regi_p2[j, i] / max_value_g, 0, 0)
    
    polygon(c(j-1, j-1, j, j), c(-i+1.5, -i+0.5, -i+0.5, -i+1.5), col=color1)
    polygon(c(num_g+j+2, num_g+j+2, num_g+j+1, num_g+j+1), c(-i+1.5, -i+0.5, -i+0.5, -i+1.5), col=color2)
  } 
}
