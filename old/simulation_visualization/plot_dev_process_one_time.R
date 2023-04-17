num_g <- 8
dev_t <- 10
alpha <- 0.2
max_value <- 1/alpha

### specify the line in "regi_genotype.txt"
row <- 2001

data <- read.table("regi_genotypes.txt")
subdata <- data[row, -1]

g0 <- subdata[1:num_g]
subdata <- subdata[-(1:num_g)]

b <- subdata[1:(num_g**2)]
b_mat <- matrix(b, nrow=num_g, ncol=num_g, byrow=T)
subdata <- subdata[-(1:(num_g**2))]

e <- subdata[1:num_g]

#postscript("tmp.eps", height=10, width=3.5)
plot(NA, xlim=c(-1, num_g), ylim=c(0, 2*(dev_t+4)), axes=F, ann=F)

for(j in 0:dev_t){
  text(-1, 2*dev_t+4.65-j, j)
  text(-1, dev_t+0.65-j, j)
}

text(-1, 2*dev_t+5.5, "env 1")
text(-1, dev_t+2, "env 2")

text((num_g-1)/2, 2*dev_t+6.5, paste("t =", t))

delta <- 1

g <- g0
g_tmp <- g

for(i in 1:num_g){
  color <- rgb(g_tmp[i]/max_value, 0, 0)
  polygon(c(i-1, i-1, i, i), c(2*dev_t+4.7-0, 2*dev_t+4-0, 2*dev_t+4-0, 2*dev_t+4.7-0), col=color)
}

for(i in 1:10){
  for(j in 1:num_g){
    effect <- 0.0
    for(k in 1:num_g){
      effect <- effect + b_mat[j, k] * g[k]
    }
    effect <- effect + e[j] * delta
    g_tmp[j] <- g[j] + (0.5 + 0.5 * tanh(effect)) - alpha * g[j]
  }
  g <- g_tmp
  
  for(j in 1:num_g){
    color <- rgb(g_tmp[j]/max_value, 0, 0)
    polygon(c(j-1, j-1, j, j), c(2*dev_t+4.7-i, 2*dev_t+4-i, 2*dev_t+4-i, 2*dev_t+4.7-i), col=color)
  }
}
#dev.off()

delta <- 0

g <- g0
g_tmp <- g

for(i in 1:num_g){
  color <- rgb(g_tmp[i]/max_value, 0, 0)
  polygon(c(i-1, i-1, i, i), c(dev_t+0.7-0, dev_t-0, dev_t-0, dev_t+0.7-0), col=color)
}

for(i in 1:10){
  for(j in 1:num_g){
    effect <- 0.0
    for(k in 1:num_g){
      effect <- effect + b_mat[j, k] * g[k]
    }
    effect <- effect + e[j] * delta
    g_tmp[j] <- g[j] + (0.5 + 0.5 * tanh(effect)) - alpha * g[j]
  }
  g <- g_tmp
  
  for(j in 1:num_g){
    color <- rgb(g_tmp[j]/max_value, 0, 0)
    polygon(c(j-1, j-1, j, j), c(dev_t+0.7-i, dev_t-i, dev_t-i, dev_t+0.7-i), col=color)
  }
}

#dev.off()

print(data[row, 1])
