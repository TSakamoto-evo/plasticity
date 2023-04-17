num_g <- 8

data <- read.table("regi_genotypes.txt")
max_value_g <- 5
max_value_b <- max(abs(data[,(num_g+2):(1+num_g+num_g**2)]))
max_value_e <- max(abs(data[,(num_g+2+num_g**2):(1+2*num_g+num_g**2)]))

row <- 2001

subdata <- data[row, -1]

#postscript("tmp.eps", width = 3.8, height=2.9)
plot(NA, xlim=c(0, num_g + 5.0), ylim=c(0, num_g), axes=F, ann=F)

if(num_g > 0){
  for(j in 1:num_g){
    color <- rgb(subdata[j]/max_value_g, 0, 0)
    #polygon(c(j-1, j-1, j, j), c(-1.5, -0.5, -0.5, -1.5), col=color)
    polygon(c(num_g+1.5, num_g+1.5, num_g+2.5, num_g+2.5), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color)
  }
}

for(j in 0:(num_g - 1)){
  for(k in 0:(num_g - 1)){
    index <- num_g + j * num_g + k + 1
    if(subdata[index] < 0){
      color <- rgb(0, -subdata[index]/max_value_b, 0)
    }else{
      color <- rgb(subdata[index]/max_value_b, 0, 0)
    }
    polygon(c(k, k, k+1, k+1), c(num_g-j, num_g-j-1, num_g-j-1, num_g-j), col=color)
  }
}

for(j in 1:num_g){
  index <- num_g + num_g ** 2 + j
  if(subdata[index] < 0){
    color <- rgb(0, -subdata[index]/max_value_e, 0)
  }else{
    color <- rgb(subdata[index]/max_value_e, 0, 0)
  }
  polygon(c(num_g+0.25, num_g+0.25, num_g+1.25, num_g+1.25), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color)
  #polygon(c(j-1, j-1, j, j), c(-1.5, -0.5, -0.5, -1.5), col=color)
}

for(j in 1:num_g){
  index <- num_g + num_g ** 2 + num_g + j
  if(subdata[index] < 0){
    color <- rgb(0, -subdata[index]/max_value_g, 0)
  }else{
    color <- rgb(subdata[index]/max_value_g, 0, 0)
  }
  polygon(c(num_g+2.75, num_g+2.75, num_g+3.75, num_g+3.75), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color)
  #polygon(c(j-1, j-1, j, j), c(-1.5, -0.5, -0.5, -1.5), col=color)
}

for(j in 1:num_g){
  index <- num_g + num_g ** 2 + 2*num_g + j
  if(subdata[index] < 0){
    color <- rgb(0, -subdata[index]/max_value_g, 0)
  }else{
    color <- rgb(subdata[index]/max_value_g, 0, 0)
  }
  polygon(c(num_g+4, num_g+4, num_g+5, num_g+5), c(num_g-j, num_g-j+1, num_g-j+1, num_g-j), col=color)
  #polygon(c(j-1, j-1, j, j), c(-1.5, -0.5, -0.5, -1.5), col=color)
}

#dev.off()


print(data[row, 1])
print(max(abs(data[,(num_g+2):(1+num_g+num_g**2)])))
print(max(abs(data[,(num_g+2+num_g**2):(1+2*num_g+num_g**2)])))

