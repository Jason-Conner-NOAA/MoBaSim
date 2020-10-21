# To get a sample of random data that follows a log-normal distribution and has arithmetic mean of 7 and a 
# standard deviation of 75, you need to reparameterize things.  Roughly, you need to figure out what parameters 
# should go into to the normal distribution such that when you exponentiate the draws, you end up with a mean of 7 
# and a standard deviation of 75. 

m <- Data$WEIGHT[100]
s <- sqrt(var(Data$WEIGHT))
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
print(paste("location:", location))
print(paste("shape:", shape))
draws3 <- rlnorm(n=10000, location, shape)
plot(density(draws3))


library(sumfish)
library(FNN)
getSQL() 

rb <- getRacebase(2019,'EBS_SHELF')
haul <- sumHaul(rb)

atf <- filter(haul, SPECIES_CODE==10110)

nn <- get.knn(select(atf,START_LONGITUDE,START_LATITUDE),k=4)

for (i in 1:nrow(atf)) {
  nn.sd <- sqrt(var(atf[nn$nn.index[i,],]$WEIGHT))
  if (nn.sd > 0) {
    m <- atf$WEIGHT[i]
    s <- nn.wt
    location <- log(m^2 / sqrt(s^2 + m^2))
    shape <- sqrt(log(1 + (s^2 / m^2)))
    print(paste("location:", location))
    print(paste("shape:", shape))
    
    newWt <- rlnorm(n=1, location, shape)
  }
}
