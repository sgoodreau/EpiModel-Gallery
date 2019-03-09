###

diss_prob_non <- 0.01
diss_prob_chb <- 0.001
n2c_prob <- 0.03

n <- 1000
timesteps <- 5e3

states <- array(NA, dim=c(1, 3, timesteps+1))
states[,,1] <- matrix(c(n, 0, 0), nrow=1) ## Number of rels that are noncohab (N), cohbab (C), or terminated (T)

transmat <- matrix(c(
  1-diss_prob_non-n2c_prob, n2c_prob, diss_prob_non,
  0, 1-diss_prob_chb, diss_prob_chb,
  0,0, 1
), 3, 3, byrow=TRUE)

for (i in 2:(timesteps+1)) states[,,i] <- states[,,i-1] %*% transmat

plot(states[,1,1:200])
