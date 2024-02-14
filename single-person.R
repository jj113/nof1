library(binom);library(MASS);library(ldbounds)

source("./est_fun.R") # Load the main function for estimation in n-of-1

# p = probability of belonging to the 3 classes: A<B, A=B, A>B
p = c(0, 0, 1) 

sig = 3.54 # sigma

abs_delta = 2.1 #|delta|

mu = T_ent = 6.9  #\tau

power = 0.8
alpha = 0.05

#sample size calculation for 1 person
kn = single.sample(idx.p = idx.p, power = power, alpha = alpha, mu = mu, 
                   abs_delta = abs_delta, sig = sig)

cat('number of cycles per person:', kn)




