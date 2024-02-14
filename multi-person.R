library(binom);library(MASS);library(ldbounds)

source("./est_fun.R") # Load the main function for estimation in n-of-1

# p = probability of belonging to the 3 classes: A<B, A=B, A>B
p = c(0, 0.2, 0.8) 

sig = 3.54 # sigma

abs_delta = 2.1 #|delta|

mu = T_ent = 6.9  #\tau

power = 0.8
alpha = 0.05

pa = 0.96 #pa3_star
p0 = 0.46 #p0_sta
p.inconc = 0.33 # probability of inconclusive

n_pre = 1
success = F
while(success == F){
  dummy = c(0:n_pre)
  ind1 = c()
  for(i in 1:length(dummy)){
    ind1[i] = dbinom(dummy[i], size = n_pre, prob = p0)
  }
  
  cum = cumsum(ind1)
  w = min(which(cum > (1-alpha)))
  
  if(!is.infinite(w)){
    y = c(0:dummy[w])
    
    ind = c()
    for(i in 1:length(y)){
      ind[i] = dbinom(y[i], size = n_pre, prob = pa)
    }
    if( (1 - sum(ind)) < (power)){
      success = F
      n_pre = n_pre + 1
    }else{
      success = T
      break
    }
  }else{
    success = F
    n_pre = n_pre + 1
  }

}

w = dummy[w]
cat('w:\n', w)
n = (n_pre/(1-p.inconc))
cat('number of individuals:\n', n)
