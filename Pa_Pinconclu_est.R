nsim.i = 5000

n = 1; w = 2
k = rep(kn, n)

Tmean = list()
for(i in 1:nsim.i){
  Tmean[[i]] = rmultinom(n, 1, p)
}

TT = matrix(0, nrow = nsim.i, ncol = n)
class.save = matrix(0, nrow = nsim.i, ncol = n)

for(i in 1:nsim.i){
  for(nn in 1:n){
    class = which( Tmean[[i]][,nn] != 0)
    if(class == 1){
      TT[i,nn] = T_ent
      class.save[i, nn] = -1
    }
    if(class == 2){
      TT[i,nn] = 0 
      class.save[i, nn] = 0
    }
    if(class == 3){
      TT[i,nn] = T_ent
      class.save[i, nn] = 1
    }
  }
  
}

entry = list()
c_time_list = list()
for(j in 1:nsim.i){
  entry[[j]] = seq(0, 2000, length.out = n)
  c_time = list()
  for(i in 1:n){
    t = rep((2 + 0.5), (k[i] - 1))
    c_time[[i]] = cumsum(c((entry[[j]][i] + 2), t))
  }
  c_time_list[[j]] = c_time
}

interim = length(c_time_list[[1]][[1]] - 1)

ZZk = ldBounds(t = (2:(kn))/(kn), t2 = (2:(kn))/(kn), 
               iuse = 1, alpha = 0.05, sides = 1)$upper.bounds
ZZk = ZZk[length(ZZk)]
z_value = ZZk*sqrt(kn/c(1:kn))

workers = 10

cl <- makeCluster(workers) 

registerDoParallel(cl)


out <- foreach(i=1:nsim.i) %dopar%{
  run = myfun(power, alpha, method = "OBF", abs_delta = abs_delta, gamma = gamma, 
              TT = TT[i,], k = k, n = n, entry = entry[[i]], c_time = c_time_list[[i]], class.save.in = class.save[i,],
              z_value = z_value, pstar = 2)
  
  run
}

stopCluster(cl)

all = do.call(cbind, out)
stop_for = do.call(rbind, all[1,])
stop_glb = unlist(all[2,])
stop_slot = do.call(rbind, all[3,])
bin_check = do.call(rbind, all[4,])

rm.ind = NULL

if(length(rm.ind)==0){
  rm.TT = as.numeric(class.save)
  stop_for = stop_for
}else{
  rm.TT = as.numeric(class.save)[-rm.ind]
  stop_for = stop_for[-as.numeric(rm.ind)]
}

bin = ifelse(as.numeric(stop_for) == as.numeric(rm.TT), 1, 0 )

who_rm = as.numeric(TT)[rm.ind]
compare = data.frame(truth = rm.TT, stop_for = as.numeric(stop_for), 
                     decision = as.numeric(bin))


wrong = A4 = B4 = E4 = NULL
for(i in 1:nrow(compare)){
  if(compare[i,]$stop_for != 4 && compare[i,]$decision == 0){
    wrong = c(wrong, i)
  }
  if(compare[i,]$stop_for == 4 && compare[i,]$decision == 0){
    if(compare[i,]$truth == 1){
      A4 = c(A4, i)
    }
    if(compare[i,]$truth == -1){
      B4 = c(B4, i)
    }
    if(compare[i,]$truth == 0){
      E4 = c(E4, i)
    }
  }
}


