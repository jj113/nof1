# level = decimals to keep
floor_dec <- function(x, level=2) round(x - 5*10^(-level-1), level)

single.sample = function(idx.p, power, alpha, mu, abs_delta, sig){
  
  if(idx.p == 1 | idx.p == 3){
    If2 = (qnorm(power) + qnorm(1-alpha))^2/(mu-abs_delta)^2
    n.no.Rp = floor(2*sig^2*If2)
    interim = n.no.Rp
    kn.no.Rp = interim + 1
    ZZk = ldBounds(t = (2:(kn.no.Rp))/(kn.no.Rp), t2 = (2:(kn.no.Rp))/(kn.no.Rp), 
                   iuse = 1, alpha = alpha, sides = 1)
    
    Rp = sqrt(ZZk$upper.bounds[length(ZZk$upper.bounds)]/qnorm(1-alpha)) 
    kn = ceiling(n.no.Rp*Rp)
  }
  if(idx.p == 2){
    b0 = 1-power
    If2 = (qnorm(1-b0/2) + qnorm(1-alpha))^2/(0-abs_delta)^2
    n.no.Rp = floor(2*sig^2*If2)
    
    interim = n.no.Rp
    kn.no.Rp = interim + 1
    
    ZZk = ldBounds(t = (2:(kn.no.Rp))/(kn.no.Rp), t2 = (2:(kn.no.Rp))/(kn.no.Rp), 
                   iuse = 1, alpha = alpha, sides = 1)
    
    Rp = sqrt(ZZk$upper.bounds[length(ZZk$upper.bounds)]/qnorm(1-alpha)) 
    
    kn = ceiling(n.no.Rp*Rp)
  }
  return(kn)
}

myfun <- function(power, alpha, method, abs_delta, gamma, TT, k, n, n_pre, entry, c_time,class.save.in, z_value, pstar){
  
  int_time = sort(unique(unlist(c_time)))

  d = list()
  for(i in 1:n){
    d[[i]] = rnorm(k[i], mean = TT[i], sd = sqrt(2*sig^2))
  }
  

  t_hat = c()
  ci = as.list(rep(0, n))
  t_stat = c()
  
  ind_cum = matrix(0, nrow = n, ncol = max(k))
  
  interm = as.list(rep(4, n))
  stop_at = rep(0, n)
  sum_g = 1
  who_active = NULL
  who_na = NA
  SUCCESS_G = F
  for(i in 1:length(unique(int_time))){
    ck = cumsum(k)
    ppl = which(int_time[i] == unlist(c_time))
    who_all = c()
    for(j in 1:length(ppl)){
      who_all[j] = which(ck >= ppl[j])[1]
    }
    
    for(j in 1:length(ppl)){
      who = who_all[j]
      
      if(all(interm[[who]] == 4)){
        prev = length(which(ind_cum[who,] == 1))
        ind_cum[who,] = c(rep(1, (prev + 1)), rep(0, (ncol(ind_cum) - prev - 1)) ) 
      }else{
        prev = 0
      }

      if((prev >= 1) && all(interm[[who]] == 4) && (prev < k[who]) && !(who %in% who_na)){
        
        who_active = c(who_active, who)
        
        
        diff = d[[who]][1:(prev + 1)]
        t_hat = sum(diff)/length(diff)
        
        lower = t_hat - z_value[prev]*(sd(diff)/sqrt(length(diff)))
        upper = t_hat + z_value[prev]*(sd(diff)/sqrt(length(diff)))
        
        If2 = 1/((abs_delta - (t_hat))^2/(pnorm(power)^(-1) - (qnorm(1-0.025)))^2)
        
        Ik = If2*Rp(k = (k[who] - 1), alpha = alpha, beta = power, method = method)

        zl =  t_hat/(sd(diff)/sqrt(length(diff)))
        Il = (length(diff)/(k[who]))*Ik

        B = 1 - pnorm(
          (-zl*sqrt(Il) + (1.96 + abs_delta*sqrt(Ik))*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
        )
        
        C = pnorm(
          (-zl*sqrt(Il) - (1.96 + abs_delta*sqrt(Ik))*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
        )
        
        A =  pnorm(
          (-zl*sqrt(Il) + (-abs_delta*sqrt(Ik) + 1.96)*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
        ) - pnorm(
          (-zl*sqrt(Il) + (abs_delta*sqrt(Ik) - 1.96)*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
        )
        
        
        ci_curt = c(lower, upper)
        
        ci[[who]] = cbind(ci[[who]], ci_curt)
        if((ci_curt[1] >= (-1)*abs_delta) && (ci_curt[2] <= abs_delta)){
          interm[[who]] = cbind(interm[[who]], 0)
          stop_at[who] = which(int_time[i] == c_time[[who]])
          
        }else if(ci_curt[2] <= (-1)*abs_delta){
          interm[[who]] =  cbind(interm[[who]], -1)
          stop_at[who] = which(int_time[i] == c_time[[who]])
          
        }else if(ci_curt[1] >= abs_delta){
          interm[[who]] =  cbind(interm[[who]], 1)
          stop_at[who] = which(int_time[i] == c_time[[who]])
          
        }else{
          if((A < gamma) & (B < gamma) & (C < gamma)){
            interm[[who]] =  cbind(interm[[who]], 100)
            stop_at[who] = which(int_time[i] == c_time[[who]])
            
          }else{
            interm[[who]] =  cbind(interm[[who]], 4)
            stop_at[who] = which(int_time[i] == c_time[[who]])
            
          }
        }
        
      }
      
    }

    bin = bin1 = 0
    who_nstp = c()
    sum1 = 1
    c_cut = NULL
    
    for(j in 1:n){
      if((sum(interm[[j]]!=4) >= 1) | (length(interm[[j]]) == (kn))){
        bin1 = bin1 + 1
      }

      if((sum(interm[[j]]!=4) >= 1) && (sum(interm[[j]]==100) < 1)){
        bin = bin + 1
      }else{
        who_nstp[sum1] = j
        sum1 = sum1 + 1
        
      }
    }
    A = B = C = 0
    for(z in 1:n){
      if(sum(interm[[z]]!=4) >= 1 ){
        place = which(interm[[z]] != 4)
        if(interm[[z]][place] == 0){
          A = A + 1
        }else if(interm[[z]][place] == 1){
          
          if(any(class.save.in == -1)){
            idx = which(class.save.in == -1)
            if(z %in% idx){ 
              C = C + 1
              interm[[z]][place] = -1
            }else{
              B = B + 1
            }
          }else{
            B = B + 1
          }
          
        }else if(interm[[z]][place] == -1){
          C = C + 1
        }
      }
    }
    
    if(bin >= pstar){
      
      replace =  bin1 - pstar + 1

      gamma_glb = w
      lower_A = A
      lower_B = B
      lower_C = C
      
      
      sum_g = sum_g + 1
      
      if( lower_A > gamma_glb ){
        interm_g = 0
        find_id = which(sort(unique(who_active)) %in% c(1:n))
        if(length(find_id) == n){
          who_na = NA
        }else{
          who_na = c(1:n)[-(which(sort(unique(who_active)) %in% c(1:n)))]
        }
  
        SUCCESS_G = T
        
        break
      }else if(lower_B > gamma_glb){
        interm_g = 1
        find_id = which(sort(unique(who_active)) %in% c(1:n))
        if(length(find_id) == n){
          who_na = NA
        }else{
          who_na = c(1:n)[-(which(sort(unique(who_active)) %in% c(1:n)))]
        }
        SUCCESS_G = T
        break
      }else if(lower_C > gamma_glb ){
        interm_g = -1
        find_id = which(sort(unique(who_active)) %in% c(1:n))
        if(length(find_id) == n){
          who_na = NA
        }else{
          who_na = c(1:n)[-(which(sort(unique(who_active)) %in% c(1:n)))]
        }
        SUCCESS_G = T
        break
      }else if((bin>n_pre) &(lower_C <= gamma_glb)&(lower_B <= gamma_glb)&(lower_A <= gamma_glb)){
        interm_g = 4 
        break
      }else{
        interm_g = 4
      }
      
    }else{
      interm_g = 999
    }
    
    
  }
  
  
  if(SUCCESS_G == T){
    int_time = NULL
    c_time_new = list()
    sum = 1
    for(j in unique(who_active)){
      int_time = c(int_time, c_time[[j]])
      c_time_new[[sum]] = c_time[[j]]
      sum = sum + 1
    }
    
    int_time = sort(unique(int_time))
    
    t_hat = c()
    ci = as.list(rep(0, n))
    t_stat = c()
    
    ind_cum = matrix(0, nrow = n, ncol = max(k))
    
    interm = as.list(rep(4, n))
    stop_at = rep(0, n)
    sum_g = 1
    who_active = NULL
    for(i in 1:length(unique(int_time))){
      ck = cumsum(k)
      ppl = which(int_time[i] == unlist(c_time_new))
      who_all = c()
      for(j in 1:length(ppl)){
        who_all[j] = which(ck >= ppl[j])[1]
      }

      for(j in 1:length(ppl)){
        who = who_all[j]
        
        if(all(interm[[who]] == 4)){
          prev = length(which(ind_cum[who,] == 1))

          ind_cum[who,] = c(rep(1, (prev + 1)), rep(0, (ncol(ind_cum) - prev - 1)) ) 
        }else{
          prev = 0
        }

        if((prev >= 1) && all(interm[[who]] == 4) && (prev < k[who])){
          
          who_active = c(who_active, who)
          
          
          diff = d[[who]][1:(prev + 1)]
          t_hat = sum(diff)/length(diff)
          
          lower = t_hat - z_value[prev]*(sd(diff)/sqrt(length(diff)))
          upper = t_hat + z_value[prev]*(sd(diff)/sqrt(length(diff)))
          
          If2 = 1/((abs_delta - (t_hat))^2/(pnorm(power)^(-1) - (qnorm(1-0.025)))^2)
          
          Ik = If2*Rp(k = (k[who] - 1), alpha = alpha, beta = power, method = method)

          zl =  t_hat/(sd(diff)/sqrt(length(diff)))
          Il = (length(diff)/(k[who]))*Ik

          B = 1 - pnorm(
            (-zl*sqrt(Il) + (1.96 + abs_delta*sqrt(Ik))*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
          )
          
          C = pnorm(
            (-zl*sqrt(Il) - (1.96 + abs_delta*sqrt(Ik))*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
          )
          
          A =  pnorm(
            (-zl*sqrt(Il) + (-abs_delta*sqrt(Ik) + 1.96)*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
          ) - pnorm(
            (-zl*sqrt(Il) + (abs_delta*sqrt(Ik) - 1.96)*sqrt(Ik) - (Ik - Il)*t_hat)/sqrt(Ik-Il)
          )
          
          
          ci_curt = c(lower, upper)
          
          ci[[who]] = cbind(ci[[who]], ci_curt)
          if((ci_curt[1] >= (-1)*abs_delta) && (ci_curt[2] <= abs_delta)){
            interm[[who]] = cbind(interm[[who]], 0)
            stop_at[who] = which(int_time[i] == c_time_new[[who]])
            
          }else if(ci_curt[2] <= (-1)*abs_delta){
            interm[[who]] =  cbind(interm[[who]], -1)
            stop_at[who] = which(int_time[i] == c_time_new[[who]])
            
          }else if(ci_curt[1] >= abs_delta){
            interm[[who]] =  cbind(interm[[who]], 1)
            stop_at[who] = which(int_time[i] == c_time_new[[who]])
            
          }else{
            if((A < gamma) & (B < gamma) & (C < gamma)){
              interm[[who]] =  cbind(interm[[who]], 100)
              stop_at[who] = which(int_time[i] == c_time_new[[who]])
              
            }else{
              interm[[who]] =  cbind(interm[[who]], 4)
              stop_at[who] = which(int_time[i] == c_time_new[[who]])
              
            }
          }
          
        }
        
      }

    }
    
    
  }
  
  
  out = NULL
  for(j in 1:n){
    ones = sum(ind_cum[j,] == 1)
    if(length(interm[[j]]) != ones){
      out = c(out, j)
    }
  }
  
  stop_for = c()
  for(i in 1:n){
    place = which(interm[[i]] != 4)
    if(length(place)>=1){
      stop_for[i] = interm[[i]][place][1]
      
      if(any(class.save.in == -1)){
        idx = which(class.save.in == -1)
        if((i %in% idx) & (stop_for[i]==1)){ 
          stop_for[i] = -1
        }else{
          stop_for[i] = interm[[i]][place][1]
        }
      }
      
    }else{
      stop_for[i] = 4
    }
  }
  
  if(any(!is.na(who_na))){
    for(i in 1:length(who_na)){
      stop_for[who_na[i]] = 999
      
    }
  }
  
  
  if(interm_g == 999){
    bin = NA
  }

  return(list(stop_for, interm_g, stop_at, bin))
  
}

Rp = function(k, alpha, beta, method){
  if(method == "Pocock"){
    if(beta == 0.8){
      out = c(1, 1.11, 1.166, 1.202, 1.229, 1.249, 1.265, 1.279, 1.291, 1.301)
    }
    if(beta == 0.9){
      out = c(1, 1.1, 1.151, 1.183, 1.207, 1.225, 1.239, 1.252, 1.262, 1.271)
    }
  }
  if(method == "OBF"){
    if(beta == 0.8){
      #1 - 10
      out = c(1, 1.008, 1.017, 1.024, 1.028, 1.032, 1.035, 1.037, 1.038, 1.04,
              #11, 12, 15, 20
              1.041, 1.042, 1.045, 1.047)
    }      
    if(beta == 0.9){
      out = c(1, 1.007, 1.016, 1.022, 1.026, 1.030, 1.032, 1.034, 1.036, 1.037,
              #11, 12, 15, 20
              1.039, 1.040, 1.042, 1.045)
    }
  }
  if(k <= 12){
    return(out[k])
  }else if(k>12 & k <15){ # approximation
    return((out[13]+out[12])/2) 
  }else if(k == 15){ # approximation
    return(out[13]) 
  }else if(k > 15 & k < 100){# approximation
    return(out[14]) 
  }else if(k >=100){# approximation
    return(1.06) 
  }
}
