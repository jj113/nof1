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

