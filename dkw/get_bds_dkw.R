##### bounds generator for DKW simple formulation


get_bds_dkw_simp = function(alpha = 0.05, n)
{
  bds_dkw_simp = data.frame(low = rep(NA, n),
                      up = rep(NA, n))
  
  bds_dkw_simp$low = (1:n)/n - sqrt(( log(2) - log(alpha))/(2*n))
  bds_dkw_simp$up = (1:n)/n + sqrt(( log(2) - log(alpha))/(2*n))
  
  bds_dkw_simp$low[ bds_dkw_simp$low < 0] = 0
  bds_dkw_simp$up[ bds_dkw_simp$up > 1] = 1
  
  return(bds_dkw_simp)
}

lambda_uppbd_func = function(lambda, n, alpha )
{
  x = lambda
  mu = 0.4345
  temp =  -1 + (x + 1/(4*x) + (3*mu)/x + (3*mu)/(2*x^2) )/sqrt(n)
  temp = temp - (mu/2)*(4 +  (1/x^2) )/n 
  temp = temp + (mu/2)*(4*x + (1/x) )*n^(-1.5)
  temp = (temp*2*x)/(3*sqrt(n)) + 1
  temp = temp*(exp(-2*x^2))
  temp = temp - alpha/2
  return(temp)
}

get_bds_dkw_complex = function(alpha = 0.05, n)
{
  lam_alpha =  uniroot( lambda_uppbd_func, lower =  1 , upper = 10, 
                        n = n , alpha = alpha, tol = 10^-10  )$root
  
  bds_dkw_complex = data.frame(low = rep(NA, n),
                            up = rep(NA, n))
  
  bds_dkw_complex$low = (1:n)/n - lam_alpha/sqrt(n)
  bds_dkw_complex$up = (1:n)/n  + lam_alpha/sqrt(n)
  
  bds_dkw_complex$low[ bds_dkw_complex$low < 0] = 0
  bds_dkw_complex$up[ bds_dkw_complex$up > 1] = 1
  
  return(bds_dkw_complex)
}

# n = 1000
# plot(x = 1:n/n, y = rep(0,n), type = "l", main = "Comparision of simple vs. complex DKW (dotted),n = 100", ylim = c(-0.05,0.05))
# points(x = 1:n/n, y = get_bds_dkw_simp(n = n)$low - 1:n/n, type = "l", col = "red")
# points(x = 1:n/n, y = get_bds_dkw_simp(n = n)$up- 1:n/n, type = "l", col = "green")
# 
# bds_dw = get_bds_dw(n = n)
# points(x = 1:n/n, y = bds_dw$low - 1:n/n, type = "l", col = "red", lty = 2)
# points(x = 1:n/n, y = bds_dw$up - 1:n/n, type = "l", col = "green", lty = 2)
# 
# bds_dkw_c = get_bds_dkw_complex(n = n)
# points(x = 1:n/n, y = bds_dkw_c$low- 1:n/n, type = "l", col = "dark red", lty = 3)
# points(x = 1:n/n, y = bds_dkw_c$up-1:n/n, type = "l", col = "dark green", lty = 3)
#  
# 
