### read data and process #######################

set.seed(2023)

df = read.csv2(file = 'Use2.csv')
df = df[df$HUC != 0,] #delete missing HUCs

normalize <- function(x) (x - mean(x)) / sd(x)

# different to opsomer !
anc <- df$ANC / 1000 # in equivalents per liter
elevation <- df$GRID.ELEV.M 
huc <- factor(df$HUC, levels = unique(df$HUC))
lat = (df$LAT.DD.1)
long = (df$LONG.DD)

n = dim(df)[1]
ni = table(huc) # huc groups

### get spline basis #######################

# do as in Opsomer

Cfun = function(x) sum(x**2) * log(sum(x**2)**0.5)
K = 80

require('fields')

lat_grid = seq(min(lat), max(lat), length = 100)
long_grid = seq(min(long), max(long), length = 100)

candidates = make.surface.grid(list(lat_grid, long_grid))
cand = candidates
tol = 0.01
for (i in 1:length(lat_grid) ** 2) {
  res = c()
  for (k in 1:length(lat)) res[k] = sum((candidates[i,] - c(lat[k], long[k]))**2)
  if (min(res) > tol) cand[i,] = c(NA,NA)
}
cand = cand[!is.na(cand[,1]), ]
out = cover.design(cand, K)

# plot knots of spline basis (cf Opsomer)
plot(lat,long, col = 8, pch = 16, lwd = 0.2)
# points(candidates[,1], candidates[,2], col = 4, lwd = 0.1) # all candidate knots ...
points(cand[,1], cand[,2], col = 3, lwd = 0.1) # ... in neighbourhood
points(out[,1], out[,2], col = 2, lwd = 2) # selected knots

# construct random effect matrix as in Opsomer
Z1 = matrix(NA, n, K)
for(i in 1:n){
  for(j in 1:K){
    Z1[i,j] = Cfun(c(lat[i], long[i]) - out[j,])
  }
}
Z2 = matrix(NA, K, K)
for(i in 1:K){
  for(j in 1:K){
    Z2[i,j] = Cfun(out[i,] - out[j,])
  }
}
diag(Z2) = 0
E = svd(Z2)
OMEGA = E$u %*% diag( E$d ** -0.5) %*% t(E$v)
Z = Z1 %*% OMEGA


### set up model equation #######################
# normalize <- function(x)  (x - mean(x)) / sd(x)

y = anc
X = cbind(1, elevation); colnames(X) = c('intercept', 'elevation')

Hv = as.matrix(bdiag(mapply(function(n) matrix(1,n,n), ni))) # HUC RE
Hk = Z %*% t(Z) # spatial RE


### get reml estimates #######################

# compute estimates directly, nlme / lmer4 in R doesnt support these covariance matrices 
tr <- function(X) sum(diag(X))
quadform <- function(A, x) {
  return(as.matrix(Matrix::crossprod(Matrix::crossprod(A, x), x))) # = x^t A x
}

delta <- rep(var(anc), 3)
delta_old <- rep(0,3)

while (sum((delta - delta_old)^2) / norm(delta) > 1e-1) { 
  delta_old <- delta
  
  V = delta[1] * Hv + delta[2] * Hk + delta[3] * diag(1,n)
  Vi = solve(V)
  ViX <- Vi %*% X
  P = Vi - ViX %*% solve( t(X) %*% ViX ) %*% t(ViX)
  
  PHvP <- quadform(Hv, P) # P %*% Hv %*% P 
  PHkP <- quadform(Hk, P) # P %*% Hk %*% P 
  PP <- P %*% P 
  
  
  score <- 1/2 * c(
    -tr(P %*% Hv) + quadform(PHvP, y),
    -tr(P %*% Hk) + quadform(PHkP, y),
    -tr(P) + quadform(PP, y)
  )
  
  Fisher <- 1/2 * matrix(c(
    tr(PHvP %*% Hv), 
    tr(PHvP %*% Hk), 
    tr(PHvP), 
    0, 
    tr(PHkP %*% Hk), 
    tr(PHkP), 
    0, 0, 
    tr(P %*% P)
  ), 3, 3); Fisher <- Fisher + t(Fisher) - diag(diag(Fisher))
  
  
  delta <- delta + solve(Fisher) %*% score
}

# 
delta # = (26201.24 314509.16  57380.79), Opsomer: (5069.44 133736.49  32220.25)
max(delta) / min(delta) # 12.00


# V = delta[1] * Hv + delta[2] * Hk + delta[3] * diag(1,n)
# Vi = solve(V)
# ViX <- Vi %*% X
# P = Vi - ViX %*% solve( t(X) %*% ViX ) %*% t(ViX)
# sum(log(eigen(V)$values)) + sum(log(eigen(t(X) %*% ViX)$values)) + y %*% P %*% y


### get wls estimates #######################

V = delta[1] * Hv + delta[2] * Hk + delta[3] * diag(1,n); Vi = solve(V)
beta_wls = solve( t(X) %*% Vi %*% X) %*% t(X) %*% Vi %*% y

beta_wls # = (257  -1.16), Opsomer: (228.6, -0.814)
diag(solve( t(X) %*% Vi %*% X))**.5 #sd for beta_wls



### do lasso  #######################
require('glmnet')

# transform model to linear model by left-multiplication of V^(-1/2)

E = eigen(V)
Vis = E$vec %*% diag(E$val ** -0.5) %*% t(E$vec)
yi = Vis %*% y
Xi = Vis %*% X

fit <- cv.glmnet(Xi, yi, intercept = F, nfolds = 551, grouped = F)
lambda <- fit$lambda.min # ATTENTION: is lambda / n in our notation 
beta_lasso <- coef(fit$glmnet.fit, s = fit$lambda.min)[2:3]

lambda <- lambda * n # to bring in standard notation  
lambda # 14.7
sqrt(n) / 2 # 11.7 / close enough
beta_wls # = 257.688408  -1.163266
beta_lasso # =  0.000000 -1.068735



### make inference #####################

alpha <- 0.01

# get non-centrality parameter
C = as.matrix(t(X) %*% Vi %*% X / n); Ci = solve(C) # awfully confitioned... 

LAMd = list( c(1,1)*lambda*n**(-1/2) ,  c(-1,1)*lambda*n**(-1/2) )
tau = max( mapply(function(a) t(a) %*% Ci %*% a, LAMd) )
tau

# ratio of volumes
# volume of |Ax|^2 <= 1 is det(A**-1) * vol unit ball
# in the ratio of volumes, det( C^(-1/2) ) / n cancels
qchisq(1-alpha, 2, ncp = tau) / qchisq(1-alpha, 2) # gives warnings, non-centrality parameter is too large 


# just the dimension for elevation
LAMd.slice = list(c(lambda*n**(-0.5)), c(-lambda*n**(-0.5)))
tau = max( mapply(function(a) t(a)%*%Ci[2,2]%*%a, LAMd.slice) )

# ratio of volumes
qchisq(1-alpha, 1, ncp = tau)**.5 / qchisq(1-alpha, 1)**.5


### p values - difficult to interpret but to compare to opsomer ################

pchisq(beta_wls[1] ** 2 * C[1,1] * n, df = 1, lower.tail = F) # intercept insignificant
pchisq(beta_wls[2] ** 2 * C[2,2] * n, df = 1, lower.tail = F) # elevation significant

# both vectors together are significant
pchisq(n * t(beta_wls) %*% C %*% beta_wls, df = 2, lower.tail = F)
pchisq(n * t(beta_lasso) %*% C %*% beta_lasso, ncp = tau, df = 2, lower.tail = F)


# find confidence intervals 

CI_lasso = c(
  beta_lasso[2] - (qchisq(1-alpha, 1, ncp = tau) / n / C[2,2])**.5,
  beta_lasso[2] + (qchisq(1-alpha, 1, ncp = tau) / n / C[2,2])**.5
)

CI_wls = c(
  beta_wls[2] - (qchisq(1-alpha, 1) / n / C[2,2])**.5,
  beta_wls[2] + (qchisq(1-alpha, 1) / n / C[2,2])**.5
)
df = data.frame(WLS = c(CI_wls, diff(CI_wls)), 
           LASSO = c(CI_lasso, diff(CI_lasso)))
rownames(df) = c("lower bound", "upper bound", "interval length")
df
