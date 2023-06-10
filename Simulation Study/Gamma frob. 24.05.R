#### h2o package can be installed from https://cran.r-project.org/bin/windows/contrib/4.3/h2o_3.40.0.4.zip or 
### from https://cran.r-project.org/src/contrib/h2o_3.40.0.4.tar.gz


## initiation of h2o local host
library(h2o)
h2o.init(ip="localhost",port=54321,max_mem_size = "4g",nthreads = -1)

library(CVXR) # for convex optimization
# Step 1:

n=500 ## SAMPLE SIZE
alpha_shape = 1 ## shape parameter for gamma distribution which is known
p=7   ## number of parameter components
p_0=4   #non zero component

beta=vector()
for(j in 1:p)
{
  if(j<=p_0)
  {
    beta[j]= 0.5*j*((-1)^j)
  }
  else{
    beta[j]=0
  }
}
beta  #true parameter choice


# Step 2: 

library(MASS)   #for multivariate normal sampling
mu= rep(0,p)  ## mean vector with all zero components
mu
cov_matrix = matrix(NA,p,p)
for(i in 1:p)
{
  for(j in 1:p)
  {
    if(i==j)
    {
      cov_matrix[i,j] = 1
    }
    else
    {
      cov_matrix[i,j]= 0.3^(abs(i-j))
    }
  }
}
cov_matrix  #covariance matrix

# Step 3:
x = mvrnorm(n,mu,cov_matrix)  # matrix, drawn from multivariate normal with above mean and covariance structure
x
uv= as.matrix(rep(1,n))
uv
x_int1=cbind(uv,x)
x_int1
x_int=as.matrix(x_int1)
x_int  # matrix with first column all with 1

z = scale(x,T,F) ## design matrix without intercept column, is just centered which will be used as actual design matrix throuhout.
z


intercept = 0 # true intercept component is chosen deliberately as 0.

P = vector()
for (i in 1:n)
{
  P[i]= (exp(intercept+(z[i,]%*%beta)))/alpha_shape ## from log link relation and taking unknown parameter as scale parameter.
  if(is.nan(P[i])==TRUE)
  {
    i=i-1
  }
  else next
  
}
P # scale parameter vector of n components with which n gamma response variables will be drawn with known shape parameter 1.


g_inv_der_lasso = vector() ## mean vector of n gamma response variables
for(i in 1:n)
{
  g_inv_der_lasso[i] = exp(z[i,]%*%beta)
  if(is.nan(g_inv_der_lasso[i])==TRUE)
  {
    i=i-1
  }
  else next
}
g_inv_der_lasso

# Step 4:
simulation = 500 # entire iteration of data to be done
width_matrix = matrix(NA, simulation,p) #both sided CI matrix width, (upper confidence limit-lower confidence limit) for each component
count_coverage = matrix(NA,simulation,p) # both sided empirical coverage count matrix, count is 1 if true parameter value falls within both sided limit
count_coverage_right = matrix(NA,simulation,p) #right sided empirical coverage count matrix
count_coverage_left = matrix(NA,simulation,p)   # left sided empirical coverage count matrix
count_norm_matrix = matrix(NA,simulation,1)  # count coverage in terms of euclidean norm
LAMBDA = matrix(NA,simulation,1) ## A column vector of all 500 optimal lambda value from cross validation in gamma regression
for(m in 1:simulation)
{
  y= vector()
  for(i in 1:n)
  {
    y[i]= rgamma(1,alpha_shape,scale = P[i])
    if(is.nan(y[i])==TRUE)
    {
      i=i-1
    }
    else next
  }
  y  ## Gamma response vector with known shape parameter 1 and scale parameter P
  
  
  # Step 5:
  
  y_final = data.frame(y)
  x_final = data.frame(z)
  gamma_data = cbind(y,z)
  gamma_data_final = data.frame(gamma_data)
  
  
  library(h2o)
  data.gamma = as.h2o(gamma_data_final,destination_frame = "data.gamma")
  CV = h2o.glm(y=1,x=2:8,training_frame = data.gamma,nfolds = 10,standardize = FALSE,family = "gamma",link = "log",alpha = 1)
  CV #cross validation for optimal lambda in GAMMA regression.
  lambda_opt = h2o.getLambdaBest(CV)
  lambda_opt ## optimal lambda to be used later in Bootstrap
  LAMBDA[m,1] = lambda_opt
  if(is.nan(LAMBDA[m,1])== TRUE || is.na(LAMBDA[m,1])== TRUE)
  {
    m = m-1
    next
  }


  ## Next we get LASSO estimate of original parameter through CVXR using the form of Negative log likelihood in GLM 
  
  library(CVXR)
  v1 = CVXR::Variable(p)
  penalty_lasso = LAMBDA[m,1]*(p_norm(v1,1))
  l_1 = alpha_shape*(sum(z%*%v1))
  z_1 = as.matrix(y)
  l_2 = alpha_shape*sum((z_1*(exp(-(z%*%v1)))))
  obj_lasso = l_1+l_2+penalty_lasso
  prob_lasso = Problem(Minimize(obj_lasso))
  result_lasso = solve(prob_lasso)
  result_lasso
  final_lasso=  as.function(result_lasso$getValue , list("v1"))
  final_lasso1 = final_lasso(v1)
  final_lasso1 = final_lasso1[!is.nan(final_lasso1)]
  final_lasso1
  final_lasso2= as.matrix(final_lasso1)
  final_lasso2
  beta_hat_final3=t(final_lasso2)
  beta_hat_final3  ### Lasso estimate
  
  
  stat3= sqrt(n)*(beta_hat_final3-beta) #centered distribution of (lasso estimate-original parameter)*sqrt(n)
  stat3
  stat4= as.matrix(stat3)
  stat4
  stat3_norm = norm(stat4,type = "f") # euclidean norm of above to check count norm CI to get Table 3 in simulation study
  stat3_norm
  
  # Next we do, thresholding of Lasso estimator that is to be used as centering factor later in Bootstrap case as in section 4.2
  a_n = n^(-(1/3)) ## sequence value a_n is chosen as in sec 4.2
  beta_tilde = vector()
  for(i in 1:p)
  {
    if(abs(beta_hat_final3[i])>a_n)
    {
      beta_tilde[i]= beta_hat_final3[i]
    }
    else
    {
      beta_tilde[i]=0
    }
  }
  beta_tilde  #thresholded lasso estimate
  
  
  
  g_inv_der = vector()
  for(i in 1:n)
  {
    g_inv_der[i] = exp(z[i,]%*%beta_tilde)
    if(is.nan(g_inv_der[i])==TRUE)
    {
      i=i-1
    }
    else next
  }
  g_inv_der #mean vector of response variables at beta tilde
  
  
  # Step 7:
  
  B=500 #bootstrap iterations to be made
  beta_bootstrap = matrix(NA,B,p) # collect B copies of bootstrap estimate of lasso estimator as in (2.3) in section 2
  for(c in 1:B)
  {
    G_boot = matrix(NA,1,n)
    for (j in 1:n)
    {
      G_boot[1,j]= rexp(1,1)  ## extract n iid exp(1) r.v for each stage optimization in objective function.
    }
    
    G_boot  #perturbation quantities from known density Exp(1), n copies at each stage, repeat B times
    
    ## Next we do Bootstrap convex optimization of (2.3) though CVXR function
    library(CVXR)
    v2 = CVXR::Variable(p) ##declaration: bootstrap estimate wrt which objective function is optimized
    
    penalty_boot = LAMBDA[m,1]*(p_norm(v2,1)) ## penalty term in the objective function through p_norm syntax.
    zz= t(y*G_boot)  ## a col vector consisting of component wise product of y and G_boot. Also * does component wise product and return a vector
    l11 = alpha_shape*sum(zz*(exp(-(z%*%v2)))) ##  %*% returns actually a scalar product and l11 is the first component in the first term as in (2.3) 
    zz2= t(G_boot) ## col vector G_boot
    l22 = alpha_shape*sum(zz2*(z%*%v2)) ##l22 is the second component in the first term as in (2.3)
    zz3 =  as.matrix(y-g_inv_der)  ## col vector (y-mean vector of response variables at beta tilde)
    zz4 = exp(-(z%*%beta_tilde)) 
    l33 = (alpha_shape)*sum(zz3*zz4*(z%*%v2)) ## l33 is the second term as in (2.3)
    obj_boot = l11+l22+l33+penalty_boot   ## objective function
    prob_boot = Problem(Minimize(obj_boot))  ## minimisation of obj fn
    result_boot = solve(prob_boot)
    result_boot
    
    final_boot=  as.function(result_boot$getValue , list("v2"))
    final_boot1 = final_boot(v2)
    final_boot1 = final_boot1[!is.nan(final_boot1)]
    final_boot1
    final_boot2=as.matrix(final_boot1)
    final_boot2
    final_boot3 = t(final_boot2)
    final_boot3  ## Minimizer Bootstrap estimator 
    beta_bootstrap[c,]= final_boot3  ## At c-th stage we store that bootstrap estimator of p components
    
  }
  beta_bootstrap  ##This matrix is important to construct Bootstrap Percentile Intervals
  
  unit_vector = as.matrix(rep(1,B))
  beta_tilde_matrix = unit_vector%*%beta_tilde
  beta_tilde_matrix
  
  ## Tau boot is the matrix containing B copies of centered estimator vector sqrt(n)*(Bootstrap estimator - thresholded lasso estimator)
  tau_boot= matrix(NA,B,p)
  for(i in 1:B){
    for(j in 1:p){
      tau_boot[i,j]=sqrt(n)*(beta_bootstrap[i,j]-beta_tilde_matrix[i,j])
    }
  }
  tau_boot ## B copies of centered bootstrap-beta tilde distribution scale with sqrt n
  
  
  
  sorted_tau_boot = apply(tau_boot,2,sort)
  sorted_tau_boot ## independent column wise sorting of tau boot for getting critical value in percentile intervals
  
  # Step 9:
  
  alpha = 0.10 #level of significance
  lower_index = as.integer(B*(alpha/2)) #both sided
  upper_index = 1+ as.integer(B*(1-(alpha/2))) #both sided
  index_right = as.integer(B*alpha) #right sided
  index_left = 1+ as.integer(B*(1-alpha)) #left sided
  
  lower_index
  upper_index
  index_left
  index_right
  
  tau_boot_norm = matrix(NA,B,1) #euclidean norm of each row of tau boot
  for(i in 1:B){
    tau_boot_norm[i,1]= norm(as.matrix(tau_boot[i,]),type = "f")
  }
  tau_boot_norm
  
  sorted_tau_boot_norm = apply(tau_boot_norm,2,sort) #sorted
  sorted_tau_boot_norm
  
  if(stat3_norm<=sorted_tau_boot_norm[index_left,1]){
    count_norm = 1
  }else
  {
    count_norm = 0
    
  }
  count_norm_matrix[m,1]=count_norm #count if stat3norm is less than sorted tauboot norm at 0.90 level
  
  critical_value_matrix = matrix(NA,2,p)
  for(j in 1:p)
  {
    critical_value_matrix[1,j]= -(sorted_tau_boot[upper_index,j])/sqrt(n)
    critical_value_matrix[2,j]= -(sorted_tau_boot[lower_index,j])/sqrt(n)
  }
  critical_value_matrix #both sided
  
  critical_value_matrix_right = matrix(NA,1,p)
  right_confidence_limit_matrix = matrix(NA,1,p)
  for(j in 1:p){
    critical_value_matrix_right[1,j] = -(sorted_tau_boot[index_right,j])/sqrt(n)
    right_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix_right[1,j]
  }
  critical_value_matrix_right
  right_confidence_limit_matrix #right sided
  
  critical_value_matrix_left = matrix(NA,1,p)
  left_confidence_limit_matrix = matrix(NA,1,p)
  for(j in 1:p){
    critical_value_matrix_left[1,j] = -(sorted_tau_boot[index_left,j])/sqrt(n)
    left_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix_left[1,j]
  }
  critical_value_matrix_left
  left_confidence_limit_matrix #left sided
  
  
  lower_confidence_limit_matrix = matrix(NA,1,p)
  upper_confidence_limit_matrix = matrix(NA,1,p)
  for(j in 1:p)
  {
    lower_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix[1,j]
    upper_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix[2,j]
  }
  lower_confidence_limit_matrix
  upper_confidence_limit_matrix #both sided
  
  
  count_beta_right = matrix(NA,1,p)
  for(j in 1:p){
    if(beta[j]<=right_confidence_limit_matrix[1,j]){
      count_beta_right[1,j]=1
    }
    else{
      count_beta_right[1,j]=0
    }
  }
  count_beta_right
  
  count_beta_left = matrix(NA,1,p)
  for(j in 1:p){
    if(beta[j]>=left_confidence_limit_matrix[1,j]){
      count_beta_left[1,j]=1
    }
    else{
      count_beta_left[1,j]=0
    }
  }
  count_beta_left
  
  count_beta = matrix(NA,1,p)
  for(j in 1:p)
  {
    if(beta[j]>= lower_confidence_limit_matrix[1,j] && beta[j]<= upper_confidence_limit_matrix[1,j])
    {
      count_beta[1,j]= 1
    }
    else
    {
      count_beta[1,j]= 0
    }
  }
  count_beta
  
  
  width_matrix[m,] = upper_confidence_limit_matrix[1,]-lower_confidence_limit_matrix[1,]
  count_coverage[m,] = count_beta[1,]
  count_coverage_right[m,]=count_beta_right[1,]
  count_coverage_left[m,]= count_beta_left[1,]
}

width_matrix
count_coverage
count_coverage_right
count_coverage_left
count_norm_matrix
LAMBDA

average_width_matrix = matrix(NA,1,p)
for(j in 1:p)
{
  average_width_matrix[1,j] = mean(width_matrix[,j])
}
average_width_matrix


coverage_probability_right = matrix(NA,1,p)
for(j in 1:p){
  coverage_probability_right[1,j]= mean(count_coverage_right[,j])
}
coverage_probability_right ## Right sided empirical coverage probability of right sided 90% CI

coverage_probability_left = matrix(NA,1,p)
for(j in 1:p){
  coverage_probability_left[1,j]= mean(count_coverage_left[,j])
}
coverage_probability_left  ## Left sided empirical coverage probability of left sided 90% CI


coverage_probability = matrix(NA, nrow=1, ncol= p)
for(j in 1:p)
{
  coverage_probability[1,j]= mean(count_coverage[,j])
}
coverage_probability ## both sided empirical coverage probability of 90% CI

count_norm_coverage_probability = mean(count_norm_matrix[,1])
count_norm_coverage_probability ## see table 3



