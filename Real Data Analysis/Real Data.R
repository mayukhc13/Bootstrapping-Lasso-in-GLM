## The following implementation is on the data available in https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Coimbra
## We applied our method to get important factors that involve in detection of Breast Cancer.
### The number of observation n is 116
## The response variable is binary. In original data set by 1, denote "healthy control(or absence)"
## By 2, denote "patients(or presence of cancer)".
## For ease of notations we adjusted 1 and 2 respectively by 0 and 1 in the simulation. 
## There are p=9 predictor variables as in the data set.
## Here we obtain Lasso estimate of all these 9 factors and 90% Bootstrap Percentile Intervals(Both,Right,Left Sided)
### Finally based on these intervals and lasso estimates, we capture the important factors weighing to the presence of Breast Cancer
## This coincides with the findings of Patricio et al.(2018)




n=116
#load the data
library(datasets)
data(dataR2)
outcome = dataR2$y

### Steps 22 upto step 67 is optional to our simulation purpose. These were done basically to check whether the data set satisfy Logistic regression assumptions.
library(tidyverse)
library(broom)
library(mlbench)
theme_set(theme_classic())

#fit the logistic regression model
model = glm(dataR2$y~.,data = dataR2, family = binomial(link = "logit"))
model
#predict the probability of presence of breast cancer
probabilities = predict(model,type = "response")
probabilities
predicted.classes = ifelse(probabilities>0.5, "pos", "neg")
predicted.classes
head(predicted.classes)

# logistic regression diagnostics 
#(check linear relation between continuous predictors and logit of the outcomes by scatter plot)
#Bind the logit and tidying data for plot
mydata = dataR2 %>% 
  mutate(logit = log(probabilities/(1-probabilities))) %>%
   gather(key = "predictors", value = "predictor.value", -logit)
mydata

#create the scatter plot
#Smooth scatter plot shows the variables Adiponectin,Glucose,MCP.1,Resistin are all quite linerly associated
#However Age, BMI, HOMA, Insulin,Leptin are not linear and might need some transformation
ggplot(mydata, aes(logit, predictor.value)) + geom_point(size=0.5, alpha = 0.8) + 
  geom_smooth(method = "loess") + theme_bw() + facet_wrap(~predictors, scales = "free_y")

#check for influential variables by cook's distance with top 3 largest values
#note that not all outliers are influential variables
plot(model, which = 4, id.n = 3)

## to check for influential observation, we check absolute standardised residual above 3 or not
##extract model results, check for potential outliers, and plot standardised residuals
model.data = augment(model) %>% mutate(index = 1:n())
model.data %>% top_n(3, .cooksd)
ggplot(model.data, aes(index, .std.resid))+ geom_point(aes(color= outcome), alpha = 0.5) + theme_bw()

model.data %>% filter(abs(.std.resid)>3)
#There's no influential variables 

## check for multicollinearity
library(car)
car::vif(model)
## indeed Insulin and HOMA need caution


#### Here starts our main steps of the simulation.
dataf = as.data.frame(dataR2)
dataf
dataf_matrix = as.matrix(dataf)
dataf_matrix
response_y = dataf_matrix[,colnames(dataf_matrix)=="y"]
response_y
design_matrix = dataf_matrix[,colnames(dataf_matrix)!="y"]
design_matrix



### now to fit lasso
library(glmnet)
cv= cv.glmnet(x=design_matrix,y= dataf$y, family= binomial(link = "logit"), type.measure = "class", intercept = F)
cv
lambda_opt = cv$lambda.min
lambda_opt #optimal lambda

fit_model = glmnet(design_matrix,dataf$y,family = binomial(link = "logit"),alpha = 1,intercept = FALSE, lambda=lambda_opt,  penalty.factor = c(rep(1,9)))
fit_model

beta_hat = coef(fit_model)
beta_hat
beta_hat1 = as.array(beta_hat)
beta_hat2 = data.frame(beta_hat1)
beta_hat_final = beta_hat2$s0
beta_hat_final
beta_hat_final3 = beta_hat_final[-1]
beta_hat_final3 #lasso estimate of true parameter

a_n = n^(-(1/3)) #thresholding
beta_tilde = vector()
for(i in 1:9)
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
beta_tilde #thresholded lasso estimate
P_tilde = vector()
for(i in 1:n)
{
  P_tilde[i]= exp(design_matrix[i,]%*%beta_tilde)/(1+exp(design_matrix[i,]%*%beta_tilde))
}
P_tilde

B=150 #bootstrap iterations to be made
beta_bootstrap = matrix(NA,B,9) # collect B copies of bootstrap estimate of lasso
for(c in 1:B)
{
  G_boot = matrix(NA,1,n)
  for (j in 1:n)
  {
    G_boot[1,j]= rexp(1,1)
  }
  
  G_boot
  
  library(CVXR)
  v = CVXR::Variable(9) ##declaration: bootstrap estimate wrt which objective fn is optimized
  
  penalty_boot = lambda_opt*(p_norm(v,1)) ## penalty term in the objective fn
  zz= t(response_y*G_boot) ## a col vector consisting of component wise product of y and G_boot
  l1= -sum(zz*(design_matrix%*%v)) ## * does component wise product and return a vector, %*% returns actually a dot product
  zz1= as.matrix(response_y-P_tilde) ## col vector y-Ptilde
  l2= sum(zz1*(design_matrix%*%v))
  zz2= t(G_boot) ## col vector G_boot
  l3 = sum(zz2*(logistic(design_matrix%*%v))) ## logistic(z) is an atom for log(1+exp z) in cvxr and sum() does the component wise sum of that vector
  obj_boot = penalty_boot+l1+l2+l3  ## objective function
  prob_boot = Problem(Minimize(obj_boot)) ## minimisation of obj fn
  resultB = solve(prob_boot)
  resultB
  
  
  final_boot = resultB$getValue(v)
  final_boot1 = as.matrix(final_boot)
  final_boot2 = t(final_boot1)
  final_boot2
  beta_bootstrap[c,]= final_boot2 ## bootstrap estimator of lasso at each stage
}
beta_bootstrap

unit_vector = as.matrix(rep(1,B))
beta_tilde_matrix = unit_vector%*%beta_tilde
beta_tilde_matrix
tau_boot= matrix(NA,B,9)
for(i in 1:B){
  for(j in 1:9){
    tau_boot[i,j]=sqrt(n)*(beta_bootstrap[i,j]-beta_tilde_matrix[i,j])
  }
}
tau_boot ## B copies of centered bootstrap-beta tilde distribution scale with sqrt n

sorted_tau_boot = apply(tau_boot,2,sort)
sorted_tau_boot ## independent col wise sorting of tau boot for getting critical value

alpha = 0.10 #level of significance
lower_index = as.integer(B*(alpha/2)) #both sided
index_right = as.integer(B*alpha) #right sided
index_left = 1+ as.integer(B*(1-alpha)) #left sided
upper_index = 1+ as.integer(B*(1-(alpha/2))) #both sided
lower_index
upper_index
index_left
index_right

critical_value_matrix = matrix(NA,2,9)
for(j in 1:9)
{
  critical_value_matrix[1,j]= -(sorted_tau_boot[upper_index,j])/sqrt(n)
  critical_value_matrix[2,j]= -(sorted_tau_boot[lower_index,j])/sqrt(n)
}
critical_value_matrix

critical_value_matrix_right = matrix(NA,1,9)
right_confidence_limit_matrix = matrix(NA,1,9)
for(j in 1:9){
  critical_value_matrix_right[1,j] = -(sorted_tau_boot[index_right,j])/sqrt(n)
  right_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix_right[1,j]
}
critical_value_matrix_right
right_confidence_limit_matrix

critical_value_matrix_left = matrix(NA,1,9)
left_confidence_limit_matrix = matrix(NA,1,9)
for(j in 1:9){
  critical_value_matrix_left[1,j] = -(sorted_tau_boot[index_left,j])/sqrt(n)
  left_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix_left[1,j]
}
critical_value_matrix_left
left_confidence_limit_matrix

lower_confidence_limit_matrix = matrix(NA,1,9)
upper_confidence_limit_matrix = matrix(NA,1,9)
for(j in 1:9)
{
  lower_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix[1,j]
  upper_confidence_limit_matrix[1,j]= beta_hat_final3[j]+critical_value_matrix[2,j]
}
lower_confidence_limit_matrix
upper_confidence_limit_matrix

