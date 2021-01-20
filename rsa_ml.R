library(haven)
library(tidyverse)
library(randomForest)
library(ggthemes)

# dataset based on:
# Generating Skilled Self-Employment in Developing Countries: Experimental Evidence 
# from Uganda (Blattman et al, 2014)
# outcome of interest is profit in last four weeks capped at 99th percentile; 
# treatment is randomised assignment to the Program (intent to treat)
# method is based on:
# Generic Machine Learning Inference on Heterogenous Treatment Effects in Randomized 
# Experiments (Chernozhukov et al, 2017)
# Note: This is merely preliminary EDA and will be refined once more of the data is 
# understood
# Final code will be run with neural nets and other forest based methods to derive 
# Best Linear Predictor

#some data cleaning

x<-read_dta('https://www.dropbox.com/s/yxgigmtcrut9fii/yop_analysis.dta?dl=1') %>% 
  filter(e2==1) %>%  #filter results only from endline survey1
  select(assigned,S_K,S_H,S_P_m,admin_cost_us,
         risk_aversion,female,urban,age,
         wealthindex,bizasset_val_real_p99_e,group_female,
         group_roster_size,avgdisteduc,grp_leader,
         group_age,group_existed,nonag_dummy,zero_hours,
         skilledtrade7da_zero,highskill7da_zero,adl,aggression_n,
         ingroup_dynamic,ingroup_hetero,acto7da_zero,D_1,D_2,D_3,
         D_4,D_5,D_6,D_7,D_8,D_9,D_10,D_11,D_12,D_13,chores7da_zero)

#compute proportion of all NA values
sum(x%>%is.na())/(ncol(x)*nrow(x))

#compute proportion of missing outcome values
sum(x$bizasset_val_real_p99_e%>%is.na())/nrow(x)

#eliminate NA values (assume missingness is random)
df<-x[which(complete.cases(x)),]

#impute outcome values using median (assume missingness is non random)
df<-x %>% mutate(bizasset_val_real_p99_e=ifelse(is.na(bizasset_val_real_p99_e),
                                                median(x$bizasset_val_real_p99_e,na.rm = T),
                                                bizasset_val_real_p99_e))

#impute NA values for other vars
df<-rfImpute(bizasset_val_real_p99_e~.,df,ntree=500,iter=6)

# create empty dataframes to store values
n_split<-100
n_group<-5
blp_coef<-data.frame(B1=1:n_split,SE_B1=1:n_split,
                     P_value_B1=1:n_split,
                     B2=1:n_split,SE_B2=1:n_split,
                     P_value_B2=1:n_split)

gate_coef<-matrix(ncol = n_group*2,nrow = n_split) %>% as.data.frame()
colnames(gate_coef)<-paste(c("G",'SE_G'), rep(1:n_group, each=2), sep = "")

#f<-which(sapply(df, class) == "factor") 
#mcol<-ncol(df[,-f])
mcol<-ncol(df)
gate_mean<-matrix(ncol = mcol*n_group,nrow = n_split)
colnames(gate_mean)<-rep(colnames(df[,]),n_group)

#split data
set.seed(55,sample.kind = 'Rounding')
for(i in 1:n_split){
  
  #randomly split data into main and auxiliary 
  random<-runif(nrow(df))
  main_ind<-which(random>0.5)
  aux_ind<-which(random<0.5)
  aux_df<-df[aux_ind,]
  main_df<-df[main_ind,]
  
  # train data on auxiliary sample
  rftreat<-randomForest(bizasset_val_real_p99_e~., data = (aux_df%>%filter(assigned==1)),
                        ntree=500,nodesize=5)  
  rfbase<-randomForest(bizasset_val_real_p99_e~., data = (aux_df%>%filter(assigned==0)), 
                       ntree=500,nodesize=5)
  
  # predict baseline and treatment outcomes on main sample
  B<-predict(rfbase,main_df)
  treat<-predict(rftreat,main_df)
  
  # specifying regression variables
  S<-treat-B #CATE: what the algorithm predicts is an individual's treatment effect
  ES<-mean(S) # the average predicted treatment effect
  p<-mean(main_df$assigned) #take mean as propensity score
  x<-S-ES #excess CATE: how far one's predicted treatment effect is from the mean
  w<-main_df$assigned-p #weighted treatment var
  
  #derive Best Linear Predictor from main sample
  blp<-lm(bizasset_val_real_p99_e~B+w+I((w*x)),data=cbind(main_df,B,S,x,w))
  blp_coef[i,]<-c(blp$coefficients[3],summary(blp)$coefficients[3:4,c(2,4)][1,],
                  blp$coefficients[4],summary(blp)$coefficients[3:4,c(2,4)][2,])
  
  
  #Group Average Treatment Effect
  qt<-quantile(S,seq(0,1,length.out = n_group+1))
  
  for(k in 1:n_group){  
    G<-ifelse(S>qt[k] & S<qt[k+1],1,0)
    gate<-lm(bizasset_val_real_p99_e~B+I(w*G),data = cbind(main_df,B,S,x,w,G))
    gate_coef[i,c((2*k)-1,2*k)]<-summary(gate)$coefficients[3,c(1,2)]
    
    # data preparation for later
    gate_mean[i,((k*mcol)-(mcol-1)):(k*mcol)]<-apply(main_df[which(G==1),],2,mean)
  }  
}  

# check distribution of coefficient
hist(blp_coef$B2,freq = F)

# obtain median values
# data for each column does not correspond to the same split
apply(blp_coef,2,median) #median for Best Linear Predictor 
apply(gate_coef,2,median) #median for Grouped Average Treatment Effect

# examining heterogeneity

for(k in 1:n_group) { 
  nam <- paste("gate_mean", k, sep = "")
  assign(nam, gate_mean[,((k*mcol)-(mcol-1)):(k*mcol)])
}

het<-matrix(ncol = mcol,nrow = n_group) %>% as.data.frame() %>%
  mutate(gate=1:n_group)
colnames(het)<-c(colnames(df[,]),'gate')

for(k in 1:n_group) { 
  het[k,(1:mcol)]<-apply(gate_mean[,((k*mcol)-(mcol-1)):(k*mcol)],2,median)
}

#choose covariates to plot
sub<-c('S_K','S_H','S_P_m','wealthindex','age','group_age','risk_aversion',
       'aggression_n','highskill7da_zero')
het_sub<-het[,c('gate',sub)] %>%
  gather("covariate", "median", -gate)

ggplot(het_sub, aes(gate, median, color = covariate)) +
  geom_line() +
  facet_wrap(~covariate,scales = "free_y") +
  theme_fivethirtyeight()

