data {
  int<lower=0> n_trees; //Amount of sampled sites
  int<lower=0> n_surveys;//Number of surveys
  int<lower=0> n_sp;//Number levels for the species random effects
  int<lower=0> n_p;//Number levels for the properties random effects
  int<lower=-1,upper=1> y[n_trees,n_surveys]; //Survival of each trees on the two surveys
  vector[n_trees] height; //predictor
  vector[n_trees] ca; 
  vector[n_trees] siteExp;
  vector[n_trees] elev;
  vector[n_trees] fc2019;//predictos
  vector[n_trees] roadP;
  vector[n_trees] propsize;
  vector[n_trees] fc2012;
  int<lower=0> sp[n_trees];
  int<lower=0> propid[n_trees];
  int<lower=0> obs; //number of observed observations from the groundtruthed data
  int<lower=0> field[obs]; //index of trees measured in the field
  int<lower=0> remote[n_trees-obs]; //index of trees measured only using remote sensing
  
}

parameters{
  real delta;
  real a;
  real b;
  real d;
  real e;
  real f;
  real g;
  real spel[n_sp];
  real propl[n_p];
  real mu_spe;
  real mu_prop;
  real gamma_spe;
  real gamma_prop;
  real<lower=0.0001,upper=0.9999> mug;
  real<lower=0,upper=1>  p;   //probability of detecting as dead an alive tree(false negative)(1-p11) classified using visual classification
  real<lower=0,upper=p>  q;   //probability of detencting as alive a dead tree (false positive)classified using visual classification (p10)
}

transformed parameters{
  real<lower=0,upper=1>  psi[n_trees]; //probability survival using visual classification data
  real spe[n_sp];
  real prop[n_p];
  real gamma;
  vector[n_trees] log_lik;
  
  //alpha=log(mua/(1-mua));
  gamma=log(mug/(1-mug));
  //kappa=log(muk/(1-muk));
  
  for(n in 1:n_sp){
    spe[n]=mu_spe+spel[n]*gamma_spe;
    
  }
  for(l in 1:n_p){
    prop[l]=mu_prop+propl[l]*gamma_prop;
  }
  
  
  for (i in 1:n_trees){
    psi[i]=inv_logit(gamma+spe[sp[i]]+prop[propid[i]]+delta *height[i]+b   *siteExp[i]+d*elev[i]+f*roadP[i]+g*propsize[i]+e*fc2012[i]+a*ca[i]);
  }
  
  for (k in 1:obs){
    log_lik[field[k]]=bernoulli_lpmf(y[field[k],1]|psi[field[k]]);
    
    log_lik[field[k]]=bernoulli_lpmf(y[field[k],2]|y[field[k],1]*p+(1-y[field[k],1])*q);
    
    //likelihood
  }
  //survival and detenction likelihood for the visuall class data
  for (s in 1:(n_trees-obs)){
    log_lik[remote[s]]= log_mix(psi[remote[s]],
                                bernoulli_lpmf(y[remote[s],2]|p),                                               bernoulli_lpmf(y[remote[s],2]|q));
  } 
  
}

model{
  
  p~beta(5,1.5);
  q~beta(1.5, 5);
  mug~beta(5,1.5);
  delta~normal(0,10);
  a~normal(0,10);
  b~normal(0,10);
  d~normal(0,10);
  e~normal(0,10);
  f~normal(0,10);
  g~normal(0,10);
  spel~normal(0,1);
  mu_spe~normal(0,1);
  gamma_spe~normal(0,1);
  propl~normal(0,1);
  mu_prop~normal(0,1);
  gamma_prop~normal(0,1);
  
  //survival likelihood
  
  target += sum(log_lik);
}