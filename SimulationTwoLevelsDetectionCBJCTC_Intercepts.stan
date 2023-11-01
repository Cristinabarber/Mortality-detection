data {
  int<lower=0> n_surveys;//Number of surveys
  int<lower=0> n_trees;//Number of trees
  int<lower=0> obs;//number of  observations in the groundtruthed data
  int<lower=-1,upper=1> y[n_trees,n_surveys]; //Survival of each trees on the two surveys
  vector[n_trees] height; //predictor
}


parameters{
  real beta;
  real alpha;
  real<lower=0.0001,upper=0.9999> mua;
  real<lower=0,upper=1>  p;   //probability of detecting as alive an alive tree(truly alive)
  real<lower=0,upper=p>  q;   //probability of detencting as alive a dead tree
  
  
}

transformed parameters{
  real<lower=0,upper=1>  psi[n_trees]; //probability survival using groundthruthed data
  //real alpha;


   
  //alpha=log(mua/(1-mua));

  for (n in 1:n_trees){
    psi[n]=inv_logit(alpha+ beta*height[n]);

    }
    
    
}

model{
  real log_psi[n_trees];
  real log1m_psi[n_trees];
  log_psi=log(psi);     //log transformation for the probability of surviving
  log1m_psi=log1m(psi); //Log transformation for the probability of not surviving
  
  p~beta(6,1);
  q~beta(1,6);
  //mua~beta(2,0.5);
  alpha~normal(0,1);
  beta~normal(0,1);
  


//survival likelihood for the field data
//survival likelihood for the groundtruthed data
for (k in 1:obs){
      y[k,1]~bernoulli(psi[k]);

      y[k,2]~bernoulli(y[k,1]*p+(1-y[k,1])*q);

//likelihood
}
//survival and detenction likelihood for the visuall class data
for (i in n_trees-obs:n_trees){
        target += log_mix(psi[i],
                          bernoulli_lpmf(y[i,2]|p),
                          bernoulli_lpmf(y[i,2]|q));

}               


}

