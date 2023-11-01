data {
int<lower=0> n_surveys;//Number of surveys
  int<lower=0> n_trees;//Number of trees
  int<lower=-1,upper=1> y[n_trees,n_surveys]; //Survival of each trees on the two surveys
  vector[n_trees] height; //predictor
}


parameters{
  
    real beta;
  real alpha;
  real<lower=0.0001,upper=0.9999> mua;
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
  //mua~beta(2,0.5);
  alpha~normal(0,1);
  beta~normal(0,1);

//survival likelihood for the remotely sensed data
  y[,2]~bernoulli_lpmf(psi);
}


