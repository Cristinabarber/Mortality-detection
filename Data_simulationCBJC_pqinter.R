library(rstan)


set.seed(2)
#total number of trees surveyed with remote sensing
n_trees <- 6229
## 329 observed trees
obs_trees <- 329
#Number of survey methods (1 = field, 2 = remote sensing)
n_surveys <- 2

# Simulate variable for the probability of detecting as alive a dead tree
#probability of false positives (e.g. detecting a tree as survived when it died)
q <- 0.27



#simulate variable for the probability of detecting as alive an alive tree
#(true positive)(p11) classified using visual classification
p <- 0.98


#Simulate variable for the probability of being alive (height m)
height <- runif( n_trees,min=50,max=3000 )
hist(height)
#standardize height
st_height <- ( height - mean(height) ) / ( 2 * sd(height ) )

#Intercept of psi0
alpha<-1
#effect of height on probability of survival
beta <-2 #
#effect of tree cover on probability of survival
#kappa<--1.2
plogis(1)
#probability of survival depends on tree height
#psi<-exp(alpha+beta*st_height+kappa*st_ndvi)/(1+exp(alpha+beta*st_height+kappa*st_ndvi))
psi <- exp( alpha + beta * st_height ) / ( 1 + exp( alpha + beta * st_height ) )

hist(psi)
plot(psi~height)
# ranging from 0.30 to 0.60 - I thought we expect psi to be heigher?

#make all true survival as 0 to start:
z <- 0
##simulate the true tree mortality
for (i in 1:n_trees){#loop over all remote sensed trees
z[i] <- rbinom( 1, 1, psi[i] ) }
#plot true survival
hist(z)
#most died in this scenario

#create observed data matrix with rows=all trees, cols=survey methods 
y <- matrix( NA, n_trees, n_surveys )
#y[,1] <- c( z[1:obs_trees], rep( -1, n_trees - obs_trees ) )
#-1 instead of NA since Stan doesn't take NAs
y[,1]<-c(z[1:obs_trees],rep(-1,n_trees-obs_trees))

##simluate mortality detection

#First we simulate the trees for which we have field data
for(tree in 1:n_trees){
  if (z[tree]==1){
    #here I want it to be a 1 (be right) with a probability of p:0.937
    y[tree,2]<-rbinom(1,1,p)}
  #here I  want it to be 1 (be wrong) with a probability of q:0.075
  else{y[tree,2]<-rbinom(1,1,q)}
}



# combine relevant data into list to import into another script:
simulation_survival_int <- list(n_trees=n_trees,n_surveys = n_surveys,
                           y = y,height=st_height, obs=obs_trees)


save(simulation_survival_int,file="simulation_survival_int.Rdata")

                 