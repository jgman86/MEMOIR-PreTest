// Implemented in Stan by Jan Goetttmann 
// M3 - Model for extended encoding in complex span tasks after Oberauer & Lewandowsky, 2018
// Questions regarding the model code to Jan.Goettmann@psychologie.uni-heidelberg.de

// Best Fitting Model for Complex Span


functions{
  // flatten_lower_tri: function that returns the lower-tri of a matrix, flattened to a vector
  
  vector flatten_lower_tri(matrix mat) {
    int n_cols = cols(mat) ;
    int n_uniq = (n_cols * (n_cols - 1)) %/% 2;
    vector[n_uniq] out ;
    int i = 1;
    for(c in 1:(n_cols-1)){
      for(r in (c+1):n_cols){
        out[i] = mat[r,c];
        i += 1;
      }
    }
    return(out) ;
  }
  
}

data {
  int <lower=0> N;  // number rownumber
  int <lower=0> K;  // categories 
  int <lower=0> J;  // Dims of Cov Matrix
  int R[K];         // number of responses per category
  int count[N,K];   // observed data
  real scale_b; // set scaling for background noise
  int retrievals;
  matrix[N,2] d_weight; // Weights for Response Category inbalances
}

parameters {
  
  // Defining vector for hyper and subject parameters 
  
  cholesky_factor_corr[J] L_Omega;
  vector<lower=0>[J] sigma;
  vector[J] hyper_pars;
  matrix[J,N] theta;
  
}


transformed parameters {
  // non-centered multivariate
  matrix[J,N] subj_pars =  (
    diag_pre_multiply( sigma, L_Omega )
    * theta
    + rep_matrix(hyper_pars,N)
    ) ;
    
    // Transform f Parameter
    real mu_f = inv_logit(hyper_pars[3]);
    row_vector[N] f = inv_logit(subj_pars[3,]);
    
   
    // activations
    real acts_IIP[N];
    real acts_IOP[N];
    real acts_DIP[N];
    real acts_DIOP[N];
    real acts_NPL[N];
    
    
    // probabilities
    vector[K] probs[N];
    real SummedActs[N];
    

    // loop over subjects and conditions to compute activations and probabilites
    
    
    for (i in 1:N){ // for each subject
      
      acts_IIP[i] = scale_b +  subj_pars[1,i] + subj_pars[2,i]; // Item in Position                      
      acts_IOP[i] = scale_b + subj_pars[2,i];        // Item in Other Position
      acts_DIP[i] = scale_b + f[i]*(subj_pars[1,i] + subj_pars[2,i]); // Distractor in Position
      acts_DIOP[i] = scale_b + f[i]*subj_pars[2,i]; // Distractor in other Position
      acts_NPL[i] = scale_b; // non presented Lure
      
      SummedActs[i] = R[1] * acts_IIP[i] + R[2] * acts_IOP[i] + d_weight[i,1] * acts_DIP[i] + d_weight[i,2] * acts_DIOP[i]+ R[5] * acts_NPL[i];
      
      probs[i,1] = (R[1] * acts_IIP[i]) ./ (SummedActs[i]);  
      probs[i,2] = (R[2] * acts_IOP[i]) ./ (SummedActs[i]);
      probs[i,3] = (d_weight[i,1] * acts_DIP[i]) ./ (SummedActs[i]);
      probs[i,4] = (d_weight[i,2] * acts_DIOP[i]) ./ (SummedActs[i]);
      probs[i,5] = (R[5] * acts_NPL[i]) ./ (SummedActs[i]);
    }
}


model {
  
  // priors for hyper parameters
  
  hyper_pars[1] ~ normal(20,10); // c
  hyper_pars[2] ~ normal(2,10); // a
  hyper_pars[3] ~ normal(0,10); // f

  // priors for covariance matrix
  L_Omega ~ lkj_corr_cholesky(2);
  sigma ~ gamma(1,0.01);
  
  
  // Loop over subjects
  
  for (i in 1:N) 
  {
    
    theta[,i] ~ normal(0,1);
    
  }
  
  
  

    for (i in 1:N) {
      // draw data from probabilities determined by MMM parms
      count[i,] ~ multinomial(probs[i,]);  
    }
  
}

generated quantities{
  
  vector[(J*(J-1))%/%2] cor_mat_lower_tri;
  int count_rep[N,K];
  
  
  cor_mat_lower_tri = flatten_lower_tri(multiply_lower_tri_self_transpose(L_Omega));
  
  
  
  
  for (i in 1:N)
 
    {
      
      count_rep[i,] = multinomial_rng(probs[i,], retrievals);
      
    }
  
}

