model {
  # Y ~  
  # prior population
  
  mu.beta0 ~ dnorm(0,0.0001)      # Intercepto
  sigma.beta0 ~ dunif(0.0001,100)
  
  mu.beta1 ~ dnorm(0,0.0001)       # EU
  sigma.beta1 ~ dunif(0.0001,100)  
  

  # learning rate feedback + 
  alpha_a.kappa    ~ dunif(0.0001, 100)
  alpha_a.mu ~ dbeta(1,1)
  alpha_a.a   <- alpha_a.mu *alpha_a.kappa 
  alpha_a.b   <- (1-alpha_a.mu) *alpha_a.kappa 
  # learning rate feedback - 
  alpha_r.kappa    ~ dunif(0.0001, 100)
  alpha_r.mu ~ dbeta(1,1)
  alpha_r.a   <- alpha_r.mu *alpha_r.kappa 
  alpha_r.b   <- (1-alpha_r.mu) *alpha_r.kappa 
  
  # alpha of prospect theory v(x) = x^alpha
  mu.rew ~ dunif(0.01,10) 
  sigma.rew ~ dunif(0.001,100) 
  
  # gamma of prospect theory pi(p)
  #mu.gamma ~ dunif(0.01,10) 
  #sigma.gamma ~ dunif(0.001,100) 
  
  
  
  
  for (SU in 1:nSUB){
    
    
    beta0[SU] ~ dnorm(mu.beta0, 1/sigma.beta0^2)
    beta1[SU] ~ dnorm(mu.beta1, 1/sigma.beta1^2)
    alpha_a[SU] ~ dbeta(alpha_a.a,alpha_a.b)
    alpha_r[SU] ~ dbeta(alpha_r.a,alpha_r.b)
    rew[SU] ~ dnorm(mu.rew,1/sigma.rew^2)T(0.01,10)
    #gamma[SU] ~ dnorm(mu.gamma,1/sigma.gamma^2)T(0.01,10)
  }
  
  p_sub[1]=0.5
  
  for (i in 2:nT) {
    
    feedback_left[i-1] = (feedback[i-1]*y[i-1]) + ((1-feedback[i-1])*(1-y[i-1]))
    
    delta_a[i] <- alpha_a[s[i]] * ((feedback_left[i-1]) - p_sub[i-1])
    delta_r[i] <- alpha_r[s[i]] * ((feedback_left[i-1]) - p_sub[i-1])
    
    p_sub[i] <- (0.5 * firstT[i]) +  ( p_sub[i-1] + delta_a[i]   )*(1-firstT[i])*(feedback[i-1]) + ( p_sub[i-1] + delta_r[i]   )*(1-firstT[i])*(1-feedback[i-1])
    
    # pi(p_sub)
    Pi.1[i] <- p_sub[i]# ^ gamma[s[i]]   /   (( p_sub[i]^gamma[s[i]]  + (1-p_sub[i])^gamma[s[i]]   )^(1/gamma[s[i]]))
    Pi.2[i] <- (1-p_sub[i])# ^ gamma[s[i]]   /   (( (1-p_sub[i])^gamma[s[i]]  + (p_sub[i])^gamma[s[i]]   )^(1/gamma[s[i]]))
    
    #
    V.1[i] <- Bd1[i] ^ rew[s[i]]
    V.2[i] <- Bd2[i] ^ rew[s[i]]
    

     U[i] <- V.1[i]  * Pi.1[i]  - V.2[i] * Pi.2[i] 
     
     p_logic[i] <- 1 /(1 + exp(-( beta0[s[i]]  + beta1[s[i]]  * U[i] )))
    
     y[i] ~ dbern(p_logic[i]) 
    
  }
  
}