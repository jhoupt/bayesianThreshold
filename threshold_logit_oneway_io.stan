data{  // Everything that must be input by the user
   int<lower=4> nTotal;                         // Total number of trials to analyze
   int<lower=2> nLevels;                        // Number of levels of Factor 
   real intensity[nTotal];                      // Intensity of the stimulus on each trial
   int<lower=0> level[nTotal];                 // Level of Factor on each trial
   int<lower=0,upper=1> correct[nTotal];        // Whether the response was correct (1) on each trial
   real<lower=0,upper=1> chancePerformance;     // Chance performance for experiment (e.g., 1/n for n-AFC)
}
parameters{  // Everything that has a prior distribution
   real mum;
   real muw;
   vector[nLevels-1] fA;
}
transformed parameters {
   real m[nLevels];
   real<lower=0> w[nLevels];
   vector[nLevels] factorA;

   if (nLevels > 2) {
      for (l in 1:(nLevels-1)) { 
         factorA[l] <- fA[l]/10;
      }
      factorA[nLevels] <- -1*sum(fA)/10;
   }
   else {
      factorA[1] <- fA[1]/10;
      factorA[2] <- -1*fA[1]/10;
   }


   for (l in 1:nLevels) {
      m[l] <- mum + factorA[l] ;
      w[l] <- muw ;
   }
   
}
model {
   real threshold[nTotal];
   real width[nTotal];
   vector[nTotal] psi;

   for (tr in 1:nTotal) {
      threshold[tr] <- m[level[tr]];
      width[tr] <- w[level[tr]];
      psi[tr] <- ( (1-chancePerformance)*inv_logit(4.4/width[tr] * (intensity[tr]-threshold[tr])) + chancePerformance);
   }

   fA ~ normal(0,1);
   muw ~ gamma(2,.5);
   correct ~ bernoulli(psi);
}
