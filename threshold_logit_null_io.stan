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
}
model {
   real threshold[nTotal];
   real width[nTotal];
   vector[nTotal] psi;

   for (tr in 1:nTotal) {
      threshold[tr] <- mum;
      width[tr] <- muw;
      psi[tr] <- ( (1-chancePerformance)*inv_logit(4.4/width[tr] * (intensity[tr]-threshold[tr])) + chancePerformance);
   }

   muw ~ gamma(2,.5);
   correct ~ bernoulli(psi);
}
