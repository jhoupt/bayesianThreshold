data{  // Everything that must be input by the user
   int<lower=0> nTotal;                         // Total number of trials to analyze
   int<lower=0> nSubj;                          // Number of unique observers
   int<lower=0> nFactorA;                       // Number of levels of Factor A
   int<lower=0> nFactorB;                       // Number of levels of Factor B
   real<lower=0> intensity[nTotal];             // Intensity of the stimulus on each trial
   int<lower=0> subject[nTotal];                // Observer on each trial
   int<lower=0> factorA[nTotal];                // Level of Factor A on each trial
   int<lower=0> factorB[nTotal];                // Level of Factor B on each trial
   int<lower=0,upper=1> correct[nTotal];        // Whether the response was correct (1) on each trial
   real<lower=0,upper=1> chancePerformance;     // Chance performance for experiment (e.g., 1/n for n-AFC)
}
parameters{  // Everything that has a prior distribution
   vector<lower=0,upper=1>[nSubj] lapse;        // Observer's lapse rate
   real<lower=0,upper=1> betaMean;              // Mean of beta prior for individual lapse rate
   real<lower=0> betaEta;                       // Precision parameter of beta prior for individual lapse rate
   real mum;
   real muw;
   real subj_alpha[nSubj];
   real factorA_alpha[nFactorA];
   real factorB_alpha[nFactorB];
   real subj_beta[nSubj];
   real factorA_beta[nFactorA];
   real factorB_beta[nFactorB];
}
transformed parameters {
   vector<lower=0,upper=1>[nTotal] psi;
   real<lower=0> lapseAlpha;
   real<lower=0> lapseBeta;
   real<lower=0> m[nSubj,nFactorA,nFactorB];
   real<lower=0> w[nSubj,nFactorA,nFactorB];

   real<lower=0> threshold[nTotal];
   real<lower=0> width[nTotal];

   for (sj in 1:nSubj) {
      for (cd in 1:nFactorA) {
         for (dy in 1:nFactorB) {
            m[sj,cd,dy] <- mum + factorA_alpha[cd] + factorB_alpha[dy] + sA[sj];
            w[sj,cd,dy] <- muw + factorA_beta[cd] + factorB_beta[dy] + sB[sj];
         }
      }
   }

   for (tr in 1:nTotal) {
      threshold[tr] <- m[subject[tr],factorA[tr],factorB[tr]];
      width[tr] <- w[subject[tr],factorA[tr],factorB[tr]];
      psi[tr] <- (1-lapse[ subject[tr] ])*( (1-chancePerformance)*inv_logit(4.4/width[tr] * (intensity[tr]-threshold[tr])) + chancePerformance) + chancePerformance*lapse[ subject[tr] ];
   }

   betaMean  <- 0.01;
   betaEta <- 100;
   lapseAlpha <- betaMean * betaEta;
   lapseBeta <- (1-betaMean) * betaEta ;
}
model {
   betaMean ~ beta(1,60);
   betaEta ~ gamma(1,0.01);
   factorA_alpha ~ normal(0,0.001);
   factorA_beta ~ normal(0,0.001);
   factorB_alpha ~ normal(0,0.001);
   factorB_beta ~ normal(0,0.001);
   sA ~ normal(0,0.001);
   sB ~ normal(0,0.001);
   lapse ~ beta(lapseAlpha,lapseBeta);
   correct ~ bernoulli(psi);
}
