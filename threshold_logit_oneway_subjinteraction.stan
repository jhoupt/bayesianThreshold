data{  // Everything that must be input by the user
   int<lower=4> nTotal;                         // Total number of trials to analyze
   int<lower=2> nSubj;                          // Number of unique observers
   int<lower=2> nLevels;                        // Number of levels of Factor 
   real intensity[nTotal];                      // Intensity of the stimulus on each trial
   int<lower=0> subject[nTotal];                // Observer on each trial
   int<lower=0> level[nTotal];                 // Level of Factor on each trial
   int<lower=0,upper=1> correct[nTotal];        // Whether the response was correct (1) on each trial
   real<lower=0,upper=1> chancePerformance;     // Chance performance for experiment (e.g., 1/n for n-AFC)
}
transformed data{
   int levelAS[nTotal];
   int nLevelAS;
   // Contrast for interaction between factor level and subject 
   for ( tr in 1:nTotal) {
      levelAB[tr] <-  (level[tr] -1) * nSubj + subject[tr]
   }
   nLevelAS <- (nLevels - 1) * (nSubj -1)
}
parameters{  // Everything that has a prior distribution
   vector<lower=0,upper=1>[nSubj] lapse;        // Observer's lapse rate
   real mum;
   real muw;
   vector[nLevels-1] fA;
   vector[nSubj-1] sA;
   vector[nSubj-1] sB;
   vector[nLevelAS] fsA;
}
transformed parameters {
   real<lower=0,upper=1> betaMean;              // Mean of beta prior for individual lapse rate
   real<lower=0> betaEta;                       // Precision parameter of beta prior for individual lapse rate

   vector<lower=0,upper=1>[nTotal] psi;

   real<lower=0> lapseAlpha;
   real<lower=0> lapseBeta;

   real m[nSubj,nLevels];
   real<lower=0> w[nSubj,nLevels];
   real threshold[nTotal];
   real<lower=0> width[nTotal];

   vector[nLevels] factorA;
   vector[nSubj] subjectA;
   vector[nSubj] subjectB;

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

   if (nSubj > 2) {
      for (l in 1:(nSubj-1)) { 
         subjectA[l] <- sA[l]/10;
         subjectB[l] <- sB[l]/10;
      }
      subjectA[nSubj] <- -1*sum(sA)/10;
      subjectB[nSubj] <- -1*sum(sB)/10;
   }
   else {
      subjectA[1] <- sA[1]/10;
      subjectA[2] <- -1*sA[1]/10;
      subjectB[1] <- sB[1]/10;
      subjectB[2] <- -1*sB[1]/10;
   }


   for (sj in 1:nSubj) {
      for (l in 1:nLevels) {
         m[sj,l] <- mum + factorA[l]  + subjectA[sj];
         w[sj,l] <- muw + subjectB[sj];
      }
   }

   for (tr in 1:nTotal) {
      threshold[tr] <- m[subject[tr],level[tr]];
      width[tr] <- w[subject[tr],level[tr]];
      psi[tr] <- (1-lapse[ subject[tr] ])*( (1-chancePerformance)*inv_logit(4.4/width[tr] * (intensity[tr]-threshold[tr])) + chancePerformance) + chancePerformance*lapse[ subject[tr] ];
   }

   betaMean  <- 0.01;
   betaEta <- 100;
   lapseAlpha <- betaMean * betaEta;
   lapseBeta <- (1-betaMean) * betaEta ;
}
model {
   //betaMean ~ beta(1,60);
   //betaEta ~ gamma(1,0.01);
   fA ~ normal(0,1);
   sA ~ normal(0,1);
   sB ~ normal(0,1);

   muw ~ gamma(2,.5);
   lapse ~ beta(lapseAlpha,lapseBeta);
   correct ~ bernoulli(psi);
}
