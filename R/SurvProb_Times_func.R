# Given estimated survival probabilities, write a function that obtains survival
# probabilities at specific times

# -- INPUT:
# Times_Specify: Times that we want to estimate the survival probability (Vector)
# Est_Times: Times that the survival probability is estimated at (Vector)
# Est_Surv: Estimated survival probabilities at time "Est_Times" (Vector)

# -- OUTPUT:
# A dataframe containing
# (i) Times_Specify: Event times
# (ii) Surv_Probs: Survival probability evaluated at "Times_Specify"

SurvProb_Times_func = function(Times_Specify,
                               Est_Times,
                               Est_Probs){

  val = sapply(Times_Specify,
               function(s){

                 val.s = 1 # Initialize

                 ind = which(Est_Times <= s)

                 if (length(ind) > 0){
                   val.s = Est_Probs[max(ind)]
                 }
                 return(val.s)
               })
  return(val)
}
