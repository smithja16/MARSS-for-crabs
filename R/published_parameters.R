## ============================================================================
##  Parameter estimates of the final model (Table S2 of the paper)
## ============================================================================
##  Used by simulate_data.R to generate the synthetic dataset, and by
##  run_example.R to check how well the parameters are recovered.
##
##  Two latent states (1 = male, 2 = female). Six observation series, in order:
##    1 commercial CPUE, male          (kg per trap lift)
##    2 commercial CPUE, female
##    3 survey recruit CPUE, male      (number per trap-night)
##    4 survey recruit CPUE, female
##    5 survey prerecruit CPUE, male, lagged 2 months
##    6 survey prerecruit CPUE, female, lagged 4 months
##  Four covariates, in order: bottom temperature (2-month lag, standardised),
##  conductivity (standardised), cos(2*pi*month/12), sin(2*pi*month/12).
## ============================================================================

published <- list(
  B  = matrix(c(0.532993, 0.221751, -0.136171, 0.638332), 2, 2),
  U  = matrix(c(0.0961587, -0.0159684), 2, 1),
  C  = matrix(c(-0.0318213, -0.00187693, -0.00847375, 0.00688591,
                 0.0275257, -0.00272841,  0.0337358, -0.0228486), 2, 4),
  Q  = matrix(c(0.00290904, -0.000785211, -0.000785211, 0.000484427), 2, 2),
  Z  = matrix(c(1, 0, 2.27562, 0, 13.1841, 0,
                0, 1, 0, 11.0679, 0, 38.7519), 6, 2),
  R  = matrix(c(0.00531264, 0, 0, 0, 0, 0,
                0, 0.00124841, 0, 0, 0, 0,
                0, 0, 0.0630355, 0, 0.0348172, 0,
                0, 0, 0, 0.0866435, 0, -0.332502,
                0, 0, 0.0348172, 0, 3.48068, 0,
                0, 0, 0, -0.332502, 0, 1.97911), 6, 6),
  x0 = matrix(c(0.293087, 0.112652), 2, 1)
)
