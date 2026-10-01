# MARSS-for-crabs

A worked example of a multivariate autoregressive state-space (MARSS) model that
integrates commercial and fishery-independent survey catch rates, and of two ways to
measure what the survey adds. It accompanies Smith, Johnson & Taylor (in prep.)
*Estimating the value of a fishery-independent survey to catch rate assessment using a multivariate time series analysis*, which applies the method to blue swimmer crab (*Portunus armatus*) in
Wallis Lake, New South Wales.

**This repository runs on synthetic data.** The commercial catch data used in the paper
are confidential. `simulate_data.R` generates a dataset with the same structure - 67
monthly observations of six series, the same months without survey data, and similar
covariates - by simulating from the parameter estimates published in Table S2 of the
paper. The workflow therefore runs end to end, but **it does not reproduce the paper's
numbers.** The real data can be requested from the corresponding author:
**james.a.smith@dpird.nsw.gov.au**.

---

## Quick start

```r
source("simulate_data.R")   # writes example_data/
source("run_example.R")     # about 10-15 minutes; prints results, writes output/
```

Requires R ≥ 4.3 with the **MARSS** (3.11.10) and **sandwich** packages.

---

## The model

Two latent states - male and female abundance available to the fishery, on the
commercial catch-rate scale - observed monthly by six series:

| Series | Units | Informs |
|---|---|---|
| Commercial CPUE, male / female | kg per trap lift | male / female state; loading fixed at 1 |
| Survey recruits (≥ 65 mm), male / female | number per trap-night | male / female state |
| Survey prerecruits (< 65 mm), male | number per trap-night | male state, **2 months later** |
| Survey prerecruits (< 65 mm), female | number per trap-night | female state, **4 months later** |

Prerecruits are a leading indicator: crabs below legal size now are the recruits of the
coming months. The lags are applied by shifting the data, because the observation
equation can only link observations and states in the same month (see
`R/model_spec.R`). Bottom temperature, conductivity and a seasonal cycle affect the
states.

## What `run_example.R` does

1. **Fits the full model**, from two starting points, and compares the estimates with
   the values the data were simulated from.
2. **Checks the estimated states against the truth.** With synthetic data the true
   abundance is known, so you can see directly how much closer the estimates get when
   the survey is included. The real analysis cannot do this.
3. **Measures survey value by forecast accuracy.** At each forecast origin both models
   are refitted and used to forecast commercial catch rates one to three months ahead.
   The commercial-only model cannot estimate its own observation error, so it is given
   the full model's estimate - a deliberate help to the baseline. Only the commercial
   series are scored, so every model is judged on the same values. A seasonal naïve
   forecast is included as a reference.
4. **Measures survey value by estimation precision.** Survey observations are removed
   from randomly chosen months and the Kalman smoother is re-run with every parameter
   held fixed, so the only thing that changes is the information available.

## Files

| File | Purpose |
|---|---|
| `simulate_data.R` | Generates the synthetic data and the true states |
| `run_example.R` | Runs the four steps above; settings at the top |
| `R/model_spec.R` | Data preparation and the model, annotated matrix by matrix |
| `R/fitting.R` | Two-stage fitting from multiple starts, and forecasting |
| `R/survey_value_cv.R` | Cross-validation, baselines and the Diebold–Mariano test |
| `R/survey_value_se.R` | The survey-removal precision test |
| `R/published_parameters.R` | Parameter estimates used for the simulation |

The paper's model selection - the choice of temperature series and lag, the error
structures, and temperature as a catchability effect - is reported in its
supplementary material and is not repeated here.

---

## Notes for adapting this to your own data

- **Fit from more than one start.** MARSS can report convergence at a local optimum.
  Fitting the same model with standardised and raw covariates is an exact
  reparameterisation that gives the optimiser two different routes.
- **Compare models with different series on forecast accuracy, not AIC.** Likelihoods
  from different observation matrices are on different scales. Even among comparable
  models, whole-model AIC can favour fit to series that are not the target of inference.
- **Give forecasts covariates on the scale the model was fitted on.** `fit_marss()`
  records which start won for this reason.
- **Gaussian errors permit negative values.** The simulated catch rates are floored at
  zero; fitted values and intervals may still dip below it.

## Use of generative AI

The code was developed with assistance from AI tools (ChatGPT, Claude), as declared in
the manuscript. All analyses were directed and checked by the authors.
