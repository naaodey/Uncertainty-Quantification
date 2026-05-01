# Prevalence Estimation Using Bayesian and Frequentist Approaches

## Overview

This project investigates prevalence estimation using a **cutoff-free framework** based on a two-component Gamma mixture model. It compares **Bayesian inference under objective priors** with **frequentist estimation via the EM algorithm**, focusing on uncertainty quantification and model performance under varying conditions.

The work builds on existing literature in diagnostic testing and epidemiological inference, emphasizing the limitations of dichotomizing continuous measurements.

---

## Motivation

Traditional approaches to prevalence estimation rely on fixed thresholds to classify individuals as positive or negative. This leads to:

- Information loss  
- Bias in prevalence estimates  
- Ignoring uncertainty in classification  

This project adopts a **mixture modeling approach** that directly models the full distribution of diagnostic measurements, avoiding arbitrary cutoffs.

---

## Methodology

### Model Framework

- Two-component Gamma mixture model:
  - One component represents infected individuals  
  - One component represents uninfected individuals  
- Prevalence parameter \( p \) governs the mixture proportion  

The likelihood is constructed using:

\[
f(y_i) = p f_1(y_i) + (1 - p) f_0(y_i)
\]

where:
- \( f_1 \), \( f_0 \) are Gamma densities for infected and uninfected populations  

---

### Inference Approaches

#### 1. Bayesian Inference
- Priors considered:
  - Uniform prior  
  - Jeffreys’ prior  
  - Maximum Data Information Prior (MDIP)  
- Posterior inference performed via **MCMC (Stan)**  

#### 2. Frequentist Approach
- EM algorithm used for maximum likelihood estimation  

---

## Simulation Study

The study evaluates performance across:

- Sample sizes: \( n = 50, 500, 1000 \) :contentReference[oaicite:0]{index=0}  
- Different prior specifications  
- Varying degrees of separation between mixture components  

Performance metrics include:
- Bias  
- Standard deviation  
- Credible/confidence intervals  

---

## Key Results

- Bayesian and EM estimates become **nearly identical as sample size increases** :contentReference[oaicite:1]{index=1}  
- Prior choice has **minimal impact on posterior location**, especially in large samples :contentReference[oaicite:2]{index=2}  
- Bayesian methods provide **better stability in small samples**  
- Increased separation between mixture components improves estimation accuracy :contentReference[oaicite:3]{index=3}  

---

## Key Insights

- Avoiding dichotomization preserves information and improves inference  
- Mixture models provide a principled way to handle classification uncertainty  
- Bayesian inference is particularly advantageous in:
  - Low-prevalence settings  
  - Small sample sizes  
  - Weak identifiability scenarios  
- Asymptotically, Bayesian and frequentist methods converge  

---

## Limitations

- Assumes correct Gamma model specification  
- Sensitive to overlap between mixture components  
- Does not include:
  - Covariates  
  - Hierarchical structure  
  - Dependence between observations  
- Mixture models may suffer from identifiability issues  

---

## Files

- `code/` – Stan models and simulation scripts  
- `report.pdf` – Full paper   

---

## Tools

- R  
- Stan (via rstan)  
- Bayesian inference techniques  

---

## References

- Gelman, A., & Carpenter, B. (2020)  
- Bouman et al. (2020)  
- Dempster et al. (1977)  
- Bernardo (1979), Jeffreys (1946), Kass & Wasserman (1996)  

---

## Future Work

- Extend to hierarchical mixture models  
- Incorporate covariates  
- Explore model robustness under misspecification  
- Apply to real-world epidemiological data  

---

## Author

Naa Odey Solomon  
MPhil in Mathematical Statistics  
