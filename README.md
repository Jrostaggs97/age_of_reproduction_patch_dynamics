# Multi-Species Competition with Maturation Delay and Patch Dynamics
*A compact guide to the model, solver, and repo structure.*

* This code adapted from the code found at: https://github.com/utrigos/Disturbance_Generated_Competitive_Coexistence_Code 
---

## What this code does

This repository simulates an **age-structured system of PDEs** describing the dynamics of \(k\) competing species that:

* delay reproduction until a species-specific age $$\tau_i$$,
* experience constant adult mortality $$\mu_i$$,
* reproduce at rate $$b_i$$ on patches that are connected via dispersal
* patches age and reset via disturbance, and
* species experience density dependence on reproduction proportional to coefficients $$\alpha_i$$.

Core workflow  

1. **Build an age grid** and the patch-age density $$\rho(a)$$.  
2. **Integrate in time** with MATLAB’s `dde23` (Runge–Kutta for delay equations).  
3. **Handle aging** in with a flux-limiter scheme (Koren or minmod).
4. **Handle demographic processes** with non-local reproduction and mortality.
5. **Post-process** total abundances, invasion scores, coexistence maps, plots.  

---

## Patch-age distribution

Patch ages evolve via a **McKendrick–von Foerster equation**

$$
\partial_t \rho(a,t) + \partial_a \rho(a,t) = -\gamma \rho(a,t), 
$$

with

$$
\rho(a,0)=\rho_0(a), 
\qquad
\rho(0,t)=\int_0^\infty \gamma\rho(s,t)\,ds .
$$

We typically assume γ is constant and work at steady state, giving the exponential density  

$$
\rho(a) = \gamma e^{-\gamma a}, \qquad a\ge0.
$$

This is exactly the weighting used in reproduction integrals and total-abundance calculations.

---

## Population Dynamics for generic species $$i$$

$$
\partial_t n_i(a,t) =
-\partial_a n_i    
+
H(a-\tau_i)\bigg[
b_i \int_{\tau_i}^{\infty} \rho(a)\,n_i(a,t-\tau_i)(1-\alpha_i\sum_{j=1}^{k} n_j(a,t-\tau_i)) \mathrm{d}s - \mu_i n_i  \bigg]
$$


Boundary and initial condition 

$$
n_i(0,t)=0,
\qquad
n_i(a,t\le0)=n_{i,0}(a).
$$

### Term-by-term cheat-sheet

| symbol / term | biological role |
|---------------|-----------------|
| $$-\partial_a n_i$$ | aging / advective transport along age axis |
| $$-\mu_i n_i$$ | constant adult mortality |
| $$b_i\int_{\tau_i} \rho(a) n_i(\cdot)\mathrm{d}s$$ | birth flux modulated by patch age |
| $$\int [...]n(1-\alpha_i\sum_j n_j) \mathrm{d}s$$ | density-dependent competition |
| $$H(a-\tau_i)$$ | enforces maturation delay $$\tau_i$$ |

---

## How the solver fits together

```text
┌─────────────────────────┐
│ 1. User parameters      │  (b, μ, τ, α, γ, grids …)
└────────────┬────────────┘
             │ build a, ρ(a)
┌────────────▼────────────┐
│ 2. Time integration     │  dde23 (RK w/ delays)
└────────────┬────────────┘
             │ RHS evaluation each step
             │ ┌──────────────────────────┐
             │ │ loop over species i      │
             │ │   • fetch delayed state  │
             │ │   • reproduction term    │
             │ │   • flux-limited aging   │
             │ │   • assemble ∂ₜn_i       │
             │ └──────────────────────────┘
┌────────────▼────────────┐
│ 3. Store n(a,t)         │
└────────────┬────────────┘
             │
┌────────────▼────────────┐
│ 4. Post-processing      │  totals, coexistence maps, plots
└─────────────────────────┘
```

### File Guide

| file                                        | purpose                                                |
| ------------------------------------------- | -------------------------------------------------------|
| `dd_reproduction_competition_dynamics.m`    | top-level wrapper: builds grids, calls `dde23`         |
| `dd_reproduction_rhs_fun.m`                 | right-hand side supplied to `dde23`                    |
| `koren_flux.m`, `minmod_flux.m`             | high-resolution flux reconstructions for aging term    |
| `generate_random_initial_conditions.m`      | helper – Gaussian initial profiles                     |
| `b_tau_coex_check.m`, `mu_tau_coex_check.m` | analytic invasion-score calculators                    |
| `full_dynamics_demo.m`                      | runnable demo producing $N_i(t)$ & age-profile plots   |
| `coex_region_demo.m`                        | runnable demo producing $\mu- \tau$ coexistence region |
