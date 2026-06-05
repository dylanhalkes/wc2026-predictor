# 2026 FIFA World Cup Predictor

A Monte Carlo tournament simulator for the 2026 FIFA World Cup, built in R. The model estimates each team's attack and defence strengths from historical match data, simulates the full 48-team tournament 50,000 times, and produces win probabilities and stage-by-stage forecasts for each competing nation.

An interactive Shiny dashboard allows users to explore predictions, simulate individual matchups, and view Dixon-Coles adjusted scoreline probability heat maps.

> **Status:** v2 predictions live (Poisson GLM + Elo ratings + Dixon-Coles correction + host nation advantage). Tournament starts 11 June 2026.

---

## Methodology

### Data
International football results from 1872 to 2026 (Kaggle - Mart Jürisoo), filtered to matches from January 2018 onwards to capture the last two World Cup qualifying cycles. Raw data comprises ~8,000 matches across 284 nations.

### Match Weighting
A two-dimensional weighting scheme is applied to each observation:

**Match type weight** - reflects competitive importance:
| Tournament | Weight |
|---|---|
| FIFA World Cup | 4.0 |
| Continental championship (Euros, Copa América, AFCON, AFC Asian Cup, Gold Cup) | 3.0 |
| World Cup qualification | 2.5 |
| Continental qualification | 2.0 |
| UEFA Nations League / CONCACAF Nations League | 1.5 |
| Friendly | 0.5 |

**Time decay weight** - recent results carry more weight via exponential decay with a two-year half-life:

$$w_t = e^{-\frac{\ln 2}{730} \cdot d}$$

where $d$ is days before the tournament start date. A match from one year ago receives approximately 71% of the weight of a match played today.

### Poisson Goals Model
Goals scored are modelled as Poisson distributed with a log-linear conditional mean:

$$\log(\lambda_{ij}) = \mu + \alpha_i + \beta_j + \gamma \cdot \text{home}_i + \delta \cdot \text{EloDiff}_{ij}$$

where:
- $\lambda_{ij}$ = expected goals for team $i$ against team $j$
- $\mu$ = intercept (baseline log goal rate)
- $\alpha_i$ = attack strength of team $i$
- $\beta_j$ = defensive weakness of team $j$
- $\gamma$ = home advantage coefficient
- $\delta$ = Elo difference coefficient
- $\text{EloDiff}_{ij}$ = (Elo rating of team $i$ − Elo rating of team $j$) / 100 at the time of the match

The model is fitted as a weighted Poisson GLM. Elo ratings are sourced from the World Football Elo Ratings dataset (Kaggle - afonsofernandescruz), covering all 48 qualified teams with annual readings from 1901 to 2026. Non-WC opponents in training data are assigned the median Elo.

**Key estimated parameters:**
- Home advantage $e^\gamma \approx 1.27\times$ goal multiplier
- Elo coefficient $\delta = -0.032$ (per 100 Elo point advantage, −3.1% expected goals after controlling for fixed effects - the fixed effects already capture most team strength, so the residual Elo signal is small but improves AIC by 243 points)

### Dixon-Coles Correction
The standard Poisson model underestimates the frequency of 0-0 and 1-1 scorelines and overestimates 0-1 and 1-0. Following Dixon & Coles (1997), we apply a correction factor $\tau$ to these four low-score cells:

$$\tau(x,y) = \begin{cases} 1 - \mu\nu\rho & x=0,\ y=0 \\ 1 + \mu\rho & x=0,\ y=1 \\ 1 + \nu\rho & x=1,\ y=0 \\ 1 - \rho & x=1,\ y=1 \\ 1 & \text{otherwise} \end{cases}$$

The parameter $\rho$ is estimated by maximising the DC log-likelihood correction holding the Poisson expected goals fixed: $\hat{\rho} = -0.051$. This increases P(0-0) and P(1-1) by approximately +0.67 percentage points each for a typical match. The correction is applied to all individual match probability calculations within the Shiny dashboard.

### Monte Carlo Simulation
The full 48-team tournament is simulated 50,000 times:

1. **Group stage** - all six round-robin matches per group simulated by drawing from $\text{Poisson}(\lambda_{ij})$. Teams ranked by points → goal difference → goals scored → random draw.
2. **Third-place selection** - all 12 third-placed teams ranked; best 8 advance.
3. **Round of 32** - bracket assembled using FIFA's official Annex C structure, with a constraint-satisfaction algorithm (minimum remaining values) for third-place assignments.
4. **Knockout rounds** - draws resolved by a Poisson-weighted penalty shootout.
5. **Host nation advantage** - USA, Canada and Mexico receive a partial home advantage of $0.5\gamma$ applied to their expected goals, reflecting crowd support and familiarity without the full home advantage (all matches are technically neutral-ground fixtures).

Standard error of probability estimates at 50,000 simulations: ±0.18 percentage points for a team with ~20% win probability.

---

## Results (v2 - Poisson GLM + Elo + Dixon-Coles + Host Advantage)

Pre-tournament predictions published **5 June 2026**, six days before kick-off.

| Rank | Team | Win % | Final % | Semi % | QF % |
|---|---|---|---|---|---|
| 1 | Spain | 18.1% | 28.0% | 41.5% | 54.1% |
| 2 | England | 14.0% | 23.7% | 36.7% | 55.6% |
| 3 | Argentina | 11.8% | 19.6% | 32.5% | 46.6% |
| 4 | Brazil | 7.9% | 15.3% | 27.6% | 41.6% |
| 5 | France | 7.9% | 15.3% | 26.2% | 45.7% |
| 6 | Portugal | 6.3% | 12.5% | 24.4% | 43.8% |
| 7 | Germany | 5.3% | 12.2% | 24.3% | 47.4% |
| 8 | Belgium | 4.3% | 9.2% | 19.3% | 44.1% |
| 9 | Morocco | 3.9% | 8.6% | 17.4% | 29.7% |
| 10 | Netherlands | 3.5% | 7.6% | 15.8% | 26.4% |

**Notable findings:**
- Spain's model-implied win probability (18.1%) reflects both the highest attack parameter among WC teams as well as elite defensive stability, combined with a favourable Group H path
- England sit second despite a lower FIFA ranking - the model rewards their high-scoring recent qualifying campaign and strong fixed-effect parameters
- Germany and Belgium show notably high QF probabilities (47.4% and 44.1%) relative to their win probabilities (5.3% and 4.3%), reflecting bracket paths that route them into stronger opposition in later rounds
- The negative Elo coefficient (δ = −0.032) is an interesting finding: after controlling for team fixed effects, Elo provides marginal additional signal. The fixed effects already capture most long-run strength; the residual Elo effect may reflect conservative game management by strong teams against weak opponents in qualifying

**v1 vs v2 comparison:**
The top five teams are stable across both model versions (within 1–2 percentage points), indicating the baseline Poisson model was already well-calibrated. The main movers in v2 are mid-tier teams with strong current Elo relative to their historical fixed effects, such as Belgium (+1.4pp win), Norway, and Ecuador - all reflecting good recent form captured by Elo but not fully by historical results alone.

---

## Planned Improvements

- **Post-tournament calibration analysis** — comparing predicted vs actual outcomes across all 104 matches to assess model accuracy
- **Dixon-Coles correction in simulation engine** — currently applied to individual match probabilities in the dashboard; extending to the Monte Carlo loop would further improve scoreline-level accuracy

---

## Project Structure

```
wc2026-predictor/
├── R/
│   ├── 01_data_prep.R       # Data loading, filtering, weighting, Elo join
│   ├── 02_model.R           # Poisson GLM + Elo fitting, Dixon-Coles estimation
│   └── 03_simulation.R      # Monte Carlo tournament simulation
├── shiny/
│   └── app.R                # Interactive Shiny dashboard
├── data/                    # Generated CSV outputs
├── plots/                   # Generated visualisations
└── renv.lock                # Reproducible package environment
```

## Reproducing the Analysis

```r
# Install dependencies
renv::restore()

# Download results.csv from Kaggle (Mart Jurisoo) and
# elo_ratings_wc2026.csv from Kaggle (afonsofernandescruz)
# Place both in the data/ folder, then run in order:
source("R/01_data_prep.R")
source("R/02_model.R")
source("R/03_simulation.R")

# Launch dashboard
shiny::runApp("shiny")
```

## Dependencies

Built with R 4.5.2. Key packages: `tidyverse`, `shiny`, `broom`. Full dependency list in `renv.lock`.

---

*Loughborough University · Mathematics with Statistics BSc*
