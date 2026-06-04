# 2026 FIFA World Cup Predictor

A Monte Carlo tournament simulator for the 2026 FIFA World Cup, built in R. The model estimates each team's attack and defence strengths from historical match data, simulates the full 48-team tournament 50,000 times, and produces win probabilities and stage-by-stage forecasts for every competing nation.

An interactive Shiny dashboard allows users to explore the predictions, simulate individual matchups, and view scoreline probability heat maps.

> **Status:** v1 predictions live (Poisson GLM baseline). Elo ratings and additional predictors being incorporated before the tournament starts on 11 June 2026.

---

## Methodology

### Data
International football results from 1872 to 2026 (Kaggle - Mart Jürisoo), filtered to matches from January 2018 onwards to capture the last two World Cup qualifying cycles. Raw data comprises ~8,000 matches across 284 nations.

### Match Weighting
Not all matches carry equal information. A two-dimensional weighting scheme is applied to each observation:

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

where $d$ is the number of days before the tournament start date. A match from one year ago receives approximately 71% of the weight of a match played today.

The total weight for each observation is the product of both factors.

### Poisson Goals Model
Goals scored are modelled as Poisson distributed with a log-linear conditional mean:

$$\log(\lambda_{ij}) = \mu + \alpha_i + \beta_j + \gamma \cdot \text{home}_i$$

where:
- $\lambda_{ij}$ = expected goals for team $i$ against team $j$
- $\mu$ = intercept (baseline log goal rate)
- $\alpha_i$ = attack strength of team $i$
- $\beta_j$ = defensive weakness of team $j$ (higher = easier to score against)
- $\gamma$ = home advantage coefficient

The model is fitted as a weighted Poisson GLM using base R's `glm()`. With ~16,000 weighted observations and 284 teams, the model estimates approximately 570 parameters.

**Estimated global parameters:**
- Intercept $\mu \approx -0.47$ → baseline goals $e^\mu \approx 0.62$ per match
- Home advantage $e^\gamma \approx 1.27$ × goal multiplier

### Monte Carlo Simulation
The full 48-team tournament is simulated 50,000 times:

1. **Group stage** - all six round-robin matches per group simulated by drawing independently from $\text{Poisson}(\lambda_{ij})$ and $\text{Poisson}(\lambda_{ji})$. Teams ranked by points → goal difference → goals scored → random draw.
2. **Third-place selection** - all 12 third-placed teams ranked; best 8 advance. Rankings use the same tiebreaker hierarchy.
3. **Round of 32** - bracket assembled using FIFA's official Annex C structure. The 8 fixed runner-up matchups are hardcoded; the 8 group-winner-vs-third-place matchups use a constraint-satisfaction algorithm (minimum remaining values) to assign third-place teams to their eligible opponent slots.
4. **Knockout rounds** - match outcomes drawn from Poisson distributions. Draws resolved by a penalty shootout, modelled as a Bernoulli trial with success probability proportional to each team's expected goals in that match.
5. **Results aggregated** - win, finalist, semi-final, quarter-final, R16, and R32 probabilities computed across all 50,000 iterations.

Standard error of probability estimates at 50,000 simulations: ±0.18 percentage points for a team with ~20% win probability.

---

## Results (v1 - Poisson GLM baseline)

Pre-tournament predictions published **4 June 2026**, one week before kick-off.

| Rank | Team | Win % | Final % | Semi % | QF % |
|---|---|---|---|---|---|
| 1 | Spain | 19.7% | 29.8% | 43.9% | 55.3% |
| 2 | England | 14.1% | 23.9% | 36.9% | 55.5% |
| 3 | Argentina | 12.9% | 20.9% | 34.3% | 47.8% |
| 4 | France | 7.9% | 15.5% | 26.3% | 45.8% |
| 5 | Brazil | 7.0% | 14.6% | 27.0% | 40.8% |
| 6 | Portugal | 5.8% | 11.3% | 23.0% | 42.0% |
| 7 | Germany | 5.1% | 11.7% | 23.7% | 46.7% |
| 8 | Morocco | 4.1% | 9.4% | 18.7% | 31.3% |
| 9 | Colombia | 3.6% | 7.9% | 17.1% | 33.1% |
| 10 | Netherlands | 3.1% | 7.3% | 15.6% | 26.4% |

**Notable findings:**
- Spain's model-implied dominance reflects both exceptional attack (highest $\alpha$ among WC teams) and elite defensive solidity, combined with a favourable Group H bracket path
- England's 14.1% win probability is notably high relative to their FIFA ranking; the model rewards their consistent high-scoring recent qualifying campaign
- Germany and Belgium show high quarter-final probabilities (46.7% and 40.8%) but comparatively low win probabilities (5.1% and 2.9%), reflecting bracket paths that route them into stronger opposition in the later rounds
- France's win probability (7.9%) sits below their FIFA world ranking (#1), as the model weights recent performance over accumulated ranking points - Spain's recent form has been more dominant

---

## Planned Improvements (v2)

- **Elo ratings prior** - incorporating eloratings.net team strength ratings to improve parameter estimates for data-sparse nations (Jordan, Haiti, Curaçao) where the current model has limited historical match data
- **Dixon-Coles correction** - low-score adjustment to better capture the frequency of 0-0 and 1-0 scorelines, which the standard Poisson model slightly underestimates
- **Calibration analysis** - post-tournament comparison of predicted vs actual outcomes across all 104 matches

---

## Project Structure

```
wc2026-predictor/
├── R/
│   ├── 01_data_prep.R       # Data loading, filtering, weighting
│   ├── 02_model.R           # Poisson GLM fitting, parameter extraction
│   └── 03_simulation.R      # Monte Carlo tournament simulation
├── shiny/
│   └── app.R                # Interactive Shiny dashboard
├── data/                    # Generated CSV outputs (gitignored raw data)
├── plots/                   # Generated visualisations
└── renv.lock                # Reproducible package environment
```

## Reproducing the Analysis

```r
# Install dependencies
renv::restore()

# Run pipeline in order
source("R/01_data_prep.R")   # requires data/results.csv from Kaggle
source("R/02_model.R")
source("R/03_simulation.R")

# Launch dashboard
shiny::runApp("shiny")
```

Raw data (`results.csv`) must be downloaded separately from [Kaggle](https://www.kaggle.com/datasets/martj42/international-football-results).

---

## Dependencies

Built with R 4.5.2. Key packages: `tidyverse`, `shiny`, `broom`. Full dependency list in `renv.lock`.

---

*Loughborough University · Mathematics with Statistics BSc*
