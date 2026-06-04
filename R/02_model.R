# ============================================================
#  World Cup 2026 Predictor
#  Step 2: Poisson Goals Model
# ============================================================
#
#  We model goals scored as Poisson distributed with a log-linear mean:
#
#    log(λ_ij) = μ + α_i + β_j + γ * home_i
#
#  where:
#    λ_ij = expected goals for team i against team j
#    μ    = intercept (log of baseline goal rate)
#    α_i  = attack strength of team i  (higher → scores more)
#    β_j  = defensive weakness of team j  (higher → concedes more)
#    γ    = home advantage coefficient
#
#  Goals ~ Poisson(λ_ij), fitted as a weighted GLM.
#  All parameters are relative to a reference team (first alphabetically),
#  whose attack and defence are fixed at 0 for identifiability.
#
#  Reference: Dixon & Coles (1997), Applied Statistics, 46(2), 265–280.
# ============================================================

library(tidyverse)
library(broom)

# Create plots directory if it doesn't exist
dir.create("plots", showWarnings = FALSE)

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

model_data <- read_csv("data/model_data.csv") %>%
  mutate(
    team     = factor(team),
    opponent = factor(opponent),
    is_home  = as.numeric(is_home)
  )

cat("── Data loaded ──\n")
cat("Rows:         ", nrow(model_data), "\n")
cat("Unique teams: ", nlevels(model_data$team), "\n")
cat("Reference team (α = β = 0):", levels(model_data$team)[1], "\n\n")

# ── 2. FIT POISSON GLM ───────────────────────────────────────────────────────
#
#  glm() with family = poisson(link = "log") fits the model above.
#  observations are weighted by total_weight (match type × time decay).
#
#  This produces ~560 parameters (attack + defence for ~280 teams).
#  Expect a few seconds to run.

cat("Fitting Poisson model — this may take a moment...\n")

poisson_model <- glm(
  goals ~ team + opponent + is_home,
  data    = model_data,
  family  = poisson(link = "log"),
  weights = total_weight
)

cat("Done.\n\n")

# ── 3. MODEL DIAGNOSTICS ─────────────────────────────────────────────────────

cat("── Model diagnostics ──\n")
cat("Null deviance:     ", round(poisson_model$null.deviance, 1),
    " (df =", poisson_model$df.null, ")\n")
cat("Residual deviance: ", round(poisson_model$deviance, 1),
    " (df =", poisson_model$df.residual, ")\n")
cat("AIC:               ", round(AIC(poisson_model), 1), "\n")

# Dispersion check: residual deviance / df should be ~1 for a well-fitting
# Poisson model. Values >1 indicate overdispersion (common in football data).
dispersion <- poisson_model$deviance / poisson_model$df.residual
cat("Dispersion:        ", round(dispersion, 3),
    ifelse(dispersion > 1.5, " ← some overdispersion (expected)\n", " ← good\n"))

# ── 4. EXTRACT PARAMETERS ────────────────────────────────────────────────────

coefs <- tidy(poisson_model)

# Global scalars
intercept <- coefs %>% filter(term == "(Intercept)") %>% pull(estimate)
home_adv  <- coefs %>% filter(term == "is_home")     %>% pull(estimate)

cat("\n── Key scalars ──\n")
cat("Intercept μ:          ", round(intercept, 4), "\n")
cat("Baseline goals exp(μ):", round(exp(intercept), 3),
    "(goals per match for reference-vs-reference)\n")
cat("Home advantage γ:     ", round(home_adv, 4), "\n")
cat("Home advantage exp(γ):", round(exp(home_adv), 3),
    "× goal multiplier at home\n\n")

# Attack parameters: team_X coefficient → α_X
attack_params <- coefs %>%
  filter(str_starts(term, "team")) %>%
  transmute(
    team   = str_remove(term, "^team"),
    attack = estimate
  )

# Defence parameters: opponent_X coefficient → β_X
# IMPORTANT: β_X is the DEFENSIVE WEAKNESS of team X:
#   high β_X → easy to score against (weak defence)
#   low β_X  → hard to score against (strong defence)
defence_params <- coefs %>%
  filter(str_starts(term, "opponent")) %>%
  transmute(
    team    = str_remove(term, "^opponent"),
    defence = estimate
  )

# ── 5. ASSEMBLE TEAM PARAMETERS ──────────────────────────────────────────────

ref_team <- levels(model_data$team)[1]

team_params <- full_join(attack_params, defence_params, by = "team") %>%
  bind_rows(tibble(team = ref_team, attack = 0, defence = 0)) %>%
  mutate(
    attack   = replace_na(attack, 0),
    defence  = replace_na(defence, 0),
    # Overall strength: good attack (high α) + strong defence (low β)
    strength = attack - defence
  ) %>%
  arrange(desc(strength))

cat("── Top 20 teams by overall strength ──\n")
team_params %>%
  select(team, attack, defence, strength) %>%
  head(20) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print()

# ── 6. EXPECTED GOALS FUNCTION ───────────────────────────────────────────────
#
#  For a match between team_a and team_b:
#    λ_a = exp(μ + α_a + β_b + γ * home_a)
#    λ_b = exp(μ + α_b + β_a + γ * home_b)
#
#  For World Cup matches: home_a = home_b = FALSE (neutral ground)
#  Exception: USA, Canada, Mexico may have a small home advantage.

expected_goals <- function(team_a, team_b,
                           home_a = FALSE, home_b = FALSE) {
  a <- team_params %>% filter(team == team_a)
  b <- team_params %>% filter(team == team_b)

  if (nrow(a) == 0) stop(paste("Team not found:", team_a))
  if (nrow(b) == 0) stop(paste("Team not found:", team_b))

  lambda_a <- exp(intercept + a$attack + b$defence + home_adv * home_a)
  lambda_b <- exp(intercept + b$attack + a$defence + home_adv * home_b)

  tibble(
    team_a   = team_a,
    team_b   = team_b,
    xg_a     = round(lambda_a, 3),
    xg_b     = round(lambda_b, 3)
  )
}

# ── 7. SPOT CHECKS ───────────────────────────────────────────────────────────

cat("\n── Expected goals spot checks (neutral ground) ──\n")
print(expected_goals("England",   "France"))
print(expected_goals("Argentina", "Brazil"))
print(expected_goals("Germany",   "Spain"))
print(expected_goals("Morocco",   "Scotland"))

# ── 8. MAP WORLD CUP TEAM NAMES ──────────────────────────────────────────────
#
#  The Kaggle dataset uses its own naming conventions which don't always
#  match official FIFA names. We create a lookup table mapping each WC team
#  to its Kaggle name. The script will flag any mismatches so you can fix them.
#
#  GROUPS (from the December 2025 draw):
#   A: Mexico, South Korea, South Africa, Czechia
#   B: Canada, Switzerland, Qatar, Bosnia and Herzegovina
#   C: Brazil, Morocco, Scotland, Haiti
#   D: United States, Paraguay, Australia, Turkey
#   E: Germany, Ecuador, Ivory Coast, Curaçao
#   F: Netherlands, Japan, Tunisia, Sweden
#   G: Belgium, Egypt, Iran, New Zealand
#   H: Spain, Uruguay, Saudi Arabia, Cape Verde
#   I: France, Senegal, Norway, Iraq
#   J: Argentina, Algeria, Austria, Jordan
#   K: Portugal, Colombia, Uzbekistan, DR Congo
#   L: England, Croatia, Ghana, Panama

wc_teams <- tribble(
  ~group, ~team,
  # Group A
  "A", "Mexico",
  "A", "South Korea",
  "A", "South Africa",
  "A", "Czech Republic",         # Kaggle uses Czech Republic, not Czechia
  # Group B
  "B", "Canada",
  "B", "Switzerland",
  "B", "Qatar",
  "B", "Bosnia and Herzegovina", # CHECK: may be "Bosnia-Herzegovina" in data
  # Group C
  "C", "Brazil",
  "C", "Morocco",
  "C", "Scotland",
  "C", "Haiti",
  # Group D
  "D", "United States",
  "D", "Paraguay",
  "D", "Australia",
  "D", "Turkey",                 # Kaggle uses Turkey, not Türkiye
  # Group E
  "E", "Germany",
  "E", "Ecuador",
  "E", "Ivory Coast",
  "E", "Curaçao",                # CHECK: may be "Curaçao" in data
  # Group F
  "F", "Netherlands",
  "F", "Japan",
  "F", "Tunisia",
  "F", "Sweden",
  # Group G
  "G", "Belgium",
  "G", "Egypt",
  "G", "Iran",
  "G", "New Zealand",
  # Group H
  "H", "Spain",
  "H", "Uruguay",
  "H", "Saudi Arabia",
  "H", "Cape Verde",
  # Group I
  "I", "France",
  "I", "Senegal",
  "I", "Norway",
  "I", "Iraq",
  # Group J
  "J", "Argentina",
  "J", "Algeria",
  "J", "Austria",
  "J", "Jordan",
  # Group K
  "K", "Portugal",
  "K", "Colombia",
  "K", "Uzbekistan",
  "K", "DR Congo",               # CHECK: may be "Congo DR" in data
  # Group L
  "L", "England",
  "L", "Croatia",
  "L", "Ghana",
  "L", "Panama"
)

# Check which team names are found in our model
cat("\n── World Cup team name check ──\n")
wc_teams <- wc_teams %>%
  mutate(in_model = team %in% team_params$team)

found   <- sum(wc_teams$in_model)
missing <- wc_teams %>% filter(!in_model) %>% pull(team)

cat("Found in model:", found, "/ 48\n")

if (length(missing) > 0) {
  cat("\nNOT FOUND — fix the spelling to match the Kaggle dataset:\n")
  cat(paste(" ✗", missing, collapse = "\n"), "\n")
  cat("\nRun this to find the correct spelling:\n")
  cat('  team_params %>% filter(str_detect(team, "Bosnia|Curacao|Congo")) %>% pull(team)\n')
} else {
  cat("All 48 teams found! ✓\n")
}

# ── 9. SAVE OUTPUTS ──────────────────────────────────────────────────────────

write_csv(team_params, "data/team_params.csv")

model_scalars <- tibble(intercept = intercept, home_adv = home_adv)
write_csv(model_scalars, "data/model_scalars.csv")

write_csv(wc_teams, "data/wc_teams.csv")

cat("\n✓ Saved: data/team_params.csv\n")
cat("✓ Saved: data/model_scalars.csv\n")
cat("✓ Saved: data/wc_teams.csv\n")

# ── 10. VISUALISE ────────────────────────────────────────────────────────────
#
#  Attack (α) vs Defence strength (-β) scatter for WC teams.
#  We flip the defence axis so "better" = right and up on both axes:
#   top-right  = elite (strong attack AND defence)  e.g. France, Spain
#   top-left   = attack-heavy (scores a lot, concedes a lot)
#   bottom-right = defensive (hard to score against, doesn't score much)
#   bottom-left  = weaker teams

wc_params <- team_params %>%
  filter(team %in% wc_teams$team) %>%
  left_join(wc_teams, by = "team")

p <- wc_params %>%
  ggplot(aes(x = -defence, y = attack, label = team)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", alpha = 0.6) +
  geom_point(aes(colour = strength), size = 3, alpha = 0.85) +
  geom_text(size = 2.6, hjust = -0.12, vjust = 0.4, check_overlap = TRUE) +
  scale_colour_gradient2(
    low      = "#3182bd",
    mid      = "#deebf7",
    high     = "#e6550d",
    midpoint = median(wc_params$strength),
    name     = "Overall\nstrength"
  ) +
  labs(
    title    = "2026 World Cup — Team Strengths",
    subtitle = "Estimated from weighted international results (2018–2026) via Poisson GLM",
    x        = "Defensive strength  →  (harder to score against)",
    y        = "Attack strength  ↑  (scores more)",
    caption  = "Parameters relative to baseline | Both axes: higher = better"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    plot.caption  = element_text(colour = "grey50", size = 8),
    legend.position = "right"
  )

print(p)
ggsave("plots/team_strengths.png", p, width = 13, height = 8, dpi = 150)
cat("\n✓ Saved: plots/team_strengths.png\n")
