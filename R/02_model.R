# ============================================================
#  World Cup 2026 Predictor
#  Step 2: Poisson Goals Model  (v2 - Elo + Dixon-Coles)
# ============================================================
#
#  Model: log(lambda_ij) = mu + alpha_i + beta_j + gamma*home + delta*elo_diff
#
#  Additionally estimates the Dixon-Coles rho correction parameter,
#  which adjusts the probability of 0-0, 0-1, 1-0 and 1-1 scorelines
#  to better match observed football data.
# ============================================================

library(tidyverse)
library(broom)

dir.create("plots", showWarnings = FALSE)

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

model_data <- read_csv("data/model_data.csv", show_col_types = FALSE) %>%
  mutate(
    team     = factor(team),
    opponent = factor(opponent),
    is_home  = as.numeric(is_home)
  )

cat("Loaded:", nrow(model_data), "rows |", nlevels(model_data$team), "teams\n\n")

# ── 2. FIT POISSON GLM WITH ELO ──────────────────────────────────────────────

cat("Fitting Poisson model...\n")

poisson_model <- glm(
  goals ~ team + opponent + is_home + elo_diff,
  data    = model_data,
  family  = poisson(link = "log"),
  weights = total_weight
)

cat("Done.\n\n")

# ── 3. DIAGNOSTICS ───────────────────────────────────────────────────────────

cat("── Diagnostics ──\n")
cat("Residual deviance:", round(poisson_model$deviance, 1),
    "| AIC:", round(AIC(poisson_model), 1), "\n")
cat("Dispersion:",
    round(poisson_model$deviance / poisson_model$df.residual, 3), "\n\n")

# ── 4. EXTRACT COEFFICIENTS ──────────────────────────────────────────────────

coefs <- tidy(poisson_model)

intercept <- coefs %>% filter(term == "(Intercept)") %>% pull(estimate)
home_adv  <- coefs %>% filter(term == "is_home")     %>% pull(estimate)
elo_coef  <- coefs %>% filter(term == "elo_diff")    %>% pull(estimate)

cat("── Key coefficients ──\n")
cat("Intercept mu:              ", round(intercept, 4), "\n")
cat("Home advantage exp(gamma): ", round(exp(home_adv), 3), "x\n")
cat("Elo coefficient delta:     ", round(elo_coef, 4), "\n")
cat("Interpretation: per 100 Elo pt advantage,",
    round((exp(elo_coef) - 1) * 100, 1), "% more expected goals\n\n")

# ── 5. EXTRACT TEAM PARAMETERS ───────────────────────────────────────────────

attack_params <- coefs %>%
  filter(str_starts(term, "team")) %>%
  transmute(team = str_remove(term, "^team"), attack = estimate)

defence_params <- coefs %>%
  filter(str_starts(term, "opponent")) %>%
  transmute(team = str_remove(term, "^opponent"), defence = estimate)

ref_team <- levels(model_data$team)[1]

team_params <- full_join(attack_params, defence_params, by = "team") %>%
  bind_rows(tibble(team = ref_team, attack = 0, defence = 0)) %>%
  mutate(
    attack   = replace_na(attack,  0),
    defence  = replace_na(defence, 0),
    strength = attack - defence
  ) %>%
  arrange(desc(strength))

cat("── Top 15 teams by overall strength ──\n")
team_params %>%
  select(team, attack, defence, strength) %>%
  head(15) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print()

# ── 6. ESTIMATE DIXON-COLES RHO ──────────────────────────────────────────────
#
#  The Poisson model underestimates 0-0 and 1-1 scorelines and
#  overestimates 0-1 and 1-0. The DC correction factor tau adjusts these
#  four low-score cells:
#
#    tau(0,0) = 1 - mu*nu*rho   (increases when rho < 0)
#    tau(0,1) = 1 + mu*rho      (decreases when rho < 0)
#    tau(1,0) = 1 + nu*rho      (decreases when rho < 0)
#    tau(1,1) = 1 - rho         (increases when rho < 0)
#    tau(x,y) = 1               (all other scorelines unchanged)
#
#  We estimate rho by maximising the DC log-likelihood correction,
#  holding the Poisson expected goals fixed.
#
#  Reference: Dixon & Coles (1997), Applied Statistics, 46(2), 265-280.

cat("\n── Estimating Dixon-Coles rho ──\n")

# Tau correction function
tau_dc <- function(x, y, mu, nu, rho) {
  case_when(
    x == 0 & y == 0 ~ 1 - mu * nu * rho,
    x == 0 & y == 1 ~ 1 + mu * rho,
    x == 1 & y == 0 ~ 1 + nu * rho,
    x == 1 & y == 1 ~ 1 - rho,
    TRUE ~ 1
  )
}

# Add Poisson fitted values to model_data
fit_data <- model_data %>%
  mutate(lambda = fitted(poisson_model))

# Reconstruct match-level pairs: need both teams' goals and lambdas
# for the same match to compute tau(goals_a, goals_b, lambda_a, lambda_b)
match_pairs <- fit_data %>%
  select(date, tournament, team, opponent, goals, lambda, total_weight) %>%
  inner_join(
    fit_data %>%
      select(date, tournament,
             team_b = team, opp_b = opponent,
             goals_opp = goals, lambda_opp = lambda),
    by = c("date"       = "date",
           "tournament" = "tournament",
           "team"       = "opp_b",
           "opponent"   = "team_b")
  )

cat("Match pairs reconstructed:", nrow(match_pairs), "\n")

# Negative DC log-likelihood (to minimise over rho)
neg_loglik_dc <- function(rho) {
  taus <- tau_dc(
    match_pairs$goals, match_pairs$goals_opp,
    match_pairs$lambda, match_pairs$lambda_opp,
    rho
  )
  -sum(log(pmax(taus, 1e-10)) * match_pairs$total_weight)
}

dc_opt <- optimize(neg_loglik_dc, interval = c(-0.5, 0))
rho_dc  <- dc_opt$minimum

cat("Estimated rho:", round(rho_dc, 4), "\n")
cat("Interpretation:\n")
cat("  Vs standard Poisson, DC increases P(0-0) and P(1-1)\n")
cat("  and decreases P(0-1) and P(1-0)\n\n")

# Illustrate effect for a representative match (xg_a=1.2, xg_b=1.0)
cat("── DC correction example: xG_a=1.2, xG_b=1.0 ──\n")
demo <- expand.grid(goals_a = 0:3, goals_b = 0:3) %>%
  mutate(
    prob_pois = dpois(goals_a, 1.2) * dpois(goals_b, 1.0),
    tau       = tau_dc(goals_a, goals_b, 1.2, 1.0, rho_dc),
    prob_dc   = prob_pois * tau,
    change    = paste0(ifelse(prob_dc > prob_pois, "+", ""),
                       round((prob_dc - prob_pois) * 100, 2), "pp")
  ) %>%
  filter(goals_a <= 1, goals_b <= 1)
print(demo %>% select(goals_a, goals_b, prob_pois, prob_dc, change))

# ── 7. ADD CURRENT ELO TO WC TEAM PARAMS ────────────────────────────────────

wc_teams    <- read_csv("data/wc_teams.csv",   show_col_types = FALSE)
current_elo <- read_csv("data/current_elo.csv", show_col_types = FALSE)

wc_params <- wc_teams %>%
  left_join(team_params %>% select(team, attack, defence, strength), by = "team") %>%
  left_join(current_elo, by = "team") %>%
  arrange(group, team)

missing_elo <- wc_params %>% filter(is.na(current_elo)) %>% pull(team)
if (length(missing_elo) > 0) {
  median_elo <- median(wc_params$current_elo, na.rm = TRUE)
  wc_params  <- wc_params %>%
    mutate(current_elo = replace_na(current_elo, median_elo))
  cat("Teams with no Elo (assigned median):", paste(missing_elo, collapse=", "), "\n")
} else {
  cat("All 48 WC teams have current Elo. OK\n\n")
}

# ── 8. SPOT CHECKS ───────────────────────────────────────────────────────────

expected_goals <- function(team_a, team_b) {
  a <- wc_params %>% filter(team == team_a)
  b <- wc_params %>% filter(team == team_b)
  xg_a <- exp(intercept + a$attack + b$defence +
               elo_coef * (a$current_elo - b$current_elo) / 100)
  xg_b <- exp(intercept + b$attack + a$defence +
               elo_coef * (b$current_elo - a$current_elo) / 100)
  tibble(team_a=team_a, team_b=team_b, xg_a=round(xg_a,3), xg_b=round(xg_b,3))
}

cat("── Spot checks ──\n")
print(expected_goals("England",   "France"))
print(expected_goals("Argentina", "Brazil"))
print(expected_goals("Morocco",   "Scotland"))

# ── 9. VISUALISE ─────────────────────────────────────────────────────────────

p <- wc_params %>%
  ggplot(aes(x = -defence, y = attack, label = team)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = strength), size = 3, alpha = 0.85) +
  geom_text(size = 2.6, hjust = -0.12, vjust = 0.4, check_overlap = TRUE) +
  scale_colour_gradient2(low = "#3182bd", mid = "#deebf7", high = "#e6550d",
                         midpoint = median(wc_params$strength, na.rm=TRUE),
                         name = "Strength") +
  labs(title    = "2026 World Cup - Team Strengths (Poisson + Elo + Dixon-Coles)",
       subtitle = "Weighted international results 2018-2026",
       x = "Defensive strength (harder to score against)",
       y = "Attack strength (scores more)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plots/team_strengths_v2.png", p, width = 13, height = 8, dpi = 150)

# ── 10. SAVE ─────────────────────────────────────────────────────────────────

write_csv(team_params, "data/team_params.csv")
write_csv(wc_params,   "data/wc_params.csv")

model_scalars <- tibble(
  intercept = intercept,
  home_adv  = home_adv,
  elo_coef  = elo_coef,
  rho_dc    = rho_dc       # Dixon-Coles correction parameter
)
write_csv(model_scalars, "data/model_scalars.csv")

cat("\n-- Saved: data/team_params.csv\n")
cat("-- Saved: data/wc_params.csv\n")
cat("-- Saved: data/model_scalars.csv  (includes rho_dc =", round(rho_dc, 4), ")\n")
cat("-- Saved: plots/team_strengths_v2.png\n")
