# ============================================================
#  World Cup 2026 Predictor
#  Step 2: Poisson Goals Model  (v2 - with Elo ratings)
# ============================================================
#
#  Model: log(lambda_ij) = mu + alpha_i + beta_j + gamma*home + delta*elo_diff
#
#  elo_diff = (team_elo - opp_elo) / 100  at the year before the match.
#  Non-WC opponents use median Elo (elo_diff ~= 0 for those matches).
#  delta tells us: per 100 Elo point advantage, how many more goals
#  are expected, over and above what the fixed effects already capture?
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

cat("Loaded:", nrow(model_data), "rows |", nlevels(model_data$team), "teams\n")
cat("elo_diff column present:", "elo_diff" %in% colnames(model_data), "\n\n")

# ── 2. FIT POISSON GLM WITH ELO ──────────────────────────────────────────────

cat("Fitting Poisson model with Elo predictor...\n")

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
cat("Intercept mu:               ", round(intercept, 4), "\n")
cat("Home advantage exp(gamma):  ", round(exp(home_adv), 3), "x\n")
cat("Elo coefficient delta:      ", round(elo_coef, 4), "\n")
cat("Interpretation: per 100 Elo point advantage,",
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

# ── 6. ADD CURRENT ELO TO WC TEAM PARAMS ────────────────────────────────────

wc_teams   <- read_csv("data/wc_teams.csv",   show_col_types = FALSE)
current_elo <- read_csv("data/current_elo.csv", show_col_types = FALSE)

wc_params <- wc_teams %>%
  left_join(team_params %>% select(team, attack, defence, strength), by = "team") %>%
  left_join(current_elo, by = "team") %>%
  arrange(group, team)

missing_params <- wc_params %>% filter(is.na(attack)) %>% pull(team)
missing_elo    <- wc_params %>% filter(is.na(current_elo)) %>% pull(team)

if (length(missing_params) > 0)
  cat("WARNING - missing model params:", paste(missing_params, collapse = ", "), "\n")
if (length(missing_elo) > 0) {
  cat("Teams with no current Elo (check name):", paste(missing_elo, collapse = ", "), "\n")
  median_elo <- median(wc_params$current_elo, na.rm = TRUE)
  wc_params <- wc_params %>% mutate(current_elo = replace_na(current_elo, median_elo))
} else {
  cat("All 48 WC teams have current Elo. OK\n")
}

cat("\n── WC teams: current Elo (top 10) ──\n")
wc_params %>%
  arrange(desc(current_elo)) %>%
  select(group, team, attack, defence, current_elo) %>%
  head(10) %>%
  print()

# ── 7. SPOT CHECKS ───────────────────────────────────────────────────────────

expected_goals <- function(team_a, team_b) {
  a <- wc_params %>% filter(team == team_a)
  b <- wc_params %>% filter(team == team_b)
  xg_a <- exp(intercept + a$attack + b$defence +
               elo_coef * (a$current_elo - b$current_elo) / 100)
  xg_b <- exp(intercept + b$attack + a$defence +
               elo_coef * (b$current_elo - a$current_elo) / 100)
  tibble(team_a = team_a, team_b = team_b,
         xg_a = round(xg_a, 3), xg_b = round(xg_b, 3))
}

cat("\n── Spot checks (neutral ground, with Elo) ──\n")
print(expected_goals("England",   "France"))
print(expected_goals("Argentina", "Brazil"))
print(expected_goals("Germany",   "Spain"))
print(expected_goals("Morocco",   "Scotland"))

# ── 8. VISUALISE ─────────────────────────────────────────────────────────────

p <- wc_params %>%
  ggplot(aes(x = -defence, y = attack, label = team)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = strength), size = 3, alpha = 0.85) +
  geom_text(size = 2.6, hjust = -0.12, vjust = 0.4, check_overlap = TRUE) +
  scale_colour_gradient2(low = "#3182bd", mid = "#deebf7", high = "#e6550d",
                         midpoint = median(wc_params$strength, na.rm = TRUE),
                         name = "Strength") +
  labs(title    = "2026 World Cup - Team Strengths (v2: Poisson + Elo)",
       subtitle = "Poisson GLM with Elo covariate, weighted 2018-2026",
       x = "Defensive strength (harder to score against)",
       y = "Attack strength (scores more)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p)
ggsave("plots/team_strengths_v2.png", p, width = 13, height = 8, dpi = 150)

# ── 9. SAVE ──────────────────────────────────────────────────────────────────

write_csv(team_params, "data/team_params.csv")
write_csv(wc_params,   "data/wc_params.csv")
write_csv(
  tibble(intercept = intercept, home_adv = home_adv, elo_coef = elo_coef),
  "data/model_scalars.csv"
)

cat("\n✓ Saved: data/team_params.csv\n")
cat("✓ Saved: data/wc_params.csv\n")
cat("✓ Saved: data/model_scalars.csv  (includes elo_coef)\n")
cat("✓ Saved: plots/team_strengths_v2.png\n")
