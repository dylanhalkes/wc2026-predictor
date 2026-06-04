# ============================================================
#  World Cup 2026 Predictor
#  Step 1: Data Loading & Cleaning
# ============================================================
#
#  Data source: Kaggle - "International football results from 1872 to 2024"
#  by Mart Jürisoo. Download results.csv and place in the /data folder.
#  https://www.kaggle.com/datasets/martj42/international-football-results
#
#  What this script does:
#   1. Loads and briefly explores the raw data
#   2. Filters to a relevant date range (2018 onwards)
#   3. Assigns match type weights (World Cup matches > friendlies)
#   4. Applies time decay (recent results matter more)
#   5. Converts to long format ready for Poisson regression
#   6. Saves cleaned output to data/model_data.csv
# ============================================================

library(tidyverse)

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

results_raw <- read_csv("data/results.csv")

cat("── Raw data overview ──\n")
cat("Rows:          ", nrow(results_raw), "\n")
cat("Columns:       ", ncol(results_raw), "\n")
cat("Date range:    ", format(min(results_raw$date)), "to",
                       format(max(results_raw$date)), "\n")
cat("Unique teams:  ", n_distinct(c(results_raw$home_team,
                                    results_raw$away_team)), "\n\n")

glimpse(results_raw)

# ── 2. DATE FILTER ───────────────────────────────────────────────────────────
#
#  We keep from 2018-01-01 onwards. This captures:
#   - The 2018 and 2022 World Cup cycles (qualifying + tournament)
#   - Enough data per team to estimate attack/defence strengths reliably
#   - Not so much old data that it distorts estimates for current squads
#
#  You can experiment with this cutoff later — try 2016 or 2020 and see
#  how your final team strength estimates change.

results_filtered <- results_raw %>%
  filter(date >= as.Date("2018-01-01")) %>%
  mutate(date = as.Date(date))

cat("── After date filter (2018 onwards) ──\n")
cat("Matches remaining:", nrow(results_filtered), "\n\n")

# Preview the tournament types present — useful to check before weighting
cat("── Tournament types in data ──\n")
results_filtered %>%
  count(tournament, sort = TRUE) %>%
  print(n = 40)

# ── 3. MATCH TYPE WEIGHTS ────────────────────────────────────────────────────
#
#  Not all matches are equally informative. A team fielding reserves in a
#  pre-tournament friendly tells us very little about their true strength.
#  We assign a weight to each match type so the model reflects this.
#
#  Weight scale (rough guide):
#   4.0 = FIFA World Cup match (highest signal)
#   3.0 = Continental championship (Euros, Copa América, AFCON etc.)
#   2.5 = World Cup qualifier
#   2.0 = Continental qualifier / other competitive
#   1.5 = UEFA Nations League (competitive but not major tournament)
#   1.0 = Other (catch-all for minor tournaments)
#   0.5 = Friendly (lowest signal)

results_weighted <- results_filtered %>%
  mutate(
    match_weight = case_when(
      tournament == "FIFA World Cup"
      ~ 4.0,
      str_detect(tournament, regex("world cup qualification|world cup qualifier",
                                   ignore_case = TRUE))
      ~ 2.5,
      tournament %in% c("UEFA Euro",
                        "Copa América",
                        "African Cup of Nations",   # fixed
                        "AFC Asian Cup",
                        "Gold Cup",                  # fixed
                        "OFC Nations Cup")
      ~ 3.0,
      str_detect(tournament, regex("qualification|qualifier",
                                   ignore_case = TRUE))
      ~ 2.0,
      str_detect(tournament, regex("nations league",
                                   ignore_case = TRUE))
      ~ 1.5,
      tournament == "Friendly"
      ~ 0.5,
      TRUE
      ~ 1.0
    )
  )

# Sanity check: make sure every tournament got mapped sensibly
cat("\n── Match type weight distribution ──\n")
results_weighted %>%
  count(tournament, match_weight, sort = TRUE) %>%
  print(n = 40)

# ── 4. TIME DECAY ────────────────────────────────────────────────────────────
#
#  Squads, managers, and form all change over time. We want recent results
#  to carry more weight than old ones. We use exponential decay:
#
#    time_weight = exp(-log(2) / half_life * days_ago)
#
#  With a half-life of 730 days (2 years):
#   - A match from today:       weight = 1.00
#   - A match from 1 year ago:  weight ≈ 0.71
#   - A match from 2 years ago: weight ≈ 0.50
#   - A match from 4 years ago: weight ≈ 0.25
#
#  The total weight for each observation is: match_weight × time_weight

TOURNAMENT_START <- as.Date("2026-06-11")
HALF_LIFE_DAYS   <- 730   # experiment with this: try 365 or 1095

results_weighted <- results_weighted %>%
  mutate(
    days_ago     = as.numeric(TOURNAMENT_START - date),
    time_weight  = exp(-log(2) / HALF_LIFE_DAYS * days_ago),
    total_weight = match_weight * time_weight
  )

# Preview weight distribution
cat("\n── Total weight summary ──\n")
summary(results_weighted$total_weight)

# ── 5. CONVERT TO LONG FORMAT ────────────────────────────────────────────────
#
#  The Poisson model we build in step 2 needs data in long format:
#  one row per TEAM per MATCH (so each match becomes two rows).
#
#  Each row tells us: "team X scored Y goals against opponent Z,
#  with home advantage = TRUE/FALSE, with this total weight."
#
#  Note on home advantage: if neutral = TRUE, neither team has home
#  advantage (this is important for tournament matches, which are
#  nearly all played on neutral ground — except for the host nations
#  USA, Canada, and Mexico who may have a partial advantage).

model_data <- bind_rows(

  # Home team perspective
  results_weighted %>%
    transmute(
      date,
      tournament,
      team         = home_team,
      opponent     = away_team,
      goals        = home_score,
      is_home      = !neutral,
      total_weight
    ),

  # Away team perspective
  results_weighted %>%
    transmute(
      date,
      tournament,
      team         = away_team,
      opponent     = home_team,
      goals        = away_score,
      is_home      = FALSE,
      total_weight
    )

) %>%
  filter(!is.na(goals)) %>%   # remove unplayed fixtures
  arrange(date)

# ── 6. SANITY CHECKS ─────────────────────────────────────────────────────────

cat("\n── Final model_data summary ──\n")
cat("Total team-match rows:  ", nrow(model_data), "\n")
cat("Unique teams:           ", n_distinct(model_data$team), "\n")
cat("Date range:             ", format(min(model_data$date)), "to",
                                format(max(model_data$date)), "\n")
cat("Mean goals per row:     ", round(mean(model_data$goals), 3), "\n")
cat("Home advantage rows:    ", sum(model_data$is_home), "\n\n")

# Check a few specific teams to make sure they look right
cat("── Sample: England observations ──\n")
model_data %>%
  filter(team == "England") %>%
  arrange(desc(date)) %>%
  select(date, tournament, team, opponent, goals, is_home, total_weight) %>%
  print(n = 15)

# ── 7. SAVE ──────────────────────────────────────────────────────────────────

write_csv(model_data, "data/model_data.csv")
cat("\n✓ Saved cleaned data to data/model_data.csv\n")
cat("  Rows:", nrow(model_data), "| Columns:", ncol(model_data), "\n")
