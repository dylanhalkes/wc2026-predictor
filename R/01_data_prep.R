# ============================================================
#  World Cup 2026 Predictor
#  Step 1: Data Loading & Cleaning  (v2 - with Elo ratings)
# ============================================================

library(tidyverse)

# ── 1. LOAD RAW RESULTS ──────────────────────────────────────────────────────

results_raw <- read_csv("data/results.csv", show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

cat("── Raw data overview ──\n")
cat("Rows:", nrow(results_raw), "| Date range:",
    format(min(results_raw$date)), "to", format(max(results_raw$date)), "\n\n")

# ── 2. DATE FILTER ───────────────────────────────────────────────────────────

results_filtered <- results_raw %>%
  filter(date >= as.Date("2018-01-01"))

cat("Matches from 2018 onwards:", nrow(results_filtered), "\n\n")

# ── 3. MATCH TYPE WEIGHTS ────────────────────────────────────────────────────

results_weighted <- results_filtered %>%
  mutate(
    match_weight = case_when(
      tournament == "FIFA World Cup"
        ~ 4.0,
      str_detect(tournament, regex("world cup qualification|world cup qualifier",
                                   ignore_case = TRUE))
        ~ 2.5,
      tournament %in% c("UEFA Euro", "Copa America", "African Cup of Nations",
                        "AFC Asian Cup", "Gold Cup", "OFC Nations Cup")
        ~ 3.0,
      str_detect(tournament, regex("qualification|qualifier", ignore_case = TRUE))
        ~ 2.0,
      str_detect(tournament, regex("nations league", ignore_case = TRUE))
        ~ 1.5,
      tournament == "Friendly"
        ~ 0.5,
      TRUE ~ 1.0
    )
  )

# ── 4. TIME DECAY ────────────────────────────────────────────────────────────

TOURNAMENT_START <- as.Date("2026-06-11")
HALF_LIFE_DAYS   <- 730

results_weighted <- results_weighted %>%
  mutate(
    days_ago     = as.numeric(TOURNAMENT_START - date),
    time_weight  = exp(-log(2) / HALF_LIFE_DAYS * days_ago),
    total_weight = match_weight * time_weight
  )

# ── 5. CONVERT TO LONG FORMAT ────────────────────────────────────────────────

model_data <- bind_rows(
  results_weighted %>%
    transmute(date, tournament, team = home_team, opponent = away_team,
              goals = home_score, is_home = !neutral, total_weight),
  results_weighted %>%
    transmute(date, tournament, team = away_team, opponent = home_team,
              goals = away_score, is_home = FALSE, total_weight)
) %>%
  filter(!is.na(goals)) %>%
  arrange(date)

cat("Model data rows:", nrow(model_data),
    "| Mean goals:", round(mean(model_data$goals), 3), "\n\n")

# ── 6. JOIN ELO RATINGS ──────────────────────────────────────────────────────
#
#  Source: elo_ratings_wc2026.csv (Kaggle - afonsofernandescruz)
#  Coverage: 48 WC teams only, annual year-end snapshots 1901-2026
#            plus one live pre-tournament snapshot per team.
#
#  Strategy: for each match in year Y, look up each team's year-end
#  Elo from year Y-1. Non-WC opponents (not in the 48) are assigned
#  the median Elo (~1650) so their contribution to elo_diff is neutral.
#
#  For the simulation all 48 WC teams have Elo, so coverage is perfect
#  where it matters most.

cat("Loading Elo ratings...\n")

elo_raw <- read_csv("data/elo_ratings_wc2026.csv", show_col_types = FALSE)

cat("Elo file columns:", paste(colnames(elo_raw), collapse = ", "), "\n")
cat("Elo file rows:", nrow(elo_raw), "\n\n")

# Name mapping: Elo dataset uses modern/official names; results.csv uses
# traditional English names. Map Elo names -> results.csv names.
elo_name_map <- c(
  "Czechia"              = "Czech Republic",
  "United States"        = "United States",          # same
  "South Korea"          = "South Korea",            # same
  "DR Congo"             = "DR Congo",               # same
  "Bosnia and Herzegovina" = "Bosnia and Herzegovina", # same
  "Turkey"               = "Turkey",                 # same
  "Ivory Coast"          = "Ivory Coast",            # same
  "Curacao"              = "Curacao"                 # same
)

# Apply name mapping (only Czechia actually differs)
elo_raw <- elo_raw %>%
  mutate(team = ifelse(country %in% names(elo_name_map),
                       elo_name_map[country],
                       country))

# Year-end snapshots only (snapshot_date ends in "12-31")
elo_annual <- elo_raw %>%
  filter(str_detect(as.character(snapshot_date), "12-31")) %>%
  select(team, year, rating) %>%
  arrange(team, year)

cat("Year-end Elo snapshots:", nrow(elo_annual), "\n")
cat("Year range:", min(elo_annual$year), "to", max(elo_annual$year), "\n\n")

# Current (pre-tournament) Elo: most recent snapshot per team
current_elo <- elo_raw %>%
  group_by(team) %>%
  slice_max(as.Date(as.character(snapshot_date)), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(team, current_elo = rating)

cat("Current Elo loaded for", nrow(current_elo), "teams.\n\n")

# For each match in year Y, use year Y-1 Elo
# (year-end Elo is a pre-season proxy for the following year's matches)
model_data <- model_data %>%
  mutate(elo_lookup_year = as.integer(format(date, "%Y")) - 1L)

# Join team Elo
model_data <- model_data %>%
  left_join(
    elo_annual %>% rename(elo_team = rating),
    by = c("team" = "team", "elo_lookup_year" = "year")
  )

# Join opponent Elo
model_data <- model_data %>%
  left_join(
    elo_annual %>% rename(elo_opp = rating),
    by = c("opponent" = "team", "elo_lookup_year" = "year")
  )

# Non-WC teams get median Elo (neutral effect on elo_diff)
median_elo <- median(elo_annual$rating, na.rm = TRUE)
cat("Median Elo (used for non-WC opponents):", round(median_elo), "\n")

n_missing_team <- sum(is.na(model_data$elo_team))
n_missing_opp  <- sum(is.na(model_data$elo_opp))
cat("Rows missing team Elo:", n_missing_team, "(", round(n_missing_team/nrow(model_data)*100,1), "%)\n")
cat("Rows missing opp Elo: ", n_missing_opp,  "(", round(n_missing_opp/nrow(model_data)*100,1), "%)\n\n")

model_data <- model_data %>%
  mutate(
    elo_team = replace_na(elo_team, median_elo),
    elo_opp  = replace_na(elo_opp,  median_elo),
    elo_diff = (elo_team - elo_opp) / 100
  ) %>%
  select(-elo_lookup_year)

cat("── Elo difference summary (per 100 Elo points) ──\n")
print(summary(model_data$elo_diff))

cat("\n── England sample (most recent 8 matches) ──\n")
model_data %>%
  filter(team == "England") %>%
  arrange(desc(date)) %>%
  select(date, opponent, goals, elo_team, elo_opp, elo_diff) %>%
  head(8) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  print()

# ── 7. SAVE ──────────────────────────────────────────────────────────────────

write_csv(model_data,  "data/model_data.csv")
write_csv(current_elo, "data/current_elo.csv")

cat("\n✓ Saved data/model_data.csv  |", nrow(model_data), "rows\n")
cat("✓ Saved data/current_elo.csv |", nrow(current_elo), "teams\n")
