# ============================================================
#  World Cup 2026 Predictor
#  Step 3: Monte Carlo Tournament Simulation
# ============================================================
#
#  We simulate the full 48-team tournament N_SIM times and aggregate
#  the results into probabilities for each team at each stage.
#
#  Each simulation run:
#   1. Group stage — 6 Poisson matches per group (× 12 groups)
#   2. Third-place ranking — best 8 of 12 third-placed teams advance
#   3. Round of 32 — built from group outcomes using the official
#      FIFA Annex C bracket (fixed matchups + eligibility constraints)
#   4. Knockouts — draws resolved by a Poisson-weighted penalty coin flip
#
#  R16 onward uses sequential pairing of adjacent R32 winners.
#
#  Run time: ~30–60 seconds for N_SIM = 10,000 on a modern laptop.
#  Set N_SIM = 50000 for tighter estimates (~3–5 minutes).
# ============================================================

library(tidyverse)
set.seed(2026)

N_SIM <- 50000

# ── 1. LOAD & SETUP ──────────────────────────────────────────────────────────

team_params   <- read_csv("data/team_params.csv",   show_col_types = FALSE)
model_scalars <- read_csv("data/model_scalars.csv", show_col_types = FALSE)
wc_teams      <- read_csv("data/wc_teams.csv",      show_col_types = FALSE)

intercept <- model_scalars$intercept
home_adv  <- model_scalars$home_adv

# Join group assignments with Poisson model parameters
wc_params <- wc_teams %>%
  left_join(team_params %>% select(team, attack, defence), by = "team") %>%
  arrange(group, team)

# Sanity check: all teams have model parameters
if (any(is.na(wc_params$attack))) {
  missing <- wc_params %>% filter(is.na(attack)) %>% pull(team)
  stop("Missing model parameters for: ", paste(missing, collapse = ", "),
       "\nCheck team name spelling in wc_teams.csv vs team_params.csv.")
}

TEAMS         <- wc_params$team
GROUPS        <- wc_params$group
N             <- nrow(wc_params)   # 48
GROUP_LETTERS <- LETTERS[1:12]     # A through L

# Team indices for each group (e.g. groups_idx$A = c(1,2,3,4))
groups_idx <- setNames(
  lapply(GROUP_LETTERS, function(g) which(GROUPS == g)),
  GROUP_LETTERS
)

# Pre-compute expected goals matrix for all 48×48 matchups (neutral ground)
# lambda_mat[i,j] = E[goals by team i against team j]
# Formula: exp(μ + α_i + β_j)  where μ = intercept, α = attack, β = defence
lambda_mat <- exp(
  outer(wc_params$attack, wc_params$defence, "+") + intercept
)
diag(lambda_mat) <- 0  # team never plays itself

cat("Setup complete.\n")
cat("Teams:", N, "| Pre-computed", N*N, "expected goals values.\n\n")

# ── 2. ELIGIBLE GROUPS FOR THIRD-PLACE BRACKET ───────────────────────────────
#
#  From FIFA Annex C regulations: 8 group winners face 3rd-place teams
#  in the Round of 32. Each group winner has an eligible pool of source groups.
#
#  Source: 2026 FIFA World Cup regulations, Annex C (combinations table)
#
#  Keys  = group winner slot (A, B, D, E, G, I, K, L)
#  Values = groups from which their 3rd-place opponent can come

eligible_3rd <- list(
  A = c("C","E","F","H","I"),
  B = c("E","F","G","I","J"),
  D = c("B","E","F","I","J"),
  E = c("A","B","C","D","F"),
  G = c("A","E","H","I","J"),
  I = c("C","D","F","G","H"),
  K = c("D","E","I","J","L"),
  L = c("E","H","I","J","K")
)

# Assign 8 qualifying 3rd-place groups to 8 group winner slots.
# Uses minimum-remaining-values (MRV) heuristic to find a valid matching
# consistent with FIFA's eligibility constraints.

assign_3rd_places <- function(qualifying_groups) {
  slots    <- names(eligible_3rd)
  assigned <- setNames(rep(NA_character_, length(slots)), slots)
  remain   <- qualifying_groups

  for (iter in seq_along(slots)) {
    # Pick the most constrained unassigned slot
    n_avail <- sapply(slots, function(s) {
      if (!is.na(assigned[s])) Inf
      else length(intersect(eligible_3rd[[s]], remain))
    })
    slot  <- names(which.min(n_avail))
    avail <- intersect(eligible_3rd[[slot]], remain)

    chosen <- if (length(avail) > 0) sample(avail, 1) else remain[1]

    assigned[slot] <- chosen
    remain <- setdiff(remain, chosen)
  }
  assigned
}

# ── 3. CORE SIMULATION FUNCTIONS ─────────────────────────────────────────────

# Simulate one group (4 teams, all 6 round-robin matches).
# Returns a list: idx (ranked team indices), pts, gd, gf.
sim_group <- function(teams4) {
  pts <- gd <- gf <- integer(4L)

  for (ii in 1L:3L) {
    for (jj in (ii + 1L):4L) {
      gi <- rpois(1L, lambda_mat[teams4[ii], teams4[jj]])
      gj <- rpois(1L, lambda_mat[teams4[jj], teams4[ii]])

      gf[ii] <- gf[ii] + gi;  gf[jj] <- gf[jj] + gj
      gd[ii] <- gd[ii] + gi - gj;  gd[jj] <- gd[jj] + gj - gi

      if      (gi > gj) { pts[ii] <- pts[ii] + 3L }
      else if (gj > gi) { pts[jj] <- pts[jj] + 3L }
      else { pts[ii] <- pts[ii] + 1L; pts[jj] <- pts[jj] + 1L }
    }
  }

  ord <- order(-pts, -gd, -gf, runif(4L))
  list(idx = teams4[ord], pts = pts[ord], gd = gd[ord], gf = gf[ord])
}

# Simulate one knockout match. Draws resolved by Poisson-weighted penalty flip.
sim_ko <- function(i, j) {
  gi <- rpois(1L, lambda_mat[i, j])
  gj <- rpois(1L, lambda_mat[j, i])

  if      (gi > gj) return(i)
  else if (gj > gi) return(j)
  else {
    p_i <- lambda_mat[i, j] / (lambda_mat[i, j] + lambda_mat[j, i])
    if (runif(1L) < p_i) i else j
  }
}

# ── 4. FULL TOURNAMENT SIMULATION ────────────────────────────────────────────
#
#  Returns an integer vector of length N (one per team) where the value
#  represents the furthest round the team was eliminated at:
#
#   0 = group stage exit (or 3rd-place team that didn't advance)
#   1 = Round of 32 exit
#   2 = Round of 16 exit
#   3 = Quarter-final exit
#   4 = Semi-final exit
#   5 = Runner-up (lost the final)
#   6 = Winner

simulate_tournament <- function() {

  rr <- rep(0L, N)  # round_reached: updated as teams are eliminated

  # ── Group stage ────────────────────────────────────────────────────────────
  pos1 <- pos2 <- pos3 <- pos4    <- integer(12L)
  t3_pts <- t3_gd <- t3_gf        <- integer(12L)

  for (g in seq_len(12L)) {
    s <- sim_group(groups_idx[[g]])
    pos1[g] <- s$idx[1L]; pos2[g] <- s$idx[2L]
    pos3[g] <- s$idx[3L]; pos4[g] <- s$idx[4L]
    t3_pts[g] <- s$pts[3L]; t3_gd[g] <- s$gd[3L]; t3_gf[g] <- s$gf[3L]
  }

  # ── 8 best 3rd-placed teams ────────────────────────────────────────────────
  #  Ranked by: points → GD → GF → random tiebreak
  t3_rank  <- order(-t3_pts, -t3_gd, -t3_gf, runif(12L))
  adv_idx  <- t3_rank[1L:8L]              # which groups' 3rd place advances
  adv_grps <- GROUP_LETTERS[adv_idx]      # their group letters

  # ── R32 bracket ────────────────────────────────────────────────────────────
  slot_grp <- assign_3rd_places(adv_grps) # slot_grp["A"] = group letter of 1A's 3rd-place opp

  get_3rd <- function(slot) {
    g <- slot_grp[slot]
    if (is.na(g)) return(pos3[adv_idx[1L]])  # fallback (shouldn't fire)
    pos3[which(GROUP_LETTERS == g)]
  }

  # 16 R32 matches: official FIFA bracket (matches M73–M88)
  # Fixed (2nd vs 2nd):
  #   M73 2A–2B | M75 1F–2C | M76 1C–2F | M78 2E–2I
  #   M83 2K–2L | M84 1H–2J | M86 1J–2H | M88 2D–2G
  # Variable (1st vs best eligible 3rd):
  #   M74 1E–3? | M77 1I–3? | M79 1A–3? | M80 1L–3?
  #   M81 1D–3? | M82 1G–3? | M85 1B–3? | M87 1K–3?

  r32 <- list(
    c(pos2[ 1], pos2[ 2]),           # M73: 2A vs 2B
    c(pos1[ 5], get_3rd("E")),       # M74: 1E vs 3rd
    c(pos1[ 6], pos2[ 3]),           # M75: 1F vs 2C
    c(pos1[ 3], pos2[ 6]),           # M76: 1C vs 2F
    c(pos1[ 9], get_3rd("I")),       # M77: 1I vs 3rd
    c(pos2[ 5], pos2[ 9]),           # M78: 2E vs 2I
    c(pos1[ 1], get_3rd("A")),       # M79: 1A vs 3rd
    c(pos1[12], get_3rd("L")),       # M80: 1L vs 3rd
    c(pos1[ 4], get_3rd("D")),       # M81: 1D vs 3rd
    c(pos1[ 7], get_3rd("G")),       # M82: 1G vs 3rd
    c(pos2[11], pos2[12]),           # M83: 2K vs 2L
    c(pos1[ 8], pos2[10]),           # M84: 1H vs 2J
    c(pos1[ 2], get_3rd("B")),       # M85: 1B vs 3rd
    c(pos1[10], pos2[ 8]),           # M86: 1J vs 2H
    c(pos1[11], get_3rd("K")),       # M87: 1K vs 3rd
    c(pos2[ 4], pos2[ 7])            # M88: 2D vs 2G
  )

  # Simulate one knockout round; losers get exit_rr, return winners
  play_round <- function(matches, exit_rr) {
    winners <- integer(length(matches))
    for (m in seq_along(matches)) {
      w <- sim_ko(matches[[m]][1L], matches[[m]][2L])
      l <- setdiff(matches[[m]], w)
      rr[l] <<- exit_rr      # <<- writes to simulate_tournament's rr
      winners[m] <- w
    }
    winners
  }

  # R32 → R16 → QF → SF → Final
  r32_w <- play_round(r32, 1L)

  r16   <- lapply(seq(1L, 15L, by = 2L), function(k) c(r32_w[k], r32_w[k + 1L]))
  r16_w <- play_round(r16, 2L)

  qf    <- lapply(seq(1L, 7L,  by = 2L), function(k) c(r16_w[k], r16_w[k + 1L]))
  qf_w  <- play_round(qf, 3L)

  sf    <- lapply(seq(1L, 3L,  by = 2L), function(k) c(qf_w[k],  qf_w[k  + 1L]))
  sf_w  <- play_round(sf, 4L)

  winner    <- sim_ko(sf_w[1L], sf_w[2L])
  runner_up <- setdiff(sf_w, winner)

  rr[runner_up] <- 5L
  rr[winner]    <- 6L

  rr
}

# ── 5. MONTE CARLO LOOP ───────────────────────────────────────────────────────

cat("Running", format(N_SIM, big.mark = ","), "simulations...\n")

results_mat <- matrix(0L, nrow = N, ncol = N_SIM)
for (s in seq_len(N_SIM)) {
  results_mat[, s] <- simulate_tournament()
  if (s %% 2000L == 0L) cat(" ", s, "/", N_SIM, "\n")
}

cat("Done!\n\n")

# ── 6. AGGREGATE RESULTS ──────────────────────────────────────────────────────

results <- wc_params %>%
  select(group, team) %>%
  mutate(
    p_win     = rowMeans(results_mat == 6L),
    p_final   = rowMeans(results_mat >= 5L),
    p_semi    = rowMeans(results_mat >= 4L),
    p_quarter = rowMeans(results_mat >= 3L),
    p_r16     = rowMeans(results_mat >= 2L),
    p_r32     = rowMeans(results_mat >= 1L)
  ) %>%
  arrange(desc(p_win))

cat("── Top 20 teams by win probability ──\n")
results %>%
  head(20) %>%
  select(group, team, p_win, p_final, p_semi, p_quarter) %>%
  mutate(across(starts_with("p_"),
                ~ paste0(round(.x * 100, 1), "%"))) %>%
  print(n = 20)

# ── 7. VISUALISE ─────────────────────────────────────────────────────────────

# Plot 1: Win probability bar chart — top 20 teams
p1 <- results %>%
  head(20) %>%
  mutate(team = fct_reorder(team, p_win)) %>%
  ggplot(aes(x = team, y = p_win * 100, fill = group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  labs(
    title    = "2026 World Cup Win Probabilities",
    subtitle = paste0("Monte Carlo simulation  ·  ",
                      format(N_SIM, big.mark = ","), " iterations"),
    x        = NULL,
    y        = "Win probability (%)",
    fill     = "Group",
    caption  = "Poisson GLM fitted on weighted international results (2018–2026)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p1)
ggsave("plots/win_probabilities.png", p1, width = 10, height = 7, dpi = 150)

# Plot 2: Stacked tournament path — top 24 teams by semi-final probability
# Each bar shows the marginal probability of exiting at each specific round
path_data <- results %>%
  arrange(desc(p_semi)) %>%
  head(24) %>%
  transmute(
    team, group,
    "Winner"        = p_win,
    "Runner-up"     = p_final   - p_win,
    "Semi-final"    = p_semi    - p_final,
    "Quarter-final" = p_quarter - p_semi,
    "Round of 16"   = p_r16     - p_quarter,
    "Round of 32"   = p_r32     - p_r16
  ) %>%
  pivot_longer(
    cols      = c("Round of 32","Round of 16","Quarter-final",
                  "Semi-final","Runner-up","Winner"),
    names_to  = "stage",
    values_to = "prob"
  ) %>%
  mutate(
    stage = factor(stage,
                   levels = c("Round of 32","Round of 16","Quarter-final",
                               "Semi-final","Runner-up","Winner")),
    team  = fct_reorder(team, prob, .fun = sum)
  )

p2 <- path_data %>%
  ggplot(aes(x = team, y = prob * 100, fill = stage)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  labs(
    title    = "2026 World Cup: Tournament Path Probabilities",
    subtitle = "Top 24 teams — marginal probability of exiting at each stage",
    x        = NULL,
    y        = "Probability (%)",
    fill     = "Stage reached",
    caption  = "Poisson GLM fitted on weighted international results (2018–2026)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p2)
ggsave("plots/tournament_paths.png", p2, width = 10, height = 9, dpi = 150)

# ── 8. SAVE ───────────────────────────────────────────────────────────────────

write_csv(results, "data/simulation_results.csv")
cat("\n✓ Saved: data/simulation_results.csv\n")
cat("✓ Saved: plots/win_probabilities.png\n")
cat("✓ Saved: plots/tournament_paths.png\n")
cat("\nNext step: run 04_shiny_app.R to launch the interactive dashboard.\n")
