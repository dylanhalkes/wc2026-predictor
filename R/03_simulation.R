# ============================================================
#  World Cup 2026 Predictor
#  Step 3: Monte Carlo Tournament Simulation  (v2 - with Elo)
# ============================================================

library(tidyverse)
set.seed(2026)

N_SIM <- 50000

# ── 1. LOAD & SETUP ──────────────────────────────────────────────────────────

model_scalars <- read_csv("data/model_scalars.csv", show_col_types = FALSE)
wc_params     <- read_csv("data/wc_params.csv",     show_col_types = FALSE)
wc_teams      <- read_csv("data/wc_teams.csv",      show_col_types = FALSE)

intercept <- model_scalars$intercept
home_adv  <- model_scalars$home_adv
elo_coef  <- model_scalars$elo_coef

cat("Model scalars loaded:\n")
cat("  Intercept:", round(intercept, 4), "\n")
cat("  Home adv: ", round(exp(home_adv), 3), "x\n")
cat("  Elo coef: ", round(elo_coef, 4),
    "-> per 100 Elo pts:", round((exp(elo_coef)-1)*100,1), "% more goals\n\n")

wc_full <- wc_params %>% arrange(group, team)

# Sanity checks
if (any(is.na(wc_full$attack)))
  stop("Missing attack params - re-run 02_model.R")
if (any(is.na(wc_full$current_elo))) {
  cat("WARNING: filling missing Elo with median\n")
  wc_full <- wc_full %>%
    mutate(current_elo = replace_na(current_elo, median(current_elo, na.rm=TRUE)))
}

TEAMS         <- wc_full$team
GROUPS        <- wc_full$group
N             <- nrow(wc_full)
GROUP_LETTERS <- LETTERS[1:12]

groups_idx <- setNames(
  lapply(GROUP_LETTERS, function(g) which(GROUPS == g)),
  GROUP_LETTERS
)

# Pre-compute expected goals matrix (neutral ground, Poisson + Elo)
# lambda[i,j] = exp(mu + alpha_i + beta_j + delta * (elo_i - elo_j)/100)
cat("Building lambda matrix with Elo adjustment...\n")

elo_vec      <- wc_full$current_elo
elo_diff_mat <- outer(elo_vec, elo_vec, "-") / 100   # (elo_i - elo_j)/100

lambda_mat <- exp(
  outer(wc_full$attack, wc_full$defence, "+") + intercept +
  elo_coef * elo_diff_mat
)
diag(lambda_mat) <- 0

cat("Setup complete. N =", N, "teams.\n\n")

# ── 2. THIRD-PLACE BRACKET (unchanged) ───────────────────────────────────────

eligible_3rd <- list(
  A = c("C","E","F","H","I"), B = c("E","F","G","I","J"),
  D = c("B","E","F","I","J"), E = c("A","B","C","D","F"),
  G = c("A","E","H","I","J"), I = c("C","D","F","G","H"),
  K = c("D","E","I","J","L"), L = c("E","H","I","J","K")
)

assign_3rd_places <- function(qualifying_groups) {
  slots <- names(eligible_3rd)
  assigned <- setNames(rep(NA_character_, length(slots)), slots)
  remain <- qualifying_groups
  for (iter in seq_along(slots)) {
    n_avail <- sapply(slots, function(s) {
      if (!is.na(assigned[s])) Inf
      else length(intersect(eligible_3rd[[s]], remain))
    })
    slot <- names(which.min(n_avail))
    avail <- intersect(eligible_3rd[[slot]], remain)
    chosen <- if (length(avail) > 0) sample(avail, 1) else remain[1]
    assigned[slot] <- chosen
    remain <- setdiff(remain, chosen)
  }
  assigned
}

# ── 3. CORE FUNCTIONS (unchanged) ────────────────────────────────────────────

sim_group <- function(teams4) {
  pts <- gd <- gf <- integer(4L)
  for (ii in 1L:3L) for (jj in (ii+1L):4L) {
    gi <- rpois(1L, lambda_mat[teams4[ii], teams4[jj]])
    gj <- rpois(1L, lambda_mat[teams4[jj], teams4[ii]])
    gf[ii] <- gf[ii]+gi; gf[jj] <- gf[jj]+gj
    gd[ii] <- gd[ii]+gi-gj; gd[jj] <- gd[jj]+gj-gi
    if      (gi>gj) { pts[ii] <- pts[ii]+3L }
    else if (gj>gi) { pts[jj] <- pts[jj]+3L }
    else { pts[ii] <- pts[ii]+1L; pts[jj] <- pts[jj]+1L }
  }
  ord <- order(-pts,-gd,-gf,runif(4L))
  list(idx=teams4[ord],pts=pts[ord],gd=gd[ord],gf=gf[ord])
}

sim_ko <- function(i, j) {
  gi <- rpois(1L, lambda_mat[i,j]); gj <- rpois(1L, lambda_mat[j,i])
  if      (gi>gj) return(i)
  else if (gj>gi) return(j)
  else { p <- lambda_mat[i,j]/(lambda_mat[i,j]+lambda_mat[j,i])
         if (runif(1L)<p) i else j }
}

# ── 4. TOURNAMENT SIMULATION (unchanged) ─────────────────────────────────────

simulate_tournament <- function() {
  rr <- rep(0L, N)
  pos1 <- pos2 <- pos3 <- pos4 <- integer(12L)
  t3_pts <- t3_gd <- t3_gf <- integer(12L)

  for (g in seq_len(12L)) {
    s <- sim_group(groups_idx[[g]])
    pos1[g]<-s$idx[1L]; pos2[g]<-s$idx[2L]
    pos3[g]<-s$idx[3L]; pos4[g]<-s$idx[4L]
    t3_pts[g]<-s$pts[3L]; t3_gd[g]<-s$gd[3L]; t3_gf[g]<-s$gf[3L]
  }

  t3_rank  <- order(-t3_pts,-t3_gd,-t3_gf,runif(12L))
  adv_idx  <- t3_rank[1L:8L]
  adv_grps <- GROUP_LETTERS[adv_idx]
  slot_grp <- assign_3rd_places(adv_grps)

  get_3rd <- function(slot) {
    g <- slot_grp[slot]
    if (is.na(g)) return(pos3[adv_idx[1L]])
    pos3[which(GROUP_LETTERS==g)]
  }

  r32 <- list(
    c(pos2[1],pos2[2]),          c(pos1[5],get_3rd("E")),
    c(pos1[6],pos2[3]),          c(pos1[3],pos2[6]),
    c(pos1[9],get_3rd("I")),     c(pos2[5],pos2[9]),
    c(pos1[1],get_3rd("A")),     c(pos1[12],get_3rd("L")),
    c(pos1[4],get_3rd("D")),     c(pos1[7],get_3rd("G")),
    c(pos2[11],pos2[12]),        c(pos1[8],pos2[10]),
    c(pos1[2],get_3rd("B")),     c(pos1[10],pos2[8]),
    c(pos1[11],get_3rd("K")),    c(pos2[4],pos2[7])
  )

  play_round <- function(matches, exit_rr) {
    winners <- integer(length(matches))
    for (m in seq_along(matches)) {
      w <- sim_ko(matches[[m]][1L], matches[[m]][2L])
      rr[setdiff(matches[[m]],w)] <<- exit_rr
      winners[m] <- w
    }
    winners
  }

  r32_w <- play_round(r32, 1L)
  r16   <- lapply(seq(1L,15L,2L), function(k) c(r32_w[k],r32_w[k+1L]))
  r16_w <- play_round(r16, 2L)
  qf    <- lapply(seq(1L,7L,2L),  function(k) c(r16_w[k],r16_w[k+1L]))
  qf_w  <- play_round(qf, 3L)
  sf    <- lapply(seq(1L,3L,2L),  function(k) c(qf_w[k],qf_w[k+1L]))
  sf_w  <- play_round(sf, 4L)

  winner    <- sim_ko(sf_w[1L], sf_w[2L])
  runner_up <- setdiff(sf_w, winner)
  rr[runner_up] <- 5L; rr[winner] <- 6L
  rr
}

# ── 5. MONTE CARLO LOOP ───────────────────────────────────────────────────────

cat("Running", format(N_SIM, big.mark=","), "simulations...\n")
results_mat <- matrix(0L, nrow=N, ncol=N_SIM)
for (s in seq_len(N_SIM)) {
  results_mat[,s] <- simulate_tournament()
  if (s %% 10000L == 0L) cat(" ", s, "/", N_SIM, "\n")
}
cat("Done!\n\n")

# ── 6. AGGREGATE ─────────────────────────────────────────────────────────────

results <- wc_full %>%
  select(group, team, current_elo) %>%
  mutate(
    p_win     = rowMeans(results_mat == 6L),
    p_final   = rowMeans(results_mat >= 5L),
    p_semi    = rowMeans(results_mat >= 4L),
    p_quarter = rowMeans(results_mat >= 3L),
    p_r16     = rowMeans(results_mat >= 2L),
    p_r32     = rowMeans(results_mat >= 1L)
  ) %>%
  arrange(desc(p_win))

cat("── Top 20 teams by win probability (v2: Poisson + Elo) ──\n")
results %>%
  head(20) %>%
  select(group, team, current_elo, p_win, p_final, p_semi, p_quarter) %>%
  mutate(across(starts_with("p_"), ~paste0(round(.x*100,1),"%"))) %>%
  print(n=20)

# ── 7. PLOT ───────────────────────────────────────────────────────────────────

p1 <- results %>%
  head(20) %>%
  mutate(team = fct_reorder(team, p_win)) %>%
  ggplot(aes(x=team, y=p_win*100, fill=group)) +
  geom_col(width=0.7) +
  geom_text(aes(label=paste0(round(p_win*100,1),"%")), hjust=-0.1, size=3.3) +
  coord_flip(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(title    = "2026 World Cup Win Probabilities (v2: Poisson + Elo)",
       subtitle = paste0(format(N_SIM,big.mark=","), " Monte Carlo simulations"),
       x=NULL, y="Win probability (%)", fill="Group") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold"))

print(p1)
ggsave("plots/win_probabilities_v2.png", p1, width=10, height=7, dpi=150)

# ── 8. SAVE ───────────────────────────────────────────────────────────────────

write_csv(results, "data/simulation_results.csv")
cat("\n✓ Saved: data/simulation_results.csv\n")
cat("✓ Saved: plots/win_probabilities_v2.png\n")
