# ============================================================
#  World Cup 2026 Predictor — Shiny Dashboard
#  File: shiny/app.R
#
#  Run from the RStudio console:
#    shiny::runApp("shiny")
#
#  Requires: data/ folder with the four CSVs produced by
#            01_data_prep.R, 02_model.R, and 03_simulation.R
# ============================================================

library(shiny)
library(tidyverse)

# ── Load pre-computed data ────────────────────────────────────────────────────
# Paths relative to the shiny/ directory where this app.R lives
team_params <- read_csv("../data/team_params.csv",        show_col_types = FALSE)
model_scl   <- read_csv("../data/model_scalars.csv",      show_col_types = FALSE)
wc_teams    <- read_csv("../data/wc_teams.csv",           show_col_types = FALSE)
sim_results <- read_csv("../data/simulation_results.csv", show_col_types = FALSE)

intercept <- model_scl$intercept
home_adv  <- model_scl$home_adv

# Master team table: group + model params + simulation probabilities
wc_full <- wc_teams %>%
  left_join(team_params %>% select(team, attack, defence), by = "team") %>%
  left_join(
    sim_results %>% select(team, p_win, p_final, p_semi, p_quarter, p_r16, p_r32),
    by = "team"
  ) %>%
  arrange(group, team)

team_choices <- sort(unique(wc_full$team))

# ── Helper functions ──────────────────────────────────────────────────────────

# Analytical match probabilities (no simulation needed)
calc_match <- function(team_a, team_b) {
  a <- wc_full %>% filter(team == team_a)
  b <- wc_full %>% filter(team == team_b)

  xg_a <- exp(intercept + a$attack + b$defence)
  xg_b <- exp(intercept + b$attack + a$defence)

  # Sum over discrete Poisson distributions up to max_g goals each
  max_g  <- 15
  p_win  <- sum(sapply(0:max_g, function(g) dpois(g, xg_a) * ppois(g - 1, xg_b)))
  p_draw <- sum(dpois(0:max_g, xg_a) * dpois(0:max_g, xg_b))
  p_loss <- 1 - p_win - p_draw

  list(xg_a = xg_a, xg_b = xg_b,
       p_win = p_win, p_draw = p_draw, p_loss = p_loss)
}

# Scoreline probability grid (used for heat map)
score_grid <- function(xg_a, xg_b, max_g = 6) {
  expand.grid(goals_a = 0:max_g, goals_b = 0:max_g) %>%
    mutate(prob = dpois(goals_a, xg_a) * dpois(goals_b, xg_b))
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(

  # Page header
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', Arial, sans-serif; }
    .page-header { background: #1a1a2e; color: white; padding: 18px 24px;
                   margin: -15px -15px 20px -15px; }
    .page-header h2 { margin: 0; font-size: 22px; }
    .page-header p  { margin: 4px 0 0; color: #aaa; font-size: 13px; }
    .summary-box { background: #f7f8fa; border-radius: 8px; padding: 16px;
                   text-align: center; margin-bottom: 8px; }
    .summary-box h3 { margin: 0 0 4px; font-size: 20px; }
    .summary-box .pct { font-size: 26px; font-weight: bold; }
    .summary-box .label { font-size: 12px; color: #888; }
  "))),

  div(class = "page-header",
    h2("🏆 2026 FIFA World Cup Predictor"),
    p("Poisson GLM · Weighted international results 2018–2026 · 50,000 Monte Carlo simulations")
  ),

  tabsetPanel(

    # ── Tab 1: Win Probabilities ───────────────────────────────────────────────
    tabPanel(
      "Win Probabilities",
      br(),
      fluidRow(
        column(3,
          sliderInput("n_teams_bar", "Show top N teams:",
                      min = 5, max = 48, value = 20, step = 1)
        ),
        column(3,
          selectInput("bar_stage", "Stage:",
            choices = c("Win tournament" = "p_win",
                        "Reach final"    = "p_final",
                        "Reach semi"     = "p_semi",
                        "Reach QF"       = "p_quarter"),
            selected = "p_win")
        )
      ),
      plotOutput("win_prob_plot", height = "580px")
    ),

    # ── Tab 2: Match Simulator ─────────────────────────────────────────────────
    tabPanel(
      "Match Simulator",
      br(),
      fluidRow(
        column(4, offset = 1,
          selectInput("team_a", "Team A:", choices = team_choices, selected = "England")
        ),
        column(4, offset = 2,
          selectInput("team_b", "Team B:", choices = team_choices, selected = "France")
        )
      ),
      br(),
      uiOutput("match_summary_ui"),
      br(),
      fluidRow(
        column(5, plotOutput("xg_bar_plot",    height = "280px")),
        column(7, plotOutput("score_heatmap",  height = "320px"))
      )
    ),

    # ── Tab 3: Team Strengths scatter ─────────────────────────────────────────
    tabPanel(
      "Team Strengths",
      br(),
      plotOutput("strength_scatter", height = "620px")
    ),

    # ── Tab 4: Groups ─────────────────────────────────────────────────────────
    tabPanel(
      "Groups",
      br(),
      fluidRow(
        column(3,
          selectInput("group_sel", "Select group:",
            choices  = setNames(LETTERS[1:12], paste("Group", LETTERS[1:12])),
            selected = "A")
        )
      ),
      tableOutput("group_table")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Tab 1: Win probability bar chart ────────────────────────────────────────

  output$win_prob_plot <- renderPlot({
    stage <- input$bar_stage
    label <- switch(stage,
      p_win     = "Win probability (%)",
      p_final   = "Probability of reaching the Final (%)",
      p_semi    = "Probability of reaching the Semi-final (%)",
      p_quarter = "Probability of reaching the Quarter-final (%)")
    title <- switch(stage,
      p_win     = "2026 World Cup Win Probabilities",
      p_final   = "Probability of Reaching the Final",
      p_semi    = "Probability of Reaching the Semi-final",
      p_quarter = "Probability of Reaching the Quarter-final")

    plot_data <- sim_results %>%
      arrange(desc(.data[[stage]])) %>%
      head(input$n_teams_bar) %>%
      mutate(value = .data[[stage]] * 100,
             team  = fct_reorder(team, value))

    ggplot(plot_data, aes(x = team, y = value, fill = group)) +
      geom_col(width = 0.72) +
      geom_text(aes(label = paste0(round(value, 1), "%")),
                hjust = -0.1, size = 3.3) +
      coord_flip(clip = "off") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(title = title,
           subtitle = "50,000 Monte Carlo simulations · Poisson GLM (2018–2026)",
           x = NULL, y = label, fill = "Group") +
      theme_minimal(base_size = 12) +
      theme(plot.title    = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(colour = "grey50", size = 10))
  })

  # ── Tab 2: Match simulator ───────────────────────────────────────────────────

  match_rv <- reactive({
    req(input$team_a, input$team_b)
    validate(need(input$team_a != input$team_b, "Select two different teams."))
    calc_match(input$team_a, input$team_b)
  })

  output$match_summary_ui <- renderUI({
    m  <- match_rv()
    ta <- input$team_a
    tb <- input$team_b

    fluidRow(
      column(4, offset = 0,
        div(class = "summary-box",
          h3(ta),
          div(class = "pct", style = "color:#2166ac;",
              paste0(round(m$p_win * 100, 1), "%")),
          div(class = "label", "to win"),
          br(),
          div(style = "font-size:18px; font-weight:bold; color:#2166ac;",
              paste0("xG: ", round(m$xg_a, 2)))
        )
      ),
      column(4,
        div(class = "summary-box",
          h3("Draw"),
          div(class = "pct", style = "color:#555;",
              paste0(round(m$p_draw * 100, 1), "%")),
          div(class = "label", "probability"),
          br(),
          div(style = "font-size: 14px; color: #888;", "Neutral ground")
        )
      ),
      column(4,
        div(class = "summary-box",
          h3(tb),
          div(class = "pct", style = "color:#d73027;",
              paste0(round(m$p_loss * 100, 1), "%")),
          div(class = "label", "to win"),
          br(),
          div(style = "font-size:18px; font-weight:bold; color:#d73027;",
              paste0("xG: ", round(m$xg_b, 2)))
        )
      )
    )
  })

  output$xg_bar_plot <- renderPlot({
    m <- match_rv()
    tibble(
      team  = c(input$team_a, input$team_b),
      xg    = c(m$xg_a, m$xg_b),
      clr   = c("#2166ac", "#d73027")
    ) %>%
      ggplot(aes(x = team, y = xg, fill = team)) +
      geom_col(width = 0.45) +
      geom_text(aes(label = round(xg, 2)), vjust = -0.4, size = 5.5, fontface = "bold") +
      scale_fill_manual(values = c("#2166ac", "#d73027")) +
      scale_y_continuous(limits = c(0, max(m$xg_a, m$xg_b) * 1.35)) +
      labs(title = "Expected Goals (neutral ground)", x = NULL, y = "xG") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
  })

  output$score_heatmap <- renderPlot({
    m <- match_rv()
    score_grid(m$xg_a, m$xg_b) %>%
      ggplot(aes(x = factor(goals_b), y = factor(goals_a), fill = prob * 100)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = paste0(round(prob * 100, 1), "%")),
                size = 3.2, colour = "grey20") +
      scale_fill_gradient(low = "#f7fbff", high = "#2166ac",
                          name = "Prob (%)") +
      labs(title = "Scoreline Probabilities",
           x = paste(input$team_b, "goals"),
           y = paste(input$team_a, "goals")) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(size = 11))
  })

  # ── Tab 3: Team strength scatter ─────────────────────────────────────────────

  output$strength_scatter <- renderPlot({
    wc_full %>%
      ggplot(aes(x = -defence, y = attack, label = team,
                 colour = p_win * 100, size = p_win * 100)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey75") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75") +
      geom_point(alpha = 0.85) +
      geom_text(size = 2.8, hjust = -0.12, vjust = 0.4, check_overlap = TRUE,
                colour = "grey20") +
      scale_colour_gradient(low = "#deebf7", high = "#e6550d",
                            name = "Win prob (%)") +
      scale_size(range = c(2, 9), guide = "none") +
      labs(
        title    = "2026 World Cup — Team Strengths",
        subtitle = "Poisson GLM parameters · point size ∝ win probability",
        x        = "Defensive strength  →  (harder to score against)",
        y        = "Attack strength  ↑  (scores more)",
        caption  = "Parameters relative to baseline | Both axes: higher = better"
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(colour = "grey50", size = 9))
  })

  # ── Tab 4: Group table ────────────────────────────────────────────────────────

  output$group_table <- renderTable({
    wc_full %>%
      filter(group == input$group_sel) %>%
      arrange(desc(p_win)) %>%
      transmute(
        Team             = team,
        `Win %`          = paste0(round(p_win     * 100, 1), "%"),
        `Final %`        = paste0(round(p_final   * 100, 1), "%"),
        `Semi-final %`   = paste0(round(p_semi    * 100, 1), "%"),
        `Quarter-final %`= paste0(round(p_quarter * 100, 1), "%"),
        `R16 %`          = paste0(round(p_r16     * 100, 1), "%"),
        `Advance (R32) %`= paste0(round(p_r32     * 100, 1), "%")
      )
  }, striped = TRUE, hover = TRUE, width = "100%")
}

# ── Launch ────────────────────────────────────────────────────────────────────

shinyApp(ui, server)
