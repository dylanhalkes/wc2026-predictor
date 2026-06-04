# ============================================================
#  World Cup 2026 Predictor - Shiny Dashboard
#  v3: Poisson GLM + Elo + Dixon-Coles correction
#  File: shiny/app.R
#  Run: shiny::runApp("shiny")
# ============================================================

library(shiny)
library(tidyverse)

# ── Load data ─────────────────────────────────────────────────────────────────

model_scl   <- read_csv("../data/model_scalars.csv",      show_col_types = FALSE)
wc_params   <- read_csv("../data/wc_params.csv",          show_col_types = FALSE)
sim_results <- read_csv("../data/simulation_results.csv", show_col_types = FALSE)

intercept <- model_scl$intercept
home_adv  <- model_scl$home_adv
elo_coef  <- model_scl$elo_coef
rho_dc    <- model_scl$rho_dc    # Dixon-Coles correction

wc_full <- wc_params %>%
  left_join(
    sim_results %>% select(team, p_win, p_final, p_semi, p_quarter, p_r16, p_r32),
    by = "team"
  ) %>%
  arrange(group, team)

team_choices <- sort(unique(wc_full$team))

# ── Helper functions ──────────────────────────────────────────────────────────

# Dixon-Coles tau correction for low-score cells
tau_dc <- function(x, y, mu, nu, rho) {
  case_when(
    x == 0 & y == 0 ~ 1 - mu * nu * rho,
    x == 0 & y == 1 ~ 1 + mu * rho,
    x == 1 & y == 0 ~ 1 + nu * rho,
    x == 1 & y == 1 ~ 1 - rho,
    TRUE ~ 1
  )
}

# DC-corrected scoreline probability grid
score_grid_dc <- function(xg_a, xg_b, rho, max_g = 6) {
  expand.grid(goals_a = 0:max_g, goals_b = 0:max_g) %>%
    mutate(
      prob_raw = dpois(goals_a, xg_a) * dpois(goals_b, xg_b) *
                 tau_dc(goals_a, goals_b, xg_a, xg_b, rho),
      prob     = prob_raw / sum(prob_raw)   # normalise to sum to 1
    ) %>%
    select(goals_a, goals_b, prob)
}

# Match probabilities using DC-corrected score grid
calc_match <- function(team_a, team_b) {
  a <- wc_full %>% filter(team == team_a)
  b <- wc_full %>% filter(team == team_b)

  xg_a <- exp(intercept + a$attack + b$defence +
               elo_coef * (a$current_elo - b$current_elo) / 100)
  xg_b <- exp(intercept + b$attack + a$defence +
               elo_coef * (b$current_elo - a$current_elo) / 100)

  grid   <- score_grid_dc(xg_a, xg_b, rho_dc)
  p_win  <- sum(grid$prob[grid$goals_a >  grid$goals_b])
  p_draw <- sum(grid$prob[grid$goals_a == grid$goals_b])
  p_loss <- sum(grid$prob[grid$goals_a <  grid$goals_b])

  list(xg_a  = xg_a,  xg_b  = xg_b,
       elo_a = a$current_elo, elo_b = b$current_elo,
       p_win = p_win, p_draw = p_draw, p_loss = p_loss)
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(

  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', Arial, sans-serif; }
    .page-header { background: #1a1a2e; color: white; padding: 18px 24px;
                   margin: -15px -15px 20px -15px; }
    .page-header h2 { margin: 0; font-size: 22px; }
    .page-header p  { margin: 4px 0 0; color: #aaa; font-size: 13px; }
    .summary-box { background: #f7f8fa; border-radius: 8px; padding: 16px;
                   text-align: center; margin-bottom: 8px; }
    .summary-box h3  { margin: 0 0 4px; font-size: 20px; }
    .summary-box .pct   { font-size: 26px; font-weight: bold; }
    .summary-box .label { font-size: 12px; color: #888; }
    .elo-badge { display: inline-block; background: #e8f4fd; color: #2166ac;
                 border-radius: 4px; padding: 2px 8px; font-size: 12px;
                 font-weight: bold; margin-top: 4px; }
  "))),

  div(class = "page-header",
    h2("2026 FIFA World Cup Predictor"),
    p("Poisson GLM + Elo ratings + Dixon-Coles correction  |  50,000 Monte Carlo simulations")
  ),

  tabsetPanel(

    tabPanel("Win Probabilities",
      br(),
      fluidRow(
        column(3, sliderInput("n_teams_bar", "Show top N teams:",
                              min = 5, max = 48, value = 20, step = 1)),
        column(3, selectInput("bar_stage", "Stage:",
          choices = c("Win tournament" = "p_win",
                      "Reach final"    = "p_final",
                      "Reach semi"     = "p_semi",
                      "Reach QF"       = "p_quarter"),
          selected = "p_win"))
      ),
      plotOutput("win_prob_plot", height = "580px")
    ),

    tabPanel("Match Simulator",
      br(),
      fluidRow(
        column(4, offset = 1,
          selectInput("team_a", "Team A:", choices = team_choices, selected = "England")),
        column(4, offset = 2,
          selectInput("team_b", "Team B:", choices = team_choices, selected = "France"))
      ),
      br(),
      uiOutput("match_summary_ui"),
      br(),
      fluidRow(
        column(5, plotOutput("xg_bar_plot",   height = "280px")),
        column(7, plotOutput("score_heatmap", height = "320px"))
      )
    ),

    tabPanel("Team Strengths",
      br(),
      plotOutput("strength_scatter", height = "620px")
    ),

    tabPanel("Groups",
      br(),
      fluidRow(
        column(3, selectInput("group_sel", "Select group:",
          choices  = setNames(LETTERS[1:12], paste("Group", LETTERS[1:12])),
          selected = "A"))
      ),
      tableOutput("group_table")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Tab 1: Win probability bar chart
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
      labs(title    = title,
           subtitle = "50,000 simulations  |  Poisson GLM + Elo + Dixon-Coles",
           x = NULL, y = label, fill = "Group") +
      theme_minimal(base_size = 12) +
      theme(plot.title    = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(colour = "grey50", size = 10))
  })

  # Tab 2: Match simulator
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
      column(4,
        div(class = "summary-box",
          h3(ta),
          div(class = "elo-badge", paste0("Elo: ", round(m$elo_a))),
          br(), br(),
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
          div(class = "elo-badge",
              paste0("Gap: ", round(m$elo_a) - round(m$elo_b), " pts")),
          br(), br(),
          div(class = "pct", style = "color:#555;",
              paste0(round(m$p_draw * 100, 1), "%")),
          div(class = "label", "probability"),
          br(),
          div(style = "font-size:14px; color:#888;", "Neutral ground")
        )
      ),
      column(4,
        div(class = "summary-box",
          h3(tb),
          div(class = "elo-badge", paste0("Elo: ", round(m$elo_b))),
          br(), br(),
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
    tibble(team = c(input$team_a, input$team_b),
           xg   = c(m$xg_a, m$xg_b)) %>%
      ggplot(aes(x = team, y = xg, fill = team)) +
      geom_col(width = 0.45) +
      geom_text(aes(label = round(xg, 2)), vjust = -0.4, size = 5.5, fontface = "bold") +
      scale_fill_manual(values = setNames(c("#2166ac", "#d73027"),
                                          c(input$team_a, input$team_b))) +
      scale_y_continuous(limits = c(0, max(m$xg_a, m$xg_b) * 1.35)) +
      labs(title = "Expected Goals (neutral ground)", x = NULL, y = "xG") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
  })

  output$score_heatmap <- renderPlot({
    m <- match_rv()
    # Uses Dixon-Coles corrected probabilities
    score_grid_dc(m$xg_a, m$xg_b, rho_dc) %>%
      ggplot(aes(x = factor(goals_b), y = factor(goals_a), fill = prob * 100)) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = paste0(round(prob * 100, 1), "%")),
                size = 3.2, colour = "grey20") +
      scale_fill_gradient(low = "#f7fbff", high = "#2166ac", name = "Prob (%)") +
      labs(title   = "Scoreline Probabilities (Dixon-Coles corrected)",
           x = paste(input$team_b, "goals"),
           y = paste(input$team_a, "goals")) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
  })

  # Tab 3: Team strengths scatter
  output$strength_scatter <- renderPlot({
    wc_full %>%
      ggplot(aes(x = -defence, y = attack, label = team,
                 colour = p_win * 100, size = p_win * 100)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey75") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75") +
      geom_point(alpha = 0.85) +
      geom_text(size = 2.8, hjust = -0.12, vjust = 0.4,
                check_overlap = TRUE, colour = "grey20") +
      scale_colour_gradient(low = "#deebf7", high = "#e6550d",
                            name = "Win prob (%)") +
      scale_size(range = c(2, 9), guide = "none") +
      labs(title    = "2026 World Cup - Team Strengths",
           subtitle = "Poisson GLM + Elo  |  point size proportional to win probability",
           x = "Defensive strength (harder to score against)",
           y = "Attack strength (scores more)") +
      theme_minimal(base_size = 11) +
      theme(plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(colour = "grey50", size = 9))
  })

  # Tab 4: Group table
  output$group_table <- renderTable({
    wc_full %>%
      filter(group == input$group_sel) %>%
      arrange(desc(p_win)) %>%
      transmute(
        Team              = team,
        `Elo`             = round(current_elo),
        `Win %`           = paste0(round(p_win     * 100, 1), "%"),
        `Final %`         = paste0(round(p_final   * 100, 1), "%"),
        `Semi-final %`    = paste0(round(p_semi    * 100, 1), "%"),
        `Quarter-final %` = paste0(round(p_quarter * 100, 1), "%"),
        `R16 %`           = paste0(round(p_r16     * 100, 1), "%"),
        `Advance (R32) %` = paste0(round(p_r32     * 100, 1), "%")
      )
  }, striped = TRUE, hover = TRUE, width = "100%")
}

shinyApp(ui, server)
