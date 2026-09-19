library(shiny)
library(dplyr)
library(ggplot2)
library(grid)
library(survival)

# this code will plot TSSPC, KP curve, LLP and CI
# For a categorical var with multiple category, it will plot TSSPC for each category, LLP and CI


make_survival_probability_curve <- function(dat, tau = NULL){
  
  if(is.null(tau)){
    tau <- max(dat$time[dat$arm == "Control" & dat$status == 1], na.rm = TRUE)
  }
  
  event_times <- sort(unique(dat$time[dat$status == 1 & dat$time <= tau]))
  
  S0 <- 1
  S1 <- 1
  
  curve_points <- data.frame(x = 1, y = 1, pattern = "start")
  step_table <- data.frame()
  
  for(t in event_times){
    
    Y0 <- sum(dat$arm == "Control"   & dat$time >= t)
    Y1 <- sum(dat$arm == "Treatment" & dat$time >= t)
    
    d0 <- sum(dat$arm == "Control"   & dat$time == t & dat$status == 1)
    d1 <- sum(dat$arm == "Treatment" & dat$time == t & dat$status == 1)
    
    S0_before <- S0
    S1_before <- S1
    
    S0_after <- ifelse(Y0 > 0, S0_before * (1 - d0 / Y0), S0_before)
    S1_after <- ifelse(Y1 > 0, S1_before * (1 - d1 / Y1), S1_before)
    
    llp_term <- ifelse(
      d0 > 0 & Y0 > 0,
      S1_after * S0_before * (d0 / Y0),
      0
    )
    
    step_table <- rbind(
      step_table,
      data.frame(
        time = t,
        Y0 = Y0, d0 = d0,
        Y1 = Y1, d1 = d1,
        S0_before = S0_before,
        S0_after = S0_after,
        S1_before = S1_before,
        S1_after = S1_after,
        llp_term = llp_term
      )
    )
    
    if(d0 > 0 & d1 == 0){
      curve_points <- rbind(
        curve_points,
        data.frame(x = S0_after, y = S1_before, pattern = "horizontal")
      )
    }
    
    if(d0 == 0 & d1 > 0){
      curve_points <- rbind(
        curve_points,
        data.frame(x = S0_before, y = S1_after, pattern = "vertical")
      )
    }
    
    if(d0 > 0 & d1 > 0){
      curve_points <- rbind(
        curve_points,
        data.frame(x = S0_before, y = S1_after, pattern = "vertical polygon"),
        data.frame(x = S0_after,  y = S1_after, pattern = "horizontal polygon")
      )
    }
    
    S0 <- S0_after
    S1 <- S1_after
  }
  
  curve_points <- rbind(
    curve_points,
    data.frame(x = 0, y = tail(curve_points$y, 1), pattern = "final horizontal")
  )
  
  LLP <- sum(step_table$llp_term, na.rm = TRUE)
  
  step_table$Q1 <- rev(cumsum(rev(step_table$llp_term)))
  step_table$Q0 <- step_table$S1_after * step_table$S0_after - step_table$Q1
  
  var_treat <- sum(
    ifelse(
      step_table$d1 > 0 & step_table$Y1 > step_table$d1,
      step_table$Q1^2 * step_table$d1 /
        (step_table$Y1 * (step_table$Y1 - step_table$d1)),
      0
    ),
    na.rm = TRUE
  )
  
  var_control <- sum(
    ifelse(
      step_table$d0 > 0 & step_table$Y0 > step_table$d0,
      step_table$Q0^2 * step_table$d0 /
        (step_table$Y0 * (step_table$Y0 - step_table$d0)),
      0
    ),
    na.rm = TRUE
  )
  
  var_LLP <- var_treat + var_control
  se_LLP <- sqrt(var_LLP)
  
  eps <- 1e-8
  LLP_safe <- min(max(LLP, eps), 1 - eps)
  
  logit_LLP <- log(LLP_safe / (1 - LLP_safe))
  se_logit <- se_LLP / (LLP_safe * (1 - LLP_safe))
  
  lower_logit <- logit_LLP - 1.96 * se_logit
  upper_logit <- logit_LLP + 1.96 * se_logit
  
  lower <- exp(lower_logit) / (1 + exp(lower_logit))
  upper <- exp(upper_logit) / (1 + exp(upper_logit))
  
  list(
    tau = tau,
    curve_points = curve_points,
    step_table = step_table,
    LLP = LLP,
    variance = var_LLP,
    SE = se_LLP,
    CI_lower = lower,
    CI_upper = upper
  )
}

ui <- fluidPage(
  
  titlePanel("Two-Sample Survival Probability Curve"),
  
  sidebarLayout(
    
    sidebarPanel(
      fileInput("file", "Upload CSV file", accept = ".csv"),
      
      uiOutput("column_selectors"),
      uiOutput("arm_selectors"),
      
      selectInput(
        "status_coding",
        "Status/censor coding",
        choices = c(
          "1 = censored, 0 = event" = "censor1",
          "1 = event, 0 = censored" = "event1"
        ),
        selected = "censor1"
      ),
      
      checkboxInput(
        "stratify_curve",
        "Plot by categorical variable",
        value = FALSE
      ),
      
      conditionalPanel(
        condition = "input.stratify_curve == true",
        uiOutput("strata_selector"),
        uiOutput("category_selector")
      ),
      
      checkboxInput(
        "show_km",
        "Show Kaplan-Meier curve",
        value = FALSE
      )
    ),
    
    mainPanel(
      plotOutput("llp_plot", height = "700px"),
      
      conditionalPanel(
        condition = "input.show_km == true",
        plotOutput("km_plot", height = "500px")
      )
    )
  )
)

server <- function(input, output, session){
  
  raw_data <- reactive({
    req(input$file)
    read.csv(input$file$datapath, stringsAsFactors = FALSE)
  })
  
  output$column_selectors <- renderUI({
    req(raw_data())
    
    df <- raw_data()
    
    tagList(
      selectInput("arm_col", "Treatment arm column", choices = names(df)),
      selectInput("time_col", "Time column", choices = names(df)),
      selectInput("status_col", "Status/censor column", choices = names(df))
    )
  })
  
  output$arm_selectors <- renderUI({
    req(raw_data(), input$arm_col)
    
    arms <- sort(unique(as.character(raw_data()[[input$arm_col]])))
    arms <- arms[!is.na(arms)]
    
    tagList(
      selectInput("control_arm", "Control arm label", choices = arms),
      selectInput("treatment_arm", "Treatment arm label", choices = arms)
    )
  })
  
  output$strata_selector <- renderUI({
    req(raw_data())
    
    selectInput(
      "strata_col",
      "Categorical grouping variable",
      choices = names(raw_data())
    )
  })
  
  output$category_selector <- renderUI({
    req(raw_data(), input$strata_col)
    
    cats <- sort(unique(as.character(raw_data()[[input$strata_col]])))
    cats <- cats[!is.na(cats) & cats != ""]
    
    checkboxGroupInput(
      "selected_categories",
      "Select categories to plot",
      choices = cats,
      selected = cats
    )
  })
  
  surv_data <- reactive({
    req(
      raw_data(),
      input$arm_col,
      input$time_col,
      input$status_col,
      input$control_arm,
      input$treatment_arm
    )
    
    df <- raw_data()
    
    tmp <- df %>%
      filter(.data[[input$arm_col]] %in% c(input$control_arm, input$treatment_arm)) %>%
      mutate(
        arm = case_when(
          .data[[input$arm_col]] == input$control_arm ~ "Control",
          .data[[input$arm_col]] == input$treatment_arm ~ "Treatment"
        ),
        time = as.numeric(.data[[input$time_col]]),
        raw_status = as.numeric(.data[[input$status_col]])
      )
    
    if(input$status_coding == "censor1"){
      tmp <- tmp %>% mutate(status = 1 - raw_status)
    } else {
      tmp <- tmp %>% mutate(status = raw_status)
    }
    
    if(isTRUE(input$stratify_curve)){
      req(input$strata_col, input$selected_categories)
      
      tmp <- tmp %>%
        mutate(strata = as.character(.data[[input$strata_col]])) %>%
        filter(
          !is.na(strata),
          strata != "",
          strata %in% input$selected_categories
        )
    } else {
      tmp <- tmp %>%
        mutate(strata = "Overall")
    }
    
    tmp %>%
      filter(!is.na(time), !is.na(status), !is.na(arm)) %>%
      select(arm, time, status, strata)
  })
  
  res_by_strata <- reactive({
    
    dat <- surv_data()
    strata_levels <- unique(dat$strata)
    
    all_curves <- data.frame()
    all_results <- data.frame()
    
    for(g in strata_levels){
      
      dat_g <- dat %>% filter(strata == g)
      
      valid_group <-
        sum(dat_g$arm == "Control") > 0 &&
        sum(dat_g$arm == "Treatment") > 0 &&
        sum(dat_g$arm == "Control" & dat_g$status == 1) > 0
      
      if(valid_group){
        
        res_g <- make_survival_probability_curve(dat_g)
        
        curve_g <- res_g$curve_points %>%
          mutate(strata = g)
        
        all_curves <- rbind(all_curves, curve_g)
        
        all_results <- rbind(
          all_results,
          data.frame(
            category = g,
            LLP = res_g$LLP,
            SE = res_g$SE,
            CI_lower = res_g$CI_lower,
            CI_upper = res_g$CI_upper,
            tau = res_g$tau
          )
        )
      }
    }
    
    list(
      curve_points = all_curves,
      results = all_results
    )
  })
  
  km_data <- reactive({
    dat <- surv_data()
    
    if(isTRUE(input$stratify_curve)){
      fit <- survfit(Surv(time, status) ~ arm + strata, data = dat)
    } else {
      fit <- survfit(Surv(time, status) ~ arm, data = dat)
    }
    
    s <- summary(fit)
    
    data.frame(
      time = s$time,
      surv = s$surv,
      group = s$strata
    )
  })
  
  output$llp_plot <- renderPlot({
    
    curve_dat <- res_by_strata()$curve_points
    result_dat <- res_by_strata()$results
    
    if(nrow(curve_dat) == 0){
      plot.new()
      text(
        0.5, 0.5,
        "No valid curve could be computed.\nEach selected category must have both arms and at least one control event."
      )
      return()
    }
    
    if(isTRUE(input$stratify_curve)){
      p <- ggplot(curve_dat, aes(x = x, y = y, color = strata, linetype = strata))
      legend_title <- input$strata_col
      plot_title <- "Two-Sample Survival Probability Curve by Category"
      
      label_text <- paste0(
        result_dat$category,
        ": LLP = ", round(result_dat$LLP, 3),
        " (95% CI: ",
        round(result_dat$CI_lower, 3), ", ",
        round(result_dat$CI_upper, 3), ")",
        collapse = "\n"
      )
      
    } else {
      p <- ggplot(curve_dat, aes(x = x, y = y))
      legend_title <- NULL
      plot_title <- "Two-Sample Survival Probability Curve"
      
      label_text <- paste0(
        "LLP = ", round(result_dat$LLP, 3),
        "\n95% CI: (",
        round(result_dat$CI_lower, 3), ", ",
        round(result_dat$CI_upper, 3), ")"
      )
    }
    
    p +
      annotate("rect",
               xmin = 0, xmax = 1,
               ymin = 0, ymax = 1,
               fill = "grey90") +
      
      annotate("polygon",
               x = c(0, 1, 1),
               y = c(0, 0, 1),
               fill = "grey75") +
      
      geom_abline(intercept = 0, slope = 1,
                  colour = "grey55", linewidth = 0.8) +
      
      geom_path(linewidth = 0.9) +
      
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      
      labs(
        title = plot_title,
        x = paste0("Survival for Control Arm: ", input$control_arm),
        y = paste0("Survival for Treatment Arm: ", input$treatment_arm),
        color = legend_title,
        linetype = legend_title
      ) +
      
      annotate("text",
               x = 0.18, y = 0.82,
               label = "Increasing benefit\nfrom treatment",
               size = 5) +
      
      annotate("segment",
               x = 0.14, y = 0.88,
               xend = 0.01, yend = 1.00,
               linetype = "dashed",
               linewidth = 0.7,
               arrow = arrow(type = "closed",
                             length = unit(0.10, "inches"))) +
      
      annotate("text",
               x = 0.78, y = 0.18,
               label = "Increasing benefit\nfrom control",
               size = 5) +
      
      annotate("segment",
               x = 0.86, y = 0.14,
               xend = 1.00, yend = 0.02,
               linetype = "dashed",
               linewidth = 0.7,
               arrow = arrow(type = "closed",
                             length = unit(0.10, "inches"))) +
      
      annotate("text",
               x = 0.47, y = 0.45,
               label = "No treatment benefit",
               angle = 36.8,
               size = 5) +
      
      annotate("text",
               x = 0.30, y = 0.10,
               label = label_text,
               size = 4,
               fontface = "bold") +
      
      theme_classic(base_size = 16) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_blank()
      )
  })
  
  
  output$km_plot <- renderPlot({
    
    req(input$show_km)
    
    ggplot(km_data(), aes(x = time, y = surv, color = group, linetype = group)) +
      geom_step(linewidth = 1.1) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title = "Kaplan-Meier Survival Curves",
        x = "Time",
        y = "Survival probability",
        color = "Group",
        linetype = "Group"
      ) +
      theme_classic(base_size = 16) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        legend.position = "right"
      )
  })
}

shinyApp(ui, server)