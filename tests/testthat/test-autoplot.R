test_that("Autoplot", {
  library(tidyverse)

  asreml_model_random <- readRDS(test_path("fixtures/asreml_model_random.rds"))

  # Extract the random effects
  extract_var_comps.asreml <- function(model){
    var_comps <- model$vparameters |> 
    tibble() 
    var_comps$component <- names(asreml_model_random$vparameters)
    var_comps <- var_comps |> 
    rename(value = 1) |> 
    mutate(percent = value / sum(value) * 100)
  }

  var_comps <- extract_var_comps.asreml(asreml_model_random)

  # Bar plot
  bar_plot <- function(var_comps){
  ggplot(var_comps, aes(x = component, y = value)) +
    geom_bar(stat = "identity") +
    labs(title = "Variance Components", x = "Component", y = "Value") +
    theme_minimal()
  }

  # Stacked bar plot
  stacked_bar_plot <- function(var_comps){
  ggplot(var_comps, aes(x = "", y = value, fill = component)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = paste0(round(percent, 1), "%")), 
              position = position_stack(vjust = 0.5)) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()
        ) 
  }

  # Pie chart
  pie_plot <- function(var_comps){
  ggplot(var_comps, aes(x = "", y = value, fill = component)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = paste0(round(percent, 1), "%")), 
              position = position_stack(vjust = 0.5)) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()
        )
  }

  # Donut chart
  # Compute percentages
  donut_plot <- function(var_comps){
    var_comps <- var_comps |> 
      mutate(fraction = percent/100,
           ymax = cumsum(fraction), # Compute the cumulative percentages (top of each rectangle)
           ymin = c(0, head(ymax, n=-1)))   # Compute the bottom of each rectangle
  

  ggplot(var_comps, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=component)) +
      geom_rect() +
      coord_polar(theta="y") + 
      xlim(c(2, 4)) +
      geom_label(aes(x=3.5, y=(ymin + ymax)/2, label=paste0(round(percent, 1), "%")), 
             color="black", fill=NA, label.size=0) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()
        )
        }

  autoplot    
})
