# @export
extract_var_comps.asreml <- function(model) {
  tibble(
    component = names(model$vparameters),
    value     = model$vparameters
  ) |>
    mutate(
      percent = value / sum(value) * 100
    ) |>
    structure(class = c("var_comps", "tbl_df", "tbl"))
}

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

  # Pie plot
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
    # Donut plot
    donut_plot <- function(var_comps) {
        df <- var_comps |>
            mutate(
                fraction = percent / 100,
                ymax = cumsum(fraction),
                ymin = lag(ymax, default = 0)
            )

        ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = component)) +
            geom_rect() +
            coord_polar(theta = "y") +
            xlim(2, 4) +
            geom_label(
                aes(
                    x = 3.5,
                    y = (ymin + ymax) / 2,
                    label = paste0(round(percent, 1), "%")
                ),
                fill = NA,
                label.size = 0
            ) +
            theme_void()
            }

plot_var_comps <- function(var_comps, type = c("bar", "stacked", "pie", "donut")) {
    type <- match.arg(type)

    switch(type,
        bar = bar_plot(var_comps),
        stacked = stacked_bar_plot(var_comps),
        pie = pie_plot(var_comps),
        donut = donut_plot(var_comps)
    )
}

#' export
autoplot.asreml <- function(object, type = "bar", ...) {
  var_comps <- extract_var_comps.asreml(object)
  plot_var_comps(var_comps, type = type)
}

#' export
autoplot.var_comps <- function(object, type = "bar", ...) {
  plot_var_comps(object, type = type)
}
