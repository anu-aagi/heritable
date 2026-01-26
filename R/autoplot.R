#' Extract variance component values method
#' @export
extract_var_comps <- function(model){
  UseMethod("extract_var_comps")
}

#' Extract variance component values from asreml model object
#' @noRd
#' @keywords internal
#' @export
extract_var_comps.asreml <- function(model) {
  dplyr::tibble(
    component = names(model$vparameters),
    value     = model$vparameters
  ) |>
    dplyr::mutate(
      percent = value / sum(value) * 100
    )
}

#' Bar plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal
bar_plot <- function(var_comps){
  ggplot(var_comps, aes(x = component, y = value, fill = component)) +
    geom_bar(stat = "identity") +
    theme_minimal()
}

#' Stacked bar plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_minimal theme element_blank position_stack
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

#' Pie plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar geom_text theme_minimal theme element_blank position_stack
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

#' Donut plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_rect coord_polar geom_label theme_void xlim
donut_plot <- function(var_comps) {
  df <- var_comps |>
    dplyr::mutate(
      fraction = percent / 100,
      ymax = cumsum(fraction),
      ymin = dplyr::lag(ymax, default = 0)
    )

  ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = component)) +
    geom_rect() +
    coord_polar(theta = "y") +
    xlim(2, 4) +
    geom_label(
      aes(x = 3.5, y = (ymin + ymax) / 2, label = paste0(round(percent, 1), "%")),
      fill = NA,label.size = 0
    ) +
    theme_void()
}

#' Plot variance components
#' @noRd
#' @keywords internal
plot_var_comps <- function(var_comps, type = c("bar", "stacked", "pie", "donut")) {
  type <- match.arg(type)

  switch(type,
         bar = bar_plot(var_comps),
         stacked = stacked_bar_plot(var_comps),
         pie = pie_plot(var_comps),
         donut = donut_plot(var_comps)
  )
}

#' @export
#' @noRd
#' @importFrom ggplot2 autoplot
autoplot.asreml <- function(object, type = "bar", ...) {
  var_comps <- extract_var_comps.asreml(object)
  plot_var_comps(var_comps, type = type)
}


