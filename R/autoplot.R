# Check if random effects are iid
check_ran_iid.heritable <- function(model) {
  # if not wrapped then iid?
}

#' Extract variance component values method
#' @export
extract_var_comps <- function(heritable) {
  UseMethod("extract_var_comps")
}

#' Extract variance component values from asreml model object
#' @noRd
#' @export
extract_var_comps.heritable <- function(heritable) {
  model <- attr(heritable, "model")
  target <- attr(heritable, "target")

  #TODO add check if iid before extraction

  dplyr::tibble(
    component = names(model$vparameters),
    value     = model$vparameters
  ) |>
    dplyr::mutate(
      percent = value / sum(value) * 100,
      fill_key = stringr::str_detect(component, target)
    )
}

# Make fill scale dynamically for non target variance components

make_fill_scale <- function(var_comps,
                            target_col = "#1f78b4",
                            grey_start = "grey80",
                            grey_end   = "grey40") {

  target <- var_comps |>
    dplyr::filter(fill_key) |>
    dplyr::pull(component) |>
    unique()

  non_targets <- var_comps |>
    dplyr::filter(!fill_key) |>
    dplyr::pull(component) |>
    unique()

  grey_cols <- grDevices::colorRampPalette(
    c(grey_start, grey_end)
  )(length(non_targets))

  scale_fill_manual(
    values = c(
      setNames(grey_cols, non_targets),
      setNames(target_col, target)
    ),
    breaks = c(non_targets, target),
    labels = c(non_targets, target)
  )
}

#' Bar plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal scale_fill_manual
bar_plot <- function(var_comps) {
  ggplot(var_comps, aes(x = component, y = value, fill = component)) +
    geom_bar(stat = "identity") +
    make_fill_scale(var_comps) +
    theme_minimal()
}

#' Stacked bar plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_minimal theme element_blank position_stack
stacked_bar_plot <- function(var_comps) {
  ggplot(var_comps, aes(x = "", y = value, fill = component)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = paste0(round(percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    make_fill_scale(var_comps) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

#' Pie plot
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar geom_text theme_minimal theme element_blank position_stack
pie_plot <- function(var_comps) {
  ggplot(var_comps, aes(x = "", y = value, fill = component)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = paste0(round(percent, 1), "%")),
      position = position_stack(vjust = 0.5)
    ) +
    make_fill_scale(var_comps) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
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
    make_fill_scale(var_comps) +
    xlim(2, 4) +
    geom_label(
      aes(x = 3.5, y = (ymin + ymax) / 2, label = paste0(round(percent, 1), "%")),
      fill = NA, label.size = 0
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
autoplot.heritable <- function(object, type = "bar", ...) {
  var_comps <- extract_var_comps.heritable(object)
  plot_var_comps(var_comps, type = type)
}

