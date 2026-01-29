test_that("compare models", {

  # model <- list(fixed = list(), random = list())
  # model$random[["asreml"]] <- readRDS(test_path("fixtures/asreml_model_random.rds"))
  # model$random[["lme4"]] <- readRDS(test_path("fixtures/lmer_model_random.rds"))
  #
  # blup <- list("asreml" = broom.asreml::tidy(model$random[["asreml"]], "random"),
  #              "lme4" = broom.mixed::tidy(model$random[["lme4"]], effects = "ran_vals"))
  #
  # blue <- list("asreml" = broom.asreml::tidy(model$random[["asreml"]], "fixed") |>
  #                mutate(term = str_remove(term, "_")),
  #              "lme4" = broom.mixed::tidy(model$random[["lme4"]], effects = "fixed"))
  #
  # vc <- list("asreml" = broom.asreml::tidy(model$random[["asreml"]], "vcomp") |>
  #              dplyr::mutate(term = stringr::str_replace(term, "units!R", "Residual")),
  #              "lme4" = broom.mixed::tidy(model$random[["lme4"]], effects = "ran_pars") |>
  #                dplyr::mutate(estimate = estimate^2))
  #
  # plot_hist <- function(df) {
  #   df |>
  #     ggplot2::ggplot(ggplot2::aes(estimate_asreml - estimate_lme4)) +
  #     ggplot2::geom_histogram() +
  #     ggplot2::geom_vline(xintercept = 0, color = "red")
  # }
  #
  # df <- dplyr::full_join(blup$asreml, blup$lme4, by = c("group", "level"), suffix = c("_asreml", "_lme4"))
  # plot_hist(df)
  # waldo::compare(df$estimate_asreml, df$estimate_lme4, tolerance = 1e-4)
  #
  # df <- dplyr::full_join(blue$asreml, blue$lme4, by = c("term"), suffix = c("_asreml", "_lme4"))
  # plot_hist(df)
  # waldo::compare(df$estimate_asreml, df$estimate_lme4, tolerance = 1e-4)
  #
  # df <- dplyr::full_join(vc$asreml, vc$lme4, by = c("term" = "group"), suffix = c("_asreml", "_lme4"))
  # plot_hist(df)
  # waldo::compare(df$estimate_asreml, df$estimate_lme4, tolerance = 1e-4)


})
