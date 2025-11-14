test_that("compare models", {

  library(tidyverse)
  model <- list(fixed = list(), random = list())
  model$random[["asreml"]] <- asreml::asreml(yield ~ rep,
                                            random = ~ gen + rep:block,
                                            data = agridat::john.alpha,
                                            trace = FALSE)
  model$random[["lme4"]] <- lme4::lmer(yield ~ rep + (1 | gen) + (1 | rep:block), data = agridat::john.alpha)

  blup <- list("asreml" = broom.asreml::tidy(model$random[["asreml"]], "random"),
               "lme4" = broom.mixed::tidy(model$random[["lme4"]], effects = "ran_vals"))

  blue <- list("asreml" = broom.asreml::tidy(model$random[["asreml"]], "fixed") |>
                 mutate(term = str_remove(term, "_")),
               "lme4" = broom.mixed::tidy(model$random[["lme4"]], effects = "fixed"))

  vc <- list("asreml" = broom.asreml::tidy(model$random[["asreml"]], "vcomp") |>
               mutate(term = str_replace(term, "units!R", "Residual")),
               "lme4" = broom.mixed::tidy(model$random[["lme4"]], effects = "ran_pars") |>
                 mutate(estimate = estimate^2))

  plot_hist <- function(df) {
    df |>
      ggplot(aes(estimate_asreml - estimate_lme4)) +
      geom_histogram() +
      geom_vline(xintercept = 0, color = "red")
  }

  df <- full_join(blup$asreml, blup$lme4, by = c("group", "level"), suffix = c("_asreml", "_lme4"))
  plot_hist(df)
  waldo::compare(df$estimate_asreml, df$estimate_lme4, tolerance = 1e-4)

  df <- full_join(blue$asreml, blue$lme4, by = c("term"), suffix = c("_asreml", "_lme4"))
  plot_hist(df)
  waldo::compare(df$estimate_asreml, df$estimate_lme4, tolerance = 1e-4)

  df <- full_join(vc$asreml, vc$lme4, by = c("term" = "group"), suffix = c("_asreml", "_lme4"))
  plot_hist(df)
  waldo::compare(df$estimate_asreml, df$estimate_lme4, tolerance = 1e-4)


})
