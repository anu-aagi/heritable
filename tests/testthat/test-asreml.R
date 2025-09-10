test_that("heritability for asreml works", {

  fit <- asreml::asreml(yield ~ rep,
                        random = ~gen + rep:block,
                        data = agridat::john.alpha)

  expect_equal(H2(fit, target = "gen"), 0.8090841)

})
