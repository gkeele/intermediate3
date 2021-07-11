context('Tmem68 mediation scan LOD scores')

test_that("equals to saved.version",{
  
  example(mediation_test)
  testthat::expect_equal(med_test$best$LR_annotation, as.vector(med_test$best$LR_mediator))
 })
