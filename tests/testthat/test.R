context('Tmem68 mediation scan LOD scores')

test_that("equals to saved.version",{

  data("Tmem68", package = "Tmem68")
  med <- mediation_scan(target=Tmem68$target,
                        mediator=Tmem68$mediator,
                        driver=Tmem68$driver,
                        annotation=Tmem68$annotation,
                        covar=Tmem68$covar,
                        method="double-lod-diff",
                        minN = 1)

  data("Tmem68.lod", package = "Tmem68")

  # Should agree except at Tmem68
  mm <- match("Tmem68", Tmem68$annotation$symbol)
  
  m <- match(rownames(med), colnames(Tmem68$mediator))
  expect_equal(med$lod[-mm], Tmem68.lod[m][-mm])
})
