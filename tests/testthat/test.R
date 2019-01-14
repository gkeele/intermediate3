context('Tmem68 mediation scan LOD scores')

test_that("equals to saved.version",{

  data("Tmem68")
  med <- mediation_scan(target=Tmem68$target,
                        mediator=Tmem68$mediator,
                        driver=Tmem68$qtl.geno,
                        annotation=Tmem68$annotation,
                        covar=Tmem68$covar,
                        method="double-lod-diff",
                        minN = 1)

  data("Tmem68.lod")

  m <- match(rownames(med), colnames(Tmem68$mediator))
  expect_equal(med$lod, Tmem68.lod[m])
})
