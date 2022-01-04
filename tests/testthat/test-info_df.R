test_that("duplicated taxa in build info", {
  species <- c(" Abies_pinsapo", "Abies_alba ", "Pinus_pinaster",
             "Achyllea_milleifolium", "Rubus_plicatus", "Abies",
             "Achillum_millefolium", "    Rubus  albiflorus ", "Abies_pinsapo")

  expect_error(build.info(species, mode = "list"), "duplicated")
})

test_that("suffix in species", {
  species <- c("Abies", "Rubus", "Pinus", "Aa")
  tree <- rtree(5)
  tree$tip.label <- c("Rubus_albiflorus", "Rubus", "Pinus_sp.", "Aa_sp.", "Aa_sp2.")

  info <- build.info(species, tree)

  expect_identical(info$taxon, c("Abies_sp.", "Rubus_sp.", "Pinus_sp2.", "Aa_sp3.",
                                 "Rubus_albiflorus", "Pinus_sp.", "Aa_sp.", "Aa_sp2."))
})


test_that("genus TRUE for build.info", {
  tree <- rtree(3)
  # Test not genus in species list
  species <- c("Abies", "Rubus", "Pinus", "Aa_sp")
  tree$tip.label <- c("Abies", "Rubus", "Pinus")
  expect_error(build.info(species, tree,  genus = TRUE), ".*only genera.*")
  # Test not genus in tree tip label
  species <- c("Abies", "Rubus", "Pinus", "Aa")
  tree$tip.label <- c("Abies", "Rubus", "Pinus_sp")
  expect_error(build.info(species, tree,  genus = TRUE), ".*only genera.*")

})
