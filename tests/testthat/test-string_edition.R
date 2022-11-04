test_that("remove_spaces in species names", {
  species <- c(" Abies_pinsapo", "Abies_alba ", "Pinus_pinaster",
             "Achyllea_milleifolium", "Abies",
             "Achillum_millefolium", "    Rubus  albiflorus ")

  sp <- remove_spaces(species, "_")

  expect_true(all(!grepl("\\s", sp)))
})
