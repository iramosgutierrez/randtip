test_that("duplicated taxa in build info", {
    species <- c(" Abies_pinsapo", "Abies_alba ", "Pinus_pinaster",
                 "Achyllea_milleifolium", "Rubus_plicatus", "Abies",
                 "Achillum_millefolium", "    Rubus  albiflorus ", "Abies_pinsapo")

    expect_error(build_info(species, mode = "list"), "duplicated")
})

test_that("suffix in species", {
    species <- c("Abies", "Rubus", "Pinus", "Aa")
    tree <- rtree(5)
    tree$tip.label <- c("Rubus_albiflorus", "Rubus", "Pinus_sp.",
                        "Aa_sp.", "Aa_sp2.")

    info <- build_info(species, tree)

    expect_identical(info$taxon, c("Abies_sp.", "Rubus_sp.", "Pinus_sp2.", "Aa_sp3.",
                                 "Rubus_albiflorus", "Rubus", "Pinus_sp.",
                                 "Aa_sp.", "Aa_sp2."))
})


test_that("genus TRUE for build.info", {
    tree <- rtree(3)
    # Test not genus in species list
    species <- c("Abies", "Rubus", "Pinus", "Aa_sp")
    tree$tip.label <- c("Abies", "Rubus", "Pinus")
    expect_error(build_info(species, tree,  genus = TRUE), ".*genera only.*")
    # Test not genus in tree tip label
    species <- c("Abies", "Rubus", "Pinus", "Aa")
    tree$tip.label <- c("Abies", "Rubus", "Pinus_sp")
    expect_error(build_info(species, tree,  genus = TRUE), ".*genera only.*")

})

test_that("bad input for build.info", {
    species <- as.matrix(c("Abies", "Rubus", "Pinus", "Aa"))
    expect_error(build_info(species), "character vector")

    species <- cbind(
                    c("Abies", "Rubus", "Pinus", "Aa"),
                    c("Abies", "Rubus", "Pinus", "Aa")
                     )

    expect_error(build_info(as.data.frame(species)), "character vector")


})


test_that("The following genera were detected as ambiguous or missing", {
  catspecies <- c("Lynx_lynx",
                  "Panthera_uncia",
                  "Panthera_onca",
                  "Felis_catus",
                  "Puma_concolor",
                  "Lynx_canadensis",
                  "Panthera_tigris",
                  "Panthera_leo",
                  "Felis_silvestris")

  expect_message(build_info(catspecies, cats), "ambiguous")
})

test_that("MDCC search", {
  input <- info2input(mythology$info.list, mythology$back.tree)
  expect_equal(input[input$taxon=="Grindylowia_yorkii","MDCC" ], "Aquatia")

  info <- edit_info(mythology$info.list, "Grindylowia_yorkii", "order", NA )
  input <- info2input(info, mythology$back.tree)
  expect_equal(input[input$taxon=="Grindylowia_yorkii","MDCC" ], "Paradoxanimalia")

})


