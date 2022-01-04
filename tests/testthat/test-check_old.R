source("../../testing/old_code/dataframes.R")
source("../../testing/old_code/edit.info.R")
source("../../testing/old_code/custom.branch.R")
source("../../testing/old_code/utils.R")
source("../../testing/old_code/rand.tip.R")


tree25 <- read.tree("../../data-raw/25tree.tre")
species25 <- c("Abies_pinsapo", "Abies_alba ", "Pinus_pinaster",
            "Achyllea_milleifolium", "Rubus_plicatus", "Abies",
            "Achillum_millefolium")

myth.sp.list <- mythology$sp.list
myth.back.info <- mythology$info.backbone
myth.back.tree <- mythology$back.tree

#TODO mode backbone
test_that("equal to build.info_old mode list", {

  m <- "list"
  suppressWarnings({
    info_new <<- build.info(species25, tree25, mode = m)
    info_old <<- build.info_old(species25, tree25, mode = m)
  })
  expect_identical(
    info_new,
    info_old
  )
})

test_that("equal to build.info old mode backbone", {
  myth.info.noranks_new <<- build.info(myth.sp.list, tree = myth.back.tree,
                                       find.ranks = FALSE, mode = "backbone")
  myth.info.noranks_old <<- build.info_old(myth.sp.list, tree = myth.back.tree,
                                           find.ranks = FALSE, mode = "backbone")
  expect_identical(
    myth.info.noranks_new,
    myth.info.noranks_old
  )
})

test_that("equal check.info to old", {
  expect_identical(
    check.info(info_new, tree25),
    check.info_old(info_old, tree25)
  )
  expect_identical(
    check.info(myth.back.info, myth.back.tree),
    check.info_old(myth.back.info, myth.back.tree)
  )
})


test_that("equal to info2input old", {

  myth.input_new <<- info2input(myth.back.info, myth.back.tree)
  myth.input_old <<- info2input_old(myth.back.info, myth.back.tree)

  expect_identical(
    myth.input_new,
    myth.input_old
  )
})

test_that("equal to rand.tip old", {

  set.seed(1)
  new.tree_new <- rand.tip(myth.input_new, myth.back.tree,
                           forceultrametric = TRUE,
                           prune = FALSE)
  set.seed(1)
  new.tree_old <- rand.tip_old(myth.input_old, myth.back.tree,
                               forceultrametric = TRUE,
                               prune = FALSE)

  expect_identical(
    new.tree_new,
    new.tree_old
  )
})
