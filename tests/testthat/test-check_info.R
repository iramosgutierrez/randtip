
test_that("There may be misspelling errors", {

  expect_message(check_info(mythology$info.list, mythology$back.tree), "TYPO")
  expect_message(check_info(mythology$info.list, mythology$back.tree), "same taxon")
  expect_message(check_info(mythology$info.list, mythology$back.tree), "not ultrametric")
  
})

test_that("Phyleticity of groups", {
  
  checkinfo <- check_info(mythology$info.list, mythology$back.tree)
  expect_equal(as.character(checkinfo[checkinfo$taxon=="Gorgona_medusii","Typo" ]), "TRUE")
  expect_equal(checkinfo[checkinfo$taxon=="Gorgona_medusii","Typo.names" ], "Gorgona_medusi")

    
})

test_that("Edit info", {
  
  info2 <- edit_info(mythology$info.list, "Gorgona_medusii", "taxon", "Gorgona_medusi" )
  checkinfo <- check_info(info2, mythology$back.tree)
  expect_equal(as.character(checkinfo[checkinfo$taxon=="Gorgona_medusi","Typo" ]), "FALSE")
})

test_that("Edit tree", {
    tree2 <- edit_tree(mythology$back.tree, "Gorgona_medusi", "Gorgona_medusii" )
  checkinfo <- check_info(mythology$info.list, tree2)
  expect_equal(as.character(checkinfo[checkinfo$taxon=="Gorgona_medusii","Typo" ]), "FALSE")
})
