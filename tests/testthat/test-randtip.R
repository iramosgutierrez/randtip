

test_that("Randtip functioning", {
  
input <- info2input(mythology$info.list, mythology$back.tree)
expect_message(rand.tip(input, mythology$back.tree)  , "not ultrametric")

seed <- sample(1:999999, 1)

set.seed(seed)
tree1 <- rand.tip(input, mythology$back.tree, forceultrametric = T)

set.seed(seed)
tree2 <- rand.tip(input, mythology$back.tree, forceultrametric = T)


expect_equal(tree1, tree2)
})



test_that("Randtip functioning", {
  
  input <- info2input(mythology$info.list, mythology$back.tree)

  seed <- sample(1:999999, 1)
  
  set.seed(seed)
  tree1 <- rand.tip(input, mythology$back.tree, forceultrametric = T, rand.type = "polytomy")
  
  set.seed(seed)
  tree2 <- rand.tip(input, mythology$back.tree, forceultrametric = T, rand.type = "polytomy")
  
  
  expect_equal(tree1, tree2)
})
