test_that("tcrdist3 prepares paired ImmunData from distinct receptors", {
  annotations <- tcrdist3_test_annotations()
  annotations <- rbind(annotations[1:2, ], annotations)
  idata <- make_tcrdist3_test_idata(annotations)
  config <- conf_tcrdist3(
    chains = c("alpha", "beta"),
    cdr3_cols = c(alpha = "junction_aa", beta = "junction_aa"),
    v_gene_cols = c(alpha = "v_call", beta = "v_call"),
    j_gene_cols = c(alpha = "j_call", beta = "j_call")
  )

  prepared <- tcrdist3_prepare_immundata(idata, config) |>
    dplyr::arrange(input_id) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(
    prepared,
    data.frame(
      input_id = c(10L, 20L),
      cdr3_a_aa = c("CAAA", "CCCC"),
      v_a_gene = c("TRAV1", "TRAV2"),
      j_a_gene = c("TRAJ1", "TRAJ2"),
      cdr3_b_aa = c("CBBB", "CDDD"),
      v_b_gene = c("TRBV1", "TRBV2"),
      j_b_gene = c("TRBJ1", "TRBJ2")
    )
  )
})

test_that("tcrdist3 selects one chain from paired ImmunData receptors", {
  idata <- make_tcrdist3_test_idata()

  prepared <- tcrdist3_prepare_immundata(
    idata,
    conf_tcrdist3(chains = "beta")
  ) |>
    dplyr::arrange(input_id) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(
    prepared,
    data.frame(
      input_id = c(10L, 20L),
      cdr3_b_aa = c("CBBB", "CDDD"),
      v_b_gene = c("TRBV1", "TRBV2"),
      j_b_gene = c("TRBJ1", "TRBJ2")
    )
  )
})
