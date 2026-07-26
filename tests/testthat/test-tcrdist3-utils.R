test_that("tcrdist3 prepares paired ImmunData from distinct receptors", {
  annotations <- data.frame(
    imd_receptor_id = c(10L, 10L, 10L, 10L, 20L, 20L),
    locus = c("TRA", "TRB", "TRA", "TRB", "TRA", "TRB"),
    junction_aa = c("CAAA", "CBBB", "CAAA", "CBBB", "CCCC", "CDDD"),
    v_call = c("TRAV1", "TRBV1", "TRAV1", "TRBV1", "TRAV2", "TRBV2"),
    j_call = c("TRAJ1", "TRBJ1", "TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2")
  ) |>
    duckplyr::as_duckdb_tibble()
  idata <- immundata::ImmunData$new(
    schema = list(
      features = c("v_call", "j_call", "junction_aa"),
      chains = c("TRA", "TRB")
    ),
    annotations = annotations
  )
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
  annotations <- data.frame(
    imd_receptor_id = c(10L, 10L, 20L, 20L),
    locus = c("TRA", "TRB", "TRA", "TRB"),
    junction_aa = c("CAAA", "CBBB", "CCCC", "CDDD"),
    v_call = c("TRAV1", "TRBV1", "TRAV2", "TRBV2"),
    j_call = c("TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2")
  ) |>
    duckplyr::as_duckdb_tibble()
  idata <- immundata::ImmunData$new(
    schema = list(
      features = c("v_call", "j_call", "junction_aa"),
      chains = c("TRA", "TRB")
    ),
    annotations = annotations
  )

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

