tcrdist3_test_annotations <- function() {
  data.frame(
    imd_receptor_id = rep(c(10L, 20L), each = 2L),
    locus = rep(c("TRA", "TRB"), 2L),
    junction_aa = c("CAAA", "CBBB", "CCCC", "CDDD"),
    v_call = c("TRAV1", "TRBV1", "TRAV2", "TRBV2"),
    j_call = c("TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2")
  )
}


make_tcrdist3_test_idata <- function(annotations = tcrdist3_test_annotations()) {
  annotations <- annotations |>
    duckplyr::as_duckdb_tibble()

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = c("v_call", "j_call", "junction_aa"),
      chains = c("TRA", "TRB")
    ),
    annotations = annotations
  )
}
