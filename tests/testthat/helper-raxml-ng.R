raxml_test_annotations <- function() {
  data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
    clone_sequence_alignment = c("ACGT", "ACGA")
  )
}


make_raxml_test_idata <- function(
  annotations = raxml_test_annotations(),
  chains = "IGH"
) {
  annotations$sequence <- gsub(
    "-",
    "",
    annotations$clone_sequence_alignment
  )
  annotations <- annotations |>
    duckplyr::as_duckdb_tibble()

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = "sequence",
      chains = chains
    ),
    annotations = annotations
  )
}


raxml_test_additional_tips <- function() {
  data.frame(
    family_id = c("family-a", "family-b"),
    tip_id = c("family-a__germline", "family-b__germline"),
    aligned_sequence = c("ACGG", "TCGA")
  ) |>
    duckplyr::as_duckdb_tibble()
}
