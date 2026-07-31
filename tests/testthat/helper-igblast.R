make_igblast_test_immundata <- function(sequence = c(
                                         "acgtryswkmbdhvn",
                                         "ACGTRYSWKMBDHVN",
                                         NA_character_
                                       ),
                                       chain_id = c(
                                         "chain-1",
                                         "chain-2",
                                         "chain-3"
                                       )) {
  annotations <- duckplyr::as_duckdb_tibble(data.frame(
    imd_receptor_id = paste0("receptor-", seq_along(sequence)),
    imd_barcode = paste0("barcode-", seq_along(sequence)),
    imd_chain_id = chain_id,
    imd_n_chains = 1L,
    sequence = sequence
  ))

  immundata::ImmunData$new(
    schema = "sequence",
    annotations = annotations
  )
}

read_igblast_fixture <- function() {
  path <- system.file(
    "extdata",
    "igblast",
    "human-trb.tsv",
    package = "immundapter"
  )
  stopifnot(nzchar(path))

  annotations <- duckplyr::read_csv_duckdb(
    path,
    options = list(delim = "\t", header = TRUE)
  ) |>
    dplyr::mutate(
      imd_receptor_id = source_id,
      imd_barcode = source_id,
      imd_chain_id = source_id,
      imd_n_chains = 1L
    )

  immundata::ImmunData$new(
    schema = "sequence",
    annotations = annotations
  )
}
