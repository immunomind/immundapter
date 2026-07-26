import argparse
import json

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from tcrdist.repertoire import TCRrep


CHAIN_SUFFIX = {
    "alpha": "a",
    "beta": "b",
    "gamma": "g",
    "delta": "d",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--mapping", required=True)
    return parser.parse_args()


def stable_cluster_labels(adjacency):
    _, labels = connected_components(
        csr_matrix(adjacency),
        directed=False,
        return_labels=True,
    )
    components = {}
    for row_index, label in enumerate(labels):
        components.setdefault(label, []).append(row_index)
    ordered = sorted(components.values(), key=min)

    stable = np.empty(len(labels), dtype=int)
    for stable_id, members in enumerate(ordered, start=1):
        stable[members] = stable_id
    return stable


def main():
    args = parse_args()
    with open(args.config, encoding="utf-8") as stream:
        config = json.load(stream)
    if isinstance(config["chains"], str):
        config["chains"] = [config["chains"]]

    cell_df = pd.read_csv(args.input, sep="\t", dtype={"input_id": str})
    receptor_cols = []
    for chain in config["chains"]:
        suffix = CHAIN_SUFFIX[chain]
        receptor_cols.extend(
            [
                f"cdr3_{suffix}_aa",
                f"v_{suffix}_gene",
                f"j_{suffix}_gene",
            ]
        )

    source_df = cell_df[["input_id", *receptor_cols]].copy()
    source_df["_source_index"] = np.arange(source_df.shape[0])

    if config["deduplicate"]:
        receptor_df = (
            source_df[receptor_cols]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        receptor_df["_deduplicated_index"] = np.arange(
            receptor_df.shape[0]
        )
        mapping = source_df.merge(
            receptor_df,
            on=receptor_cols,
            how="left",
            sort=False,
            validate="many_to_one",
        ).sort_values("_source_index")
    else:
        receptor_df = source_df[receptor_cols].reset_index(drop=True)
        receptor_df["_deduplicated_index"] = np.arange(
            receptor_df.shape[0]
        )
        mapping = source_df.copy()
        mapping["_deduplicated_index"] = np.arange(mapping.shape[0])

    analysis_df = receptor_df[receptor_cols].copy()
    analysis_df["count"] = np.bincount(
        mapping["_deduplicated_index"],
        minlength=analysis_df.shape[0],
    )

    tcrrep_args = dict(config.get("run_args", {}))
    tcrrep_args.update({
        "cell_df": analysis_df,
        "organism": config["organism"],
        "chains": config["chains"],
        "db_file": config["db_file"],
        "compute_distances": config["compute_distances"],
        "store_all_cdr": config["store_all_cdr"],
        "deduplicate": False,
        "cpus": config["cpus"],
        "index_cols": receptor_cols,
        "infer_index_cols": False,
    })
    repertoire = TCRrep(**tcrrep_args)
    if not config["compute_distances"]:
        repertoire.compute_distances()

    distances = None
    for chain in config["chains"]:
        chain_distances = getattr(repertoire, f"pw_{chain}")
        distances = (
            chain_distances.copy()
            if distances is None
            else distances + chain_distances
        )

    radius = config["radius"]
    deduplicated_adjacency = distances <= radius
    deduplicated_cluster_numbers = stable_cluster_labels(
        deduplicated_adjacency
    )
    source_to_deduplicated = mapping[
        "_deduplicated_index"
    ].to_numpy(dtype=int)
    cluster_numbers = deduplicated_cluster_numbers[
        source_to_deduplicated
    ]
    cluster_ids = np.array([
        f"cluster_{value:03d}" for value in cluster_numbers
    ])
    cluster_sizes = pd.Series(cluster_ids).value_counts()

    clusters = source_df[["input_id", *receptor_cols]].copy()
    clusters["tcrdist_cluster_id"] = cluster_ids
    clusters["tcrdist_cluster_size"] = [
        int(cluster_sizes[value]) for value in cluster_ids
    ]
    clusters.to_csv(args.clusters, sep="\t", index=False)

    mapping_output = mapping[[
        "input_id",
        "_deduplicated_index",
    ]].copy()
    mapping_output["deduplicated_receptor_id"] = [
        f"receptor_{value + 1:06d}"
        for value in mapping_output["_deduplicated_index"]
    ]
    mapping_output[[
        "input_id",
        "deduplicated_receptor_id",
    ]].to_csv(args.mapping, sep="\t", index=False)

    expanded_distances = distances[
        np.ix_(source_to_deduplicated, source_to_deduplicated)
    ]
    edge_rows = []
    for left in range(expanded_distances.shape[0]):
        for right in range(left + 1, expanded_distances.shape[1]):
            distance = int(expanded_distances[left, right])
            if distance <= radius:
                edge_rows.append(
                    {
                        "input_id_left": clusters.iloc[left]["input_id"],
                        "input_id_right": clusters.iloc[right]["input_id"],
                        "distance": distance,
                    }
                )
    pd.DataFrame(
        edge_rows,
        columns=["input_id_left", "input_id_right", "distance"],
    ).to_csv(args.edges, sep="\t", index=False)
    np.save(args.distances, expanded_distances)


if __name__ == "__main__":
    main()
