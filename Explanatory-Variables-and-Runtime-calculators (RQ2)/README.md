# MIPLIB Project

This repository groups everything needed to experiment with **MIPLIB 2017** instances and their SCIP feature vectors.

miplib/
├── feature_extractor/ # source & build tree of the feature-extractor CLI
├── miplib_collection_set/ # full “Collection” instance set (compressed *.mps.gz)
├── miplib_benchmark_set/ # curated “Benchmark” instance set (compressed *.mps.gz)
├── instances_lists/ # plain-text lists that index subsets of the Collection
└── miplib_features_already_calculated/ # ready-to-use feature CSVs
├── *-features_original.csv # features BEFORE trivial presolve
└── *-features_after_trivial_presolving.csv # features AFTER trivial presol


## Instance sets

| Set | # instances | What it is for |
|-----|-------------|----------------|
| **Collection** (`miplib_collection_set/`) | **1 065** | Diversity-maximised grab-bag of easy, hard and still-open problems. For **training** ML models. |
| **Benchmark** (`miplib_benchmark_set/`)   | **240**   | Clean, solver-solvable subset used by the community for fair **evaluation** and competitions. |

## Instance-list files (`instances_lists/`)

This directory contains six `*.test` files—each a plain-text list (one instance name per line) that references models in the **Collection** set:

| File | Meaning |
|------|---------|
| `benchmark-v2.test` | The 240 Benchmark instances (same content as `miplib_benchmark_set/`). |
| `collection-v1.test` | All 1 065 Collection instances. |
| `easy-v17.test` | Instances that SCIP solved within minutes in the original study. |
| `hard-v29.test` | Solvable but challenging instances—hours or more. |
| `open-v28.test` | Still **unsolved** as of MIPLIB 2017. |
| `infeasible-v7.test` | Proven infeasible models. |

## Pre-computed feature matrices

`miplib_features_already_calculated/` holds two CSV files per instance:

* **`-features_original.csv`**   – features on the original model  
* **`-features_after_trivial_presolving.csv`** – features after SCIP’s *trivial* presolve (matches MIPLIB 2017 methodology)

## How the features were generated

* **Binary:** `feature_extractor/build/feature-extractor`  
  *also installed in `$PATH`, call `feature-extractor` from anywhere.*
* **Dependency:** [SCIP](https://scipopt.org/) ≥ 8.0, installed at  
  `/opt/homebrew/bin/scip` (also on `$PATH`).

## Usage of feature-extractor
The usage is as follows
```
feature-extractor -p <filename> [-s <settings>]
```

where `filename` is the path to a MIP in a format that SCIP accepts.
This will print out a comma separated feature vector
after SCIP presolved the instance.
Optionally, settings can be passed to SCIP to influence
the SCIP presolving.
In order to to deactivate presolving, a settings file
that deactivates presolving must be passed.

In order to access the header containing the feature names,
run the "feature-header-names" executable.


Example from MIPLIB dir: feature-extractor -p miplib_benchmark_set/traininstance2.mps.gz --full-problem-name
