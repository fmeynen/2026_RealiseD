# config.R
# Single place for default paths and schema-version constants used across the
# simulation layer. Function argument defaults reference these (R evaluates
# defaults lazily at call time, so source order relative to this file does not
# matter), and run_all.R uses default_paths for its own path settings.
#
# Schema versions are folded into the hashes that name the generation and
# analysis output folders (see compute_data_generation_hash_from_spec() in
# data_generation_layer.R and build_analysis_run_hash() in artifact_store.R):
#   - data_generation_schema_version protects the raw per-replicate data
#     generation format (columns, semantics, RNG-to-output mapping).
#   - generation_manifest_schema_version protects the generation manifest's
#     structure.
#   - analysis_rng_scheme_version protects the per-replicate analysis RNG
#     stream scheme.
#   - results_schema_version protects the final per-replicate results schema.
#   - convergence_status_version protects the convergence_status mapping
#     rules.
#   - aggregation_schema_version protects the aggregation output schema.
# Bumping any one of them changes the corresponding run_hash, so cached
# generated data / analysis artifacts keyed on the old hash are no longer
# reused and must be regenerated.


# Default paths ------------------------------------------------------------------------------------------------------

default_paths <- list(
  generated = "data/processed/generated",
  results = "results/data",
  scripts = "scripts/simulation"
)


# Schema versions ------------------------------------------------------------------------------------------------------

# Increment this string whenever the raw per-replicate data generation format changes.
data_generation_schema_version <- "v4"

# Increment this string whenever the generation manifest's structure changes.
# Used across config.R, data_generation_layer.R and tests; not renamed to stay under 30 characters.
generation_manifest_schema_version <- "v2" # nolint: object_length_linter.

# Increment this string whenever the per-replicate analysis RNG stream scheme changes.
analysis_rng_scheme_version <- "lecuyer_analysis_substream_v2"

# Increment this string whenever the final results schema changes.
results_schema_version <- "v4"

# Increment this string whenever the convergence_status mapping rules change.
convergence_status_version <- "v2"

# Increment this string whenever the aggregation output schema changes.
aggregation_schema_version <- "v6"
