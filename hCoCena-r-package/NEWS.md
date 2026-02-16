# hcocena 1.2.0

## Major changes

- Added an S4-first API based on `HCoCenaExperiment` and `MultiAssayExperiment`.
- Added conversion helpers: `as_hcocena()` and `as_hcobject()`.
- Added S4 accessors: `hc_mae()`, `hc_config()`, `hc_layer_results()`, `hc_integration()`, `hc_clusters()`.
- Added `hc_get_reference_data()` as transition helper for ExperimentHub-backed reference data.

## Compatibility

- Legacy global-state functions remain available and now emit migration warnings.
- Main workflow wrappers (`hc_*`) currently orchestrate existing internals to preserve behavior.

## Fixes

- Fixed typo in correlation export path (`working_director` -> `working_directory`).
- Added optional dependency guards for `propr`, `networkD3`, and `RCy3`.
- Replaced `MCDA` usage with internal criterion normalization and weighted scoring helpers.
