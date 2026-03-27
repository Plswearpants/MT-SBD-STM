## Update and Deprecation Log

| Date | Type | Component | Change | Caller Migration | Legacy Recovery |
| --- | --- | --- | --- | --- | --- |
| 2026-03-26 | Deprecation | `MTSBD_all_slice.m` | Retired direct implementation; now a warning wrapper that forwards to `MTSBD_all_slice_modified.m` | Updated direct call sites in `scripts/run_test.m`, `scripts/MTSBD_block_realdata.m`, `scripts/MTSBD_block_realdata1.m`, `Dong_func/wrapper/runAllSlicesReal.m`, `examples/data/Ag111/Ag111_run.m` | Recover pre-retirement body from git history for `MTSBD_all_slice.m` |

### Append Rule

Add one row per update/deprecation. Keep entries chronological and include concrete file paths in the migration column.
