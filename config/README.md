## Configuration

1. Copy `config_template.json` to `config.json`.
2. Set `site_name` so it exactly matches the site label in `reference/clif_hospital_geography.csv`.
3. Set `tables_path` to the local CLIF 2.1 table directory.
4. Set `file_type` to `parquet` or `csv`.

Example:

```json
{
  "site_name": "YOUR_SITE_NAME",
  "tables_path": "C:/path/to/local/CLIF/2.1/tables",
  "file_type": "parquet"
}
```

The R workflow uses `clif_adt.hospital_id` plus `reference/clif_hospital_geography.csv` to assign each encounter to the exposure county:

- Keep patient county if it is the hospital county or an adjacent county.
- Override missing or nonlocal patient county to the admitting hospital county.

If a site has multiple hospital IDs, each hospital should have a row in `reference/clif_hospital_geography.csv`. If a hospitalization has an unknown hospital ID, the workflow uses the first row for the configured site as the default hospital county.

The local CLIF table path can also be supplied with the `CLIF_TABLES_PATH` environment variable. `CLIF_FILE_TYPE` can override `file_type`.

The workflow currently supports CLIF tables stored as `.parquet` or `.csv`. Sites do not need parquet files if they have CSV tables. `.fst` files are not supported yet.

`config.json` is ignored by git so local paths are not committed.
