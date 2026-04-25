## Configuration

1. Copy `config_template.json` to `config.json`.
2. Set `site_name` so it exactly matches the site label in `reference/clif_hospital_geography.csv`.
3. Set `tables_path` to the local CLIF 2.1 table directory.
4. Set `file_type` to `parquet` or `csv`.

The R workflow uses `clif_adt.hospital_id` plus `reference/clif_hospital_geography.csv` to assign each encounter to the exposure county:

- Keep patient county if it is the hospital county or an adjacent county.
- Override missing or nonlocal patient county to the admitting hospital county.

`config.json` is ignored by git so local paths are not committed.
