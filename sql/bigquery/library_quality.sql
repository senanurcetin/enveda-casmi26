-- NEXUS CASMI26 | BigQuery training-data QA and library profile.
--
-- Contract: a BigQuery table loaded from the official train.parquet with
-- original column names, including repeated FLOAT64 ms2_mzs and
-- ms2_normalized_intensities. Replace identifiers with render_bigquery_sql.py.
--
-- This query is read-only. Review BigQuery dry-run bytes/cost before execution.
WITH quality_by_source AS (
  SELECT
    COALESCE(NULLIF(ingest_lib, ''), 'UNKNOWN') AS source_library,
    COUNT(*) AS spectrum_count,
    COUNT(DISTINCT inchikey14) AS distinct_2d_structures,
    COUNTIF(inchikey14 IS NULL OR LENGTH(inchikey14) != 14) AS missing_structure_keys,
    COUNTIF(precursor_mz IS NULL OR precursor_mz <= 0) AS invalid_precursor_mz,
    COUNTIF(
      ms2_mzs IS NULL OR ms2_normalized_intensities IS NULL
      OR ARRAY_LENGTH(ms2_mzs) != ARRAY_LENGTH(ms2_normalized_intensities)
    ) AS mismatched_peak_arrays,
    COUNTIF(
      ms2_mzs IS NOT NULL AND ARRAY_LENGTH(ms2_mzs) = 0
    ) AS empty_spectra,
    COUNTIF(
      precursor_error_ppm IS NOT NULL AND ABS(precursor_error_ppm) > 10
    ) AS precursor_error_over_10ppm,
    COUNTIF(adduct IN (
      '[M+H]+', '[M+NH4]+', '[M-H2O+H]+', '[M-2H2O+H]+',
      '[M+Na]+', '[M+K]+', '[M-H]-', '[M-H2O-H]-',
      '[M+CH2O2-H]-', '[M+Cl]-'
    )) AS spectra_with_test_adduct,
    COUNTIF(ionization_mode = 'positive') AS positive_ion_mode_spectra
  FROM `{{PROJECT_ID}}.{{DATASET_ID}}.{{TABLE_NAME}}`
  GROUP BY 1
)
SELECT
  source_library,
  spectrum_count,
  distinct_2d_structures,
  missing_structure_keys,
  invalid_precursor_mz,
  mismatched_peak_arrays,
  empty_spectra,
  precursor_error_over_10ppm,
  ROUND(SAFE_DIVIDE(spectra_with_test_adduct, spectrum_count), 4)
    AS test_adduct_share,
  ROUND(SAFE_DIVIDE(positive_ion_mode_spectra, spectrum_count), 4)
    AS positive_ion_share
FROM quality_by_source
ORDER BY spectrum_count DESC, source_library;
