# Data card — spectra, structures and derived evidence

## Sources and access

| Source | Recorded scope | Access and evidence |
|---|---|---|
| Competition training | 2,539,608 spectra, 275,810 standard structure keys, 11 libraries | Accept applicable Kaggle terms; raw spectra are not redistributed. [Historical profile](reports/data_profile.json) |
| Public examples | 1,213 spectra / 400 molecules | All have exact stored training-content matches. [Executed full-content audit](reports/public_overlap_audit.json) |
| COCONUT / PubChemLite | Task-defined candidate knowledge; study-specific byte hashes | Retain release identifiers, source licensing and attribution. [Fresh internal input manifest](reports/advanced_study.json) |
| MassBank 2026.03 CC0 subset | 491 eligible molecules; zero novel against original training | Explicit per-record CC0 only. [Source preparation](reports/external_source_preparation.json) |
| Direct-contribution GNPS subsets | 85 molecules initially; 936 in the expanded source; seven novel in each against original training | Direct-contribution CC0 policy; imported/inferred libraries excluded. [Source methods and pins](docs/studies/external-source-validation.md) |

The GNPS source stages overlap; these are not independent replications. Competition access does not grant redistribution rights to all source records, candidate databases, third-party weights or branding. Public-source inputs are prepared separately and bound to hashes; no new laboratory acquisition is claimed.

## Unit of analysis and information boundaries

The prediction/evaluation unit is a molecule. Multiple spectra, libraries and tautomer-equivalent representations require grouping before splitting and reference exclusion. Summing per-library distinct counts double-counts shared identities. Standard keys and metric InChIKey14 are different identity definitions; neither alone establishes unseen chemical-family generalization or stereochemical identification.

The public-example audit confirms precursor, ordered m/z and normalized intensity identity for every example, consistent with adduct and ionization mode. Numerical content matching does not authenticate original acquisition provenance. Public examples support interface checks, not an independent accuracy claim.

| Executed population | Separation and retention policy | Evidence |
|---|---|---|
| Frozen retrospective V-in | All 250 queries retained; 1,710 mass-conflicting identities and 13 unparseable candidates quarantined; 1,071,863 candidate structures fingerprint-verified | [Rebuild and result report](reports/revalidation_review.md) |
| Fresh internal benchmark | 600/150/250 train/development/test molecules, one per generic family; historical V-in families excluded; selected metric/family identities do not cross supervised splits | [Population and limitations](docs/studies/fresh-training-family-holdout.md) |
| Separate purged-reference experiment | 150 selected external molecules; matching families/standard/metric identities blocked in the full training/reference corpus; 261,826 rows removed, 2,277,782 retained; new internal 300/100/150 partitions | [Purge, source amendment and complete report](reports/purged_external_review.md) |

In the purged experiment, one undefined-family training label blocks its entire declared identity under the registered outcome-blind amendment. Historical internal partition exclusions are reconstructed before filtering. Zero selected-family, standard/metric identity and exact-content overlap is checked before fresh fitting. The explicit `CASMI_DATA` snapshot takes precedence over automatic Kaggle discovery.

## Distribution and remaining risks

V-in is a deliberately common natural-product set; library, instrument, ionization mode, adduct, peak count and precursor-error distributions differ. In the fresh internal test, 141/250 queries are pluskal_ms2 and 56/250 GNPS; chemical-family balance is not source or prevalence balance. Descriptive cells below 20 are suppressed in that study.

Every selected query remains in the result denominator, including uncovered truths and empty candidate lists. Shared independent candidate knowledge remains permitted under C2 semantics, and structures can occur as training negatives for other queries. The fixed NP-score prior has unverified upstream overlap. Unknown numeric collision energy is encoded zero without NCE/percent-to-eV conversion. Source annotations are not independently authenticated laboratory labels.

## Reproduction and redistribution

Raw inputs, molecular artifacts, per-query records and weights stay outside Git. Use new immutable artifact destinations, content hashes, verified fingerprints and fresh D2 indexes. Only aggregates, source/protocol inventories and runtime versions are committed. [Data integrity](docs/guides/data-integrity.md) · [Reproduction](docs/guides/reproducibility.md) · [Model limitations](MODEL_CARD.md).
