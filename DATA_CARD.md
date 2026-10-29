# Data card — spectra, structures and derived evidence

## Sources and access

The historical competition training profile records **2,539,608 spectra**, **275,810 standard structure keys** and **11 source libraries**. The public example test contains 1,213 spectra for 400 molecules. These counts come from [the committed profile](reports/data_profile.json), not a new full-source read.

Competition spectra must be acquired through Kaggle after accepting the applicable terms. They are not redistributed in this repository. COCONUT, PubChemLite and optional structure sources must retain their own release identifiers, citations and licensing; the historical notebooks describe the source releases they used. Source code licensing does not grant rights to competition inputs, third-party models or the NEXUS identity.

## Unit of analysis and leakage

The prediction unit is a **molecule**, not a spectrum. Multiple spectra per structure, repeated acquisitions, different libraries and tautomer-equivalent representations require grouping and reference exclusion before scoring. Summing distinct structures across libraries double-counts shared identities.

The [executed Kaggle content audit](docs/studies/public-example-overlap.md) scanned all 2,539,608 training rows and confirmed exact precursor, full ordered m/z and normalized intensity matches for all 1,213 public spectra (400 molecules), also consistent in adduct and ionization mode. Two coarse-bucket candidate pairs were false matches. This supersedes the archived profile's approximate overlap evidence for the hashed input snapshot. Public examples are suitable for interface checks; they cannot establish hidden-test generalization. Stored numerical identity does not independently authenticate acquisition provenance.

A tautomer-canonical InChIKey14 groups metric-equivalent structures. It does not enforce Bemis–Murcko scaffold separation, stereochemical identification or unseen chemical-family generalization.

## Distribution and quality risks

The [completed source rebuild and model study](reports/revalidation_review.md) identified 1,710 shared identities with conflicting neutral masses (3,891 provider rows). Its explicit whole-identity quarantine also excluded 13 unparseable selected candidates, leaving **1,071,863** fingerprint-verified candidate structures. All **250** evaluation molecules remained in scoring; one standard truth key was absent from the retained pool. Labels and reference exclusions were derived independently of candidate retention. These checks establish computational consistency, not independent authentication of provider annotations.

Library instrument, ionization mode, adduct, peak count and precursor-error distributions differ. V-in is a small, deliberately common natural-product set; it can conceal decoy costs and fail to represent unknown hidden-test novelty. Pool source composition itself can leak the simulation class into ranker features.

## Fresh-study population and information boundaries

The [executed internal benchmark](reports/advanced_review.md) selects one representative per generic Murcko/acyclic connectivity family: 600 training, 150 development and 250 test molecules, excluding historical V-in families. Selected metric identities and families do not cross supervised splits. This reuses existing source libraries; it is not newly acquired external data. The 250 test queries include 141 pluskal_ms2 and 56 GNPS molecules, so family balancing does not imply source balancing. Cells below 20 are suppressed in descriptive slices.

Candidate/reference structures remain task-defined and can occur as negatives for other training queries. All query outcomes stay in denominators: C1/C2 pool truth coverage is 246/250; C3 truth coverage after generation is only 26/250. The fixed NP-score prior has unverified upstream training overlap. These boundaries prevent interpreting the result as universal unseen-chemistry performance. Only aggregates, hashes and version inventories are committed; private row details and weights remain scratch artifacts.

## Reproduction and redistribution

The [external-source study](docs/studies/external-source-validation.md) prepares separately downloaded MassBank CC0 and directly contributed GNPS standard spectra. Original record hashes, source URLs and available attribution are retained. MassBank uses explicit per-record CC0 licenses; GNPS uses its direct-contribution CC0 policy and excludes imported/inferred libraries. These legal/source distinctions are preserved in manifests. The first two audited populations leave zero and seven novel molecules; neither qualifies as the declared external validation population. Unknown collision energies are not converted from NCE/percent to eV.

Keep all raw files and generated molecular artifacts outside Git. Record content hashes with `scripts/manifest_artifacts.py`, verify candidate fingerprint provenance, and publish new stores/pools into separate immutable destinations. Revalidation requires the same queries and exclusions across models, negative controls and untouched evaluation data. Cloud warehouses are optional; the public dashboard reads committed aggregates only.

## Completed family-purged study

The authenticated family-purged study removed 261,826 of 2,539,608 training spectra and fitted eight fresh learned artifacts. On 150 internal test molecules, expanded C3 truth coverage rose from 13 to 17 (paired change +0.0267, 95% CI [+0.0067, +0.0533]). On 150 external molecules, expanded C2/C3 MRR@25 was 0.4454/0.0466; expanded-minus-legacy ranking effects remain unresolved in the joint 11-contrast intervals. External C3 truth coverage remains only 23/150. Original-snapshot novelty gates remain unsuccessful; this is generalization relative to a rebuilt purged corpus, not independent laboratory replication. [Methods, all outcomes and limitations](reports/purged_external_review.md).

The explicit CASMI_DATA snapshot takes precedence over automatic Kaggle discovery. All reference artifacts are rebuilt from filtered rows before fitting. Undefined generic-family labels quarantine the entire affected declared identity under the outcome-blind amendment; the failed first attempt remains preserved.
