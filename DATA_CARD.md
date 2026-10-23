# Data card — spectra, structures and derived evidence

## Sources and access

The historical competition training profile records **2,539,608 spectra**, **275,810 standard structure keys** and **11 source libraries**. The public example test contains 1,213 spectra for 400 molecules. These counts come from [the committed profile](reports/data_profile.json), not a new full-source read.

Competition spectra must be acquired through Kaggle after accepting the applicable terms. They are not redistributed in this repository. COCONUT, PubChemLite and optional structure sources must retain their own release identifiers, citations and licensing; the historical notebooks describe the source releases they used. Source code licensing does not grant rights to competition inputs, third-party models or the NEXUS identity.

## Unit of analysis and leakage

The prediction unit is a **molecule**, not a spectrum. Multiple spectra per structure, repeated acquisitions, different libraries and tautomer-equivalent representations require grouping and reference exclusion before scoring. Summing distinct structures across libraries double-counts shared identities.

The archived profile finds 1,213 coarse spectrum-key matches between public examples and train. That is a contamination warning, not verified exact peak/intensity identity. Use `src/casmi/validate/spectrum_overlap.py` on authorized inputs for content-level verification. Public examples are suitable for interface checks; they cannot establish hidden-test generalization.

A tautomer-canonical InChIKey14 groups metric-equivalent structures. It does not enforce Bemis–Murcko scaffold separation, stereochemical identification or unseen chemical-family generalization.

## Distribution and quality risks

Library instrument, ionization mode, adduct, peak count and precursor-error distributions differ. V-in is a small, deliberately common natural-product set; it can conceal decoy costs and fail to represent unknown hidden-test novelty. Pool source composition itself can leak the simulation class into ranker features.

## Rebuild policy

Keep all raw files and generated molecular artifacts outside Git. Record content hashes with `scripts/manifest_artifacts.py`, verify candidate fingerprint provenance, and publish new stores/pools into separate immutable destinations. Revalidation requires the same queries and exclusions across models, negative controls and untouched evaluation data. Cloud warehouses are optional; the public dashboard reads committed aggregates only.
