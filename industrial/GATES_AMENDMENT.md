# GATES AMENDMENT 1 - industrial/ slice, reference-set curation correction
Amended 2026-09-23 18:58 IST, committed standalone BEFORE v2 results. v1 artifacts
preserved in data/v1_superseded/ and results/validation.json (v1) - not re-fished.

## What happened
Boundary diagnosis of pipeline v1 (PC 82/82, GATE1 27/31=0.871 FAIL, GATE2 precision
1.000/AUROC 0.987 PASS) showed 3 of 4 GATE1 misses were NOT ring-hydroxylating
dioxygenase alpha subunits at all:
- P37334, Q46373 = biphenyl dioxygenase subunit BETA (bphE) - structural subunit,
  no Rieske center, outside the locked alpha-subunit scope;
- P77650 (HcaD) = ferredoxin--NAD+ reductase component (EC 1.18.1.3), pulled in by a
  gene-name collision (phdA synonym), outside scope.
Full-set audit found the ec:1.14.12.x queries had swept in 27 non-alpha entries
(26 beta subunits bphE/bnzB/hcaB-class + 1 reductase), including 16 that trained
v1 family profiles. The fourth miss, Q8U671 (DhaA haloalkane dehalogenase,
3-mer jaccard 0.081 to nearest training member), is a GENUINE fold-coverage
boundary miss - same boundary B20 documented: profile-transfer recall is bounded
by training coverage of the characterized canon.

## Correction (curation fix, not outcome tuning)
- Removed all 27 non-alpha entries from the reference set regardless of which side
  of the split they fell on (scope defined in GATES_LOCKED.md: alpha subunits only).
- reference_set.csv v2: 90 enzymes (HLD 42V+12T, RHD 24V+12T). HONEST NOTE: the
  RHD VERIFIED canon is 24, below the >=30/class aspiration; the >=30 total is met
  with THIN tier included (36), reported openly per the B20 convention. The
  alpha-subunit-only characterized canon is genuinely smaller than raw EC sweeps
  suggest - a real curation finding, not padding.
- split.json v2 regenerated with the SAME seed and method (61 train / 29 heldout,
  20 VERIFIED heldout), frozen before v2 mining. THIN audit: no mis-scoped entries.
- HLD set audit: no analogous contamination.

## Reporting
v1 (contaminated set) and v2 (curated set) metrics both reported in the paper;
the v1 GATE1 failure stands on record. Thresholds in v2 are frozen by the same
rules (training CV only, no post-hoc adjustment).
