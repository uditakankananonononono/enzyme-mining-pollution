# GATES LOCKED - pesticides/ slice (builder 20), project enzyme-mining-pollution
Locked 2026-09-23 13:35 IST, BEFORE any mining/results. Committed standalone per fleet rule.

Scope: organophosphate (OP) hydrolases + pyrethroid (PYR) hydrolases mined from public
metagenomes (MGnify proteins, UniMES) + reference proteomes. Amended bar (shared with B19):
mining platform + validated specificity gate, NOT another curation DB.

## Reference set (assembled, byte-locked 2026-09-23; see data/reference_set.csv)
- 39 OP entries (28 VERIFIED / 11 THIN), 37 PYR entries (19 VERIFIED / 18 THIN) = 76 total.
- Tiers: VERIFIED = Swiss-Prot reviewed OR record carries FUNCTION/catalytic/referenced
  experimental evidence for the target activity. THIN = annotation-only propagation.
- AMENDMENT to the wording reported at milestone 1: the >=30/class target is met with
  THIN-tier included; the pyrethroid characterized canon in the literature is genuinely
  small (~19 verifiable). Recall/specificity gates are measured on VERIFIED-tier held-out
  members only; THIN members are tracked and reported separately, never silently promoted.
- Split BEFORE mining: seeded (seed=20260923) 70/30 train/held-out, stratified by
  (class, family), frozen in data/split.json. Held-out sequences never used for profiles.

## GATE 1 - Reference recall >=90% (parent-set)
Pipeline (family profiles + active-site motif filters) must recover >=90% of HELD-OUT
VERIFIED enzymes hidden in the mining corpus (MGnify/UniMES hits + decoy proteomes).
Positive-control sanity first: 100% recovery of TRAINING VERIFIED enzymes, unblinded.
Operating threshold fixed by training-set cross-validation ONLY, before held-out test.
No post-hoc tuning; a failed gate is reported as a failure.

## GATE 2 - Specificity discriminator vs documented false-positive modes (parent-set)
Labeled negatives byte-locked in data/negative_set.csv:
  (a) 30 PTE-superfamily neighbors incl. PLLs (SacPox Q97VT7, SsoPox B5BLW5), E. coli PHP
      P45548, insect phosphotriesterase-related proteins, N-acetyltaurine hydrolases,
      dihydroorotases, adenosine deaminases (documented: most PTE-like hits are native
      lactonases - Afriat 2006; Hawwa 2009);
  (b/c) 31 generic bacterial carboxylesterases/lipases without pesticide evidence
      (documented: functional screens yield 6-12 esterase+ clones per 1 true pyrethroid
      hydrolase - Sys410/Pye3 screens; pNP-ester activity does NOT imply pyrethroid activity).
Discriminator emits CONFIDENT / PROMISCUOUS-TIER / REJECT. At the frozen threshold:
precision >=90% on held-out labeled VERIFIED-positives-vs-negatives AND AUROC >=0.9.
PLLs scoring CONFIDENT count as false positives (they may land PROMISCUOUS-TIER).

## GATE 3 - Novelty
Claimed novel candidates: <80% pairwise identity to any VERIFIED reference enzyme AND
absent from BRENDA/literature curation for the claimed activity (checked, not assumed).

## Reporting discipline
Honest verified/thin/missing counts per family in every report. Missing canon noted
(Est804, Cest2923, EstSt7, PytY names did not resolve to accessions in UniProt/NCBI
protein as of 2026-09-23 - reported, not padded). Failures reported as-is; no re-fishing.

## Corpus & compute
MGnify sequence-search API + UniProt/UniMES API server-side searches (2-core/2GB-safe).
Structure work: colabfold-scale on selected top candidates only.
