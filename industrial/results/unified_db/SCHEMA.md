# Unified vetted-candidate database - schema v1 (locked 2026-09-23)
One row per vetted candidate CLUSTER (representative sequence). Cross-class merge of
fleet slices; per-candidate evidence rows, no aggregation loss.
- db_id: UDB-<SLICE>-<cluster id> (stable)
- slice: pesticides | industrial (pet/ rows append when builder 19 outputs land in repo -
  as of 2026-09-23 19:00 IST pet/ is absent from the repo; noted, not padded)
- candidate_class: OP | PYR | HLD | RHD
- family_profile: calling profile HMM id
- rep_id / uniprot_acc / source / length: representative sequence identity + origin corpus
- score / neg_score / margin: frozen-model discriminator scores (margin = score - neg_score)
- novelty_jaccard_vs_canon: max 3-mer Jaccard vs VERIFIED canon (GATE3 proxy, <0.5 ~ novel)
- motif_check: class-specific catalytic motif presence (gxsxg for esterases; rieske CXXH..CXXH for RHD)
- gate3_verdict: NOVEL | HOMOLOG-OF-CANON
- cluster_size, annotation, evidence_refs: cluster cardinality, best annotation, source table path
Checksums: MANIFEST.sha256 (slice level) covers this file.
